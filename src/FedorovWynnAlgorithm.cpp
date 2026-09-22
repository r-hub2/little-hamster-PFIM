// Fedorov-Wynn exchange for D-optimal weights on discrete protocols (Fedorov 1972).
//
// Outer loop (each cycle):
//   1. Fedorov augmentation — add outsiders with Tr(F^{-1} F_i) > p(1+δ)
//   2. Wynn weight exchange on the current support (always, even if no add)
//   3. Prune near-zero weights
//   4. Kiefer–Wolfowitz check: max_i Tr(F^{-1} F_i) ≤ p(1+δ) over the full grid
//
// Directional gain for candidate i is Tr(F^{-1} F_i). The Fedorov insert weight
// α* = (φ − p) / (p (φ − 1)) equals 1 when p = 1 (complete replacement) — that
// case is handled explicitly (not rejected by a share < 1 guard).
//
// Storage: each elementary FIM is packed upper-triangle (p(p+1)/2). Carathéodory
// bound: active support must not exceed packed size; capacity is packed+1 so one
// Fedorov insert can temporarily overshoot before the prune phase.
//
// R↔C++: initial support from protdep / freqdep (1-based grid indices + weights);
// each seed must match a full sampling-grid row. Returned fisher is packed, not dense.
// Optional delta (default 1e-6) mirrors MultiplicativeAlgorithm's relative gap.
//
// Population / individual FIMs only — Bayesian prior is applied on the R side,
// not inside this kernel.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <pfim/pfim-linalg.hpp>
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace fw {

// --- Tolerances and iteration budgets ---
constexpr double k_ftol              = 1e-6;
constexpr double k_default_delta     = 1e-4;  // optimality gap: max φ ≤ p(1+δ)
constexpr double k_weight_floor      = 1e-8;
constexpr double k_line_search_floor = 1e-4;
constexpr double k_singular_det      = 1e30;
constexpr double k_cold_start_det    = -1e10;
constexpr double k_share_one_tol     = 1e-12; // treat α* ≥ 1 − tol as full replacement
constexpr int    k_max_wynn_iters    = 500;
constexpr int    k_max_outer_iters   = 100;
constexpr int    k_max_add_attempts  = 100;

// Solver outcome reported back to R via status string.
enum class SolveStatus {
  success,              // Kiefer–Wolfowitz: max φ ≤ p(1+δ) on the full grid
  incomplete,           // outer cycle budget exhausted; weights usable but not certified
  support_overflow,
  augmentation_failed,
  weight_optimization_failed,
  singular_fim
};

// Inner Wynn weight-exchange loop outcome.
enum class WynnStatus {
  converged,
  iteration_limit,
  singular_fim
};

// Result of shrinking an exchange step when the full step does not improve D.
enum class LineSearchStatus {
  improved,
  tolerance_reached,
  rejected
};

// --- PackedFim: upper-triangle storage for symmetric FIMs ---
// Stores a symmetric p×p matrix as p(p+1)/2 entries in column-major packed order
// (index k = row*(row+1)/2 + col for row >= col). Same layout R writes into
// protocols[[6]] / fisher buffers — do not reinterpret as row-major or dense.
class PackedFim {
public:
  explicit PackedFim( int param_dim = 0 )
    : param_dim_( param_dim ),
      values_( static_cast<arma::uword>( packed_size( param_dim ) ), arma::fill::zeros ) {}

  static int packed_size( int param_dim ) {
    return param_dim * ( param_dim + 1 ) / 2;
  }

  int param_dim() const { return param_dim_; }
  const arma::vec& values() const { return values_; }
  arma::vec& values() { return values_; }
  void fill_zero() { values_.zeros(); }

  // Unpack packed → dense symmetric matrix (both triangles filled).
  static void expand_into( const arma::vec& packed, int param_dim, arma::mat& dense ) {
    arma::uword k = 0;
    for ( int row = 0; row < param_dim; ++row )
      for ( int col = 0; col <= row; ++col )
        dense( row, col ) = dense( col, row ) = packed( k++ );
  }

  static arma::mat expand( const arma::vec& packed, int param_dim ) {
    arma::mat dense( param_dim, param_dim, arma::fill::zeros );
    expand_into( packed, param_dim, dense );
    return dense;
  }

  // Pack dense lower triangle → packed vector (upper mirrored by construction).
  static void compress_into( const arma::mat& dense, int param_dim, arma::vec& packed ) {
    for ( int row = 0; row < param_dim; ++row )
      for ( int col = 0; col <= row; ++col )
        packed( static_cast<arma::uword>( row * ( row + 1 ) / 2 + col ) ) = dense( row, col );
  }

  static arma::vec compress( const arma::mat& dense, int param_dim ) {
    arma::vec packed( packed_size( param_dim ), arma::fill::zeros );
    compress_into( dense, param_dim, packed );
    return packed;
  }

  // Directional gain Tr(F^{-1} F_i) without forming dense products:
  // off-diagonals counted twice (symmetry), diagonals once.
  static double trace_with_inverse( const arma::vec& fim,
                                    const arma::vec& inv_fim,
                                    int param_dim ) {
    double trace = 0.0;
    for ( int row = 0; row < param_dim; ++row ) {
      for ( int col = 0; col < row; ++col ) {
        const arma::uword k = static_cast<arma::uword>( row * ( row + 1 ) / 2 + col );
        trace += 2.0 * fim( k ) * inv_fim( k );
      }
      const arma::uword k = static_cast<arma::uword>( row * ( row + 3 ) / 2 );
      trace += fim( k ) * inv_fim( k );
    }
    return trace;
  }

private:
  int param_dim_;
  arma::vec values_;
};

// One discrete protocol: sampling times plus its packed elementary FIM F_i.
struct CandidateProtocol {
  std::vector<double> sample_times;
  PackedFim information;
};

// --- ProtocolGrid: discrete candidate set for Fedorov augmentation ---
// Holds every eligible protocol; directional_gain evaluates Tr(F^{-1} F_i).
class ProtocolGrid {
public:
  static ProtocolGrid from_matrices( const NumericMatrix& sample_times,
                                     const NumericMatrix& information_rows,
                                     int param_dim,
                                     int n_sample_times,
                                     double information_scale ) {
    const int n_candidates = sample_times.nrow();
    if ( information_rows.nrow() != n_candidates ||
         information_rows.ncol() != PackedFim::packed_size( param_dim ) )
      stop( "Fedorov-Wynn: information matrix shape does not match candidates." );
    if ( sample_times.ncol() != n_sample_times )
      stop( "Fedorov-Wynn: sample time columns mismatch." );

    ProtocolGrid grid;
    grid.candidates_.resize( static_cast<size_t>( n_candidates ) );
    for ( int row = 0; row < n_candidates; ++row ) {
      CandidateProtocol& candidate = grid.candidates_[ static_cast<size_t>( row ) ];
      candidate.information = PackedFim( param_dim );
      candidate.sample_times.assign( n_sample_times, 0.0 );
      for ( int t = 0; t < n_sample_times; ++t )
        candidate.sample_times[ static_cast<size_t>( t ) ] = sample_times( row, t );
      for ( int k = 0; k < PackedFim::packed_size( param_dim ); ++k )
        candidate.information.values()( k ) = information_rows( row, k ) * information_scale;
    }
    return grid;
  }

  int count() const { return static_cast<int>( candidates_.size() ); }
  const CandidateProtocol& candidate( int index ) const {
    return candidates_.at( static_cast<size_t>( index ) );
  }

  // Sensitivity / directional derivative: φ(i, ξ) = Tr(F(ξ)^{-1} F_i).
  // At a D-optimum, max_i φ(i, ξ) = p (= param_dim).
  double directional_gain( int candidate_index, const arma::vec& inv_information ) const {
    const CandidateProtocol& item = candidate( candidate_index );
    return PackedFim::trace_with_inverse(
      item.information.values(), inv_information, item.information.param_dim()
    );
  }

  // max_i φ(i, ξ) over the full candidate grid (Kiefer–Wolfowitz check).
  double max_directional_gain( const arma::vec& inv_information ) const {
    double best = -std::numeric_limits<double>::infinity();
    for ( int i = 0; i < count(); ++i )
      best = std::max( best, directional_gain( i, inv_information ) );
    return best;
  }

private:
  std::vector<CandidateProtocol> candidates_;
};

// Protocol currently in the design support (grid index, weight, copy of F_i).
struct ActiveProtocol {
  int candidate_index = 0;
  double weight = 0.0;
  std::vector<double> sample_times;
  PackedFim information;
};

// Projected gradient of log-det on the probability simplex (free weights only).
struct ProjectedGradient {
  std::vector<double> components;
  std::optional<int> descent_index;  // most negative component → mass to remove
  double max_abs_component = 0.0;
};

// Feasible exchange ray: w ← w + α d, α limited so one weight hits zero.
struct WynnStep {
  std::vector<double> direction;
  double step_length = 0.0;
  int entering_index = 0;
  int leaving_index  = 0;
};

// Aggregate result of one Wynn weight-optimisation run.
struct WynnOutcome {
  WynnStatus status = WynnStatus::converged;
  int inactive_weight_count = 0;
};

// --- WeightedDesign: current support, mixture FIM, Fedorov add + Wynn weights ---
// Maintains ξ = Σ w_i δ_i with F(ξ) = Σ w_i F_i, D = exp(log det(F)/p).
class WeightedDesign {
public:
  WeightedDesign( int param_dim, int max_support, double delta = k_default_delta )
    : param_dim_( param_dim ),
      max_support_( max_support ),
      active_size_( 0 ),
      delta_( delta > 0.0 && R_finite( delta ) ? delta : k_default_delta ),
      mixture_fim_( param_dim ),
      inv_mixture_fim_( param_dim ),
      determinant_( 1.0 ),
      fim_ready_( false ),
      singular_fim_( false ),
      active_( static_cast<size_t>( max_support ) ),
      scratch_weights_( static_cast<size_t>( max_support ), 0.0 ),
      dense_fim_( static_cast<arma::uword>( param_dim ),
                  static_cast<arma::uword>( param_dim ), arma::fill::zeros ),
      dense_inv_( static_cast<arma::uword>( param_dim ),
                  static_cast<arma::uword>( param_dim ), arma::fill::zeros ) {}

  int param_dim() const { return param_dim_; }
  int active_size() const { return active_size_; }
  double delta() const { return delta_; }
  bool is_singular() const { return singular_fim_; }
  const arma::vec& mixture_fim() const { return mixture_fim_.values(); }
  const arma::vec& inv_mixture_fim() const { return inv_mixture_fim_.values(); }

  const ActiveProtocol& active_protocol( int index ) const {
    return active_.at( static_cast<size_t>( index ) );
  }

  void clear_support() {
    active_size_ = 0;
    touch_fim();
  }

  bool append( int candidate_index, double weight, const CandidateProtocol& source ) {
    if ( active_size_ >= max_support_ )
      return false;
    ActiveProtocol& slot = active_.at( static_cast<size_t>( active_size_ ) );
    slot.candidate_index = candidate_index;
    slot.weight          = weight;
    slot.sample_times    = source.sample_times;
    slot.information     = source.information;
    ++active_size_;
    touch_fim();
    return true;
  }

  bool includes( int candidate_index ) const {
    for ( int i = 0; i < active_size_; ++i )
      if ( active_.at( static_cast<size_t>( i ) ).candidate_index == candidate_index )
        return true;
    return false;
  }

  // Kiefer–Wolfowitz: design is D-optimal iff max_i φ(i, ξ) ≤ p(1+δ).
  bool is_d_optimal( const ProtocolGrid& grid ) {
    update_information();
    if ( is_singular() || active_size_ <= 0 )
      return false;
    const double threshold =
      static_cast<double>( param_dim_ ) * ( 1.0 + delta_ );
    return grid.max_directional_gain( inv_mixture_fim_.values() ) <= threshold;
  }

  int count_inactive_weights() const {
    int inactive = 0;
    for ( int i = 0; i < active_size_; ++i )
      if ( weight( i ) < k_weight_floor )
        ++inactive;
    return inactive;
  }

  void prune_inactive_support() {
    // Compact active_ by dropping protocols with weight below k_weight_floor.
    int skipped = 0;
    int kept = active_size_;
    for ( int i = 0; i < active_size_; ++i ) {
      if ( weight( i ) < k_weight_floor ) {
        ++skipped;
        --kept;
      } else if ( skipped > 0 ) {
        active_.at( static_cast<size_t>( i - skipped ) ) = active_.at( static_cast<size_t>( i ) );
      }
    }
    active_size_ = kept;
    touch_fim();
  }

  // --- Mixture FIM and D-criterion ---
  // Rebuild F(ξ) = Σ w_i F_i, invert via Cholesky, store packed F^{-1}.
  // D-criterion used for comparisons: D = exp(log_det / p).
  double update_information() {
    if ( fim_ready_ )
      return determinant_;

    mixture_fim_.fill_zero();
    for ( int i = 0; i < active_size_; ++i )
      mixture_fim_.values() += weight( i ) * active_.at( static_cast<size_t>( i ) ).information.values();

    PackedFim::expand_into( mixture_fim_.values(), param_dim_, dense_fim_ );
    double log_det = 0.0;
    if ( !pfim::try_chol_inv_logdet( dense_fim_, dense_inv_, log_det ) ) {
      singular_fim_ = true;
      determinant_  = 0.0;
      inv_mixture_fim_.fill_zero();
      fim_ready_      = true;
      return determinant_;
    }

    singular_fim_             = false;
    // Store D-criterion exp(log_det/p) (same scale as MultiplicativeAlgorithm).
    // Avoids Inf overflow of raw det(F) while preserving monotonic comparisons.
    determinant_              = std::exp( log_det / static_cast<double>( param_dim_ ) );
    PackedFim::compress_into( dense_inv_, param_dim_, inv_mixture_fim_.values() );
    fim_ready_                = true;
    return determinant_;
  }

  // --- Fedorov augmentation step ---
  // Among candidates not already in support, pick max Tr(F^{-1} F_i).
  // Improving if φ > p(1+δ); insert with Fedorov weight share and rescale others.
  // When p = 1, α* = 1 (complete replacement) — see Fedorov (1972) eq. (2.3).
  bool try_add_best_candidate( const ProtocolGrid& grid ) {
    update_information();
    if ( is_singular() )
      return false;

    const double target    = static_cast<double>( param_dim_ );
    const double threshold = target * ( 1.0 + delta_ );
    double best_gain = -std::numeric_limits<double>::infinity();
    std::optional<int> best_candidate;

    for ( int candidate = 0; candidate < grid.count(); ++candidate ) {
      if ( includes( candidate ) )
        continue;
      const double gain = grid.directional_gain( candidate, inv_mixture_fim_.values() );
      if ( gain > best_gain ) {
        best_gain = gain;
        best_candidate = candidate;
      }
    }

    // Stationary for augmentation: no outsider with φ > p(1+δ).
    if ( !best_candidate.has_value() || !( best_gain > threshold ) )
      return false;

    // α* = (φ − p) / (p (φ − 1)). Denominator vanishes at φ = 1; for an improving
    // point φ > p ≥ 1 so φ = 1 only when p = 1 and not improving — already rejected.
    if ( std::abs( best_gain - 1.0 ) < 1e-14 )
      return false;
    const double share = ( best_gain - target ) / ( target * ( best_gain - 1.0 ) );
    if ( !R_finite( share ) || share <= 0.0 )
      return false;

    // α* = 1 exactly when p = 1; also clamp tiny float overshoot above 1.
    if ( share >= 1.0 - k_share_one_tol ) {
      clear_support();
      if ( !append( *best_candidate, 1.0, grid.candidate( *best_candidate ) ) )
        return false;
      update_information();
      return true;
    }

    if ( !append( *best_candidate, share, grid.candidate( *best_candidate ) ) )
      return false;
    // Rescale existing weights so Σ w = 1 after inserting share.
    for ( int i = 0; i < active_size_ - 1; ++i )
      active_.at( static_cast<size_t>( i ) ).weight *= ( 1.0 - share );
    touch_fim();
    update_information();
    return true;
  }

  // --- Wynn weight optimisation on fixed support ---
  // Projected-gradient exchange: take full feasible step, else shrink via line search;
  // if still stuck, temporarily pin a near-zero weight that remains improving.
  WynnOutcome optimize_weights( const ProtocolGrid& grid ) {
    int inactive_weights = 0;
    std::optional<int> pinned_index;
    double criterion = update_information();
    if ( is_singular() )
      return WynnOutcome{ WynnStatus::singular_fim, inactive_weights };

    // Wynn exchange: projected gradient on the simplex, line search, occasional support prune.
    for ( int iteration = 0; iteration < k_max_wynn_iters; ++iteration ) {
      Rcpp::checkUserInterrupt();
      save_weights( scratch_weights_ );
      const double previous_criterion = criterion;

      if ( !pinned_index.has_value() )
        inactive_weights = count_inactive_weights();

      ProjectedGradient gradient = compute_projected_gradient(
        grid, pinned_index, inactive_weights
      );
      // No descent direction → stationary on current free weights.
      if ( !gradient.descent_index.has_value() )
        return WynnOutcome{ WynnStatus::converged, inactive_weights };
      if ( gradient.components[ static_cast<size_t>( *gradient.descent_index ) ] == 0.0 )
        return WynnOutcome{ WynnStatus::converged, inactive_weights };
      if ( pinned_index.has_value() &&
           gradient.components[ static_cast<size_t>( *pinned_index ) ] <= 0.0 )
        return WynnOutcome{ WynnStatus::converged, inactive_weights };

      pinned_index.reset();
      const WynnStep step = build_exchange_step( gradient );
      apply_exchange_step( step );
      const double trial_criterion = update_information();

      if ( trial_criterion > criterion ) {
        // Full exchange improved D — accept and refresh inactive count.
        save_weights( scratch_weights_ );
        criterion = trial_criterion;
        inactive_weights = count_inactive_weights();
      } else {
        // Full step failed — restore and try halved steps along the same ray.
        restore_weights( scratch_weights_ );
        const LineSearchStatus search = line_search_along_step(
          step, scratch_weights_, criterion, previous_criterion, gradient.max_abs_component
        );
        if ( search == LineSearchStatus::rejected ) {
          // Pin a removable near-zero support point and continue.
          pinned_index = find_removable_support( grid, inactive_weights );
          if ( !pinned_index.has_value() )
            return WynnOutcome{ WynnStatus::converged, inactive_weights };
        }
      }
    }
    return WynnOutcome{ WynnStatus::iteration_limit, inactive_weights };
  }

private:
  double weight( int index ) const {
    return active_.at( static_cast<size_t>( index ) ).weight;
  }

  ActiveProtocol& mutable_protocol( int index ) {
    return active_.at( static_cast<size_t>( index ) );
  }

  void touch_fim() { fim_ready_ = false; singular_fim_ = false; }

  void save_weights( std::vector<double>& buffer ) const {
    for ( int i = 0; i < active_size_; ++i )
      buffer[ static_cast<size_t>( i ) ] = weight( i );
  }

  void restore_weights( const std::vector<double>& buffer ) {
    for ( int i = 0; i < active_size_; ++i )
      mutable_protocol( i ).weight = buffer[ static_cast<size_t>( i ) ];
    touch_fim();
  }

  // Treat near-zero weights as inactive unless temporarily pinned for gradient.
  bool is_active_weight( int index, std::optional<int> pinned_index ) const {
    return weight( index ) >= k_weight_floor ||
      ( pinned_index.has_value() && *pinned_index == index );
  }

  // --- Projected gradient on the simplex ---
  // g_i = Tr(F^{-1} F_i) for free weights; center by mean so Σ g = 0 on free set.
  // Most negative g_i marks the protocol from which mass should leave.
  ProjectedGradient compute_projected_gradient( const ProtocolGrid& grid,
                                                std::optional<int> pinned_index,
                                                int inactive_weights ) {
    update_information();
    ProjectedGradient result;
    result.components.assign( static_cast<size_t>( active_size_ ), 0.0 );

    double sum = 0.0;
    for ( int i = 0; i < active_size_; ++i ) {
      if ( !is_active_weight( i, pinned_index ) )
        continue;
      result.components[ static_cast<size_t>( i ) ] = grid.directional_gain(
        active_.at( static_cast<size_t>( i ) ).candidate_index,
        inv_mixture_fim_.values()
      );
      sum += result.components[ static_cast<size_t>( i ) ];
    }

    const int free_weights = active_size_ - inactive_weights;
    if ( free_weights <= 0 )
      return result;

    // Centering: g ← g − mean(g) enforces tangent-to-simplex (Σ Δw = 0).
    const double center = sum / static_cast<double>( free_weights );
    for ( int i = 0; i < active_size_; ++i ) {
      if ( !is_active_weight( i, pinned_index ) )
        continue;
      result.components[ static_cast<size_t>( i ) ] -= center;
      const double value = result.components[ static_cast<size_t>( i ) ];
      result.max_abs_component = std::max( result.max_abs_component, std::abs( value ) );
      if ( value < 0.0 &&
           ( !result.descent_index.has_value() ||
             value < result.components[ static_cast<size_t>( *result.descent_index ) ] ) )
        result.descent_index = i;
    }
    return result;
  }

  // --- Exchange step (vertex-to-face ray) ---
  // Maximal α ≥ 0 with w + α d ≥ 0; leaving_index is the weight that hits zero first.
  WynnStep build_exchange_step( const ProjectedGradient& gradient ) const {
    const int entering = *gradient.descent_index;
    WynnStep step;
    step.direction     = gradient.components;
    step.entering_index = entering;
    step.leaving_index  = entering;
    step.step_length    = -weight( entering ) / step.direction[ static_cast<size_t>( entering ) ];

    for ( int i = 0; i < active_size_; ++i ) {
      if ( step.direction[ static_cast<size_t>( i ) ] >= 0.0 )
        continue;
      const double bound = -weight( i ) / step.direction[ static_cast<size_t>( i ) ];
      if ( bound < step.step_length ) {
        step.step_length   = bound;
        step.leaving_index = i;
      }
    }
    return step;
  }

  void apply_exchange_step( const WynnStep& step ) {
    for ( int i = 0; i < active_size_; ++i )
      mutable_protocol( i ).weight += step.step_length * step.direction[ static_cast<size_t>( i ) ];
    mutable_protocol( step.leaving_index ).weight = 0.0;
    touch_fim();
  }

  // --- Backtracking line search along the exchange ray ---
  // Halve α until D improves or α ‖d‖_∞ falls below k_line_search_floor.
  LineSearchStatus line_search_along_step( const WynnStep& step,
                                           std::vector<double>& saved_weights,
                                           double& criterion,
                                           double previous_criterion,
                                           double max_direction ) {
    double shrink = 1.0;
    do {
      shrink *= 0.5;
      for ( int i = 0; i < active_size_; ++i )
        mutable_protocol( i ).weight =
          saved_weights[ static_cast<size_t>( i ) ] +
          shrink * step.step_length * step.direction[ static_cast<size_t>( i ) ];
      touch_fim();
      const double trial = update_information();
      if ( trial <= criterion )
        continue;

      save_weights( saved_weights );
      criterion = trial;
      if ( !R_finite( previous_criterion ) || previous_criterion <= 0.0 )
        return LineSearchStatus::rejected;
      // Relative change in D vs previous outer iterate: below k_ftol → done enough.
      return std::abs( ( previous_criterion - criterion ) / previous_criterion ) > k_ftol
        ? LineSearchStatus::improved
        : LineSearchStatus::tolerance_reached;
    } while ( step.step_length * shrink * max_direction > k_line_search_floor );

    restore_weights( saved_weights );
    return LineSearchStatus::rejected;
  }

  // Near-zero weight still with gain > p(1+δ) may be pinned so gradient can revive it.
  std::optional<int> find_removable_support( const ProtocolGrid& grid,
                                             int& inactive_weights ) {
    update_information();
    const double threshold =
      static_cast<double>( param_dim_ ) * ( 1.0 + delta_ );
    for ( int i = 0; i < active_size_; ++i ) {
      if ( weight( i ) > k_weight_floor )
        continue;
      const double gain = grid.directional_gain(
        active_.at( static_cast<size_t>( i ) ).candidate_index,
        inv_mixture_fim_.values()
      );
      if ( gain > threshold ) {
        inactive_weights -= 1;
        return i;
      }
    }
    return std::nullopt;
  }

  int param_dim_;
  int max_support_;
  int active_size_;
  double delta_;
  PackedFim mixture_fim_;
  PackedFim inv_mixture_fim_;
  double determinant_;
  bool fim_ready_;
  bool singular_fim_;
  std::vector<ActiveProtocol> active_;
  std::vector<double> scratch_weights_;
  arma::mat dense_fim_;
  arma::mat dense_inv_;
};

// FedorovWynnSolver: Fedorov add, Wynn weights, prune, then optimality check.
// Alternates augmentation and weight exchange until Kiefer–Wolfowitz holds or
// the outer budget is exhausted. Wynn runs even when no outsider was added so
// poor initial weights on a correct support are still rebalanced.
class FedorovWynnSolver {
public:
  FedorovWynnSolver( ProtocolGrid grid, WeightedDesign design, bool log_progress )
    : grid_( std::move( grid ) ),
      design_( std::move( design ) ),
      log_progress_( log_progress ) {}

  SolveStatus solve() {
    double determinant = design_.update_information();
    if ( design_.is_singular() )
      return SolveStatus::singular_fim;

    for ( int cycle = 0; cycle < k_max_outer_iters; ++cycle ) {
      Rcpp::checkUserInterrupt();

      // Fedorov: add outsiders with φ > p(1+δ) until none remain.
      if ( !augment_until_stationary( determinant ) )
        return design_.is_singular() ? SolveStatus::singular_fim
                                     : SolveStatus::augmentation_failed;

      // Wynn: re-optimise weights on the current support (even if no add).
      const WynnOutcome wynn = design_.optimize_weights( grid_ );
      if ( wynn.status == WynnStatus::singular_fim )
        return SolveStatus::singular_fim;
      if ( wynn.status == WynnStatus::iteration_limit )
        return SolveStatus::weight_optimization_failed;

      // Drop near-zero weights; refuse if support exceeds Carathéodory packed
      // size (capacity was packed+1 to allow one temporary overshoot).
      design_.prune_inactive_support();
      if ( design_.active_size() > PackedFim::packed_size( design_.param_dim() ) )
        return SolveStatus::support_overflow;

      determinant = design_.update_information();
      if ( design_.is_singular() )
        return SolveStatus::singular_fim;
      if ( log_progress_ )
        Rcout << "Fedorov-Wynn cycle " << ( cycle + 1 )
              << ": D = " << determinant
              << ", max_phi = "
              << grid_.max_directional_gain( design_.inv_mixture_fim() )
              << "\n";

      // Optimality: max_i φ(i, ξ) ≤ p(1+δ) on the full grid.
      if ( design_.is_d_optimal( grid_ ) )
        return SolveStatus::success;
    }
    if ( design_.is_singular() )
      return SolveStatus::singular_fim;
    return SolveStatus::incomplete;
  }

  const WeightedDesign& design() const { return design_; }

private:
  // Add outsiders with φ > p(1+δ) until none remain (or singular retries expire).
  bool augment_until_stationary( double& determinant ) {
    int attempts = 0;
    while ( attempts < k_max_add_attempts ) {
      Rcpp::checkUserInterrupt();
      if ( !design_.try_add_best_candidate( grid_ ) )
        return !design_.is_singular();
      ++attempts;
      determinant = design_.update_information();
      // Legacy cold-start: keep trying while the mixture is singular / sentinel.
      if ( !( determinant < k_cold_start_det || design_.is_singular() ) )
        continue;  // healthy add — keep scanning for further outsiders
    }
    return !( design_.is_singular() || determinant < k_cold_start_det );
  }

  ProtocolGrid grid_;
  WeightedDesign design_;
  bool log_progress_;
};

// Build initial support from 1-based R indices (protdep) and weights (freqdep).
WeightedDesign initial_design( const ProtocolGrid& grid,
                               const IntegerVector& initial_indices,
                               const NumericVector& initial_weights,
                               double delta ) {
  if ( initial_indices.size() < 2 )
    stop( "Fedorov-Wynn: initial protocol indices are incomplete." );

  const int start_size = initial_indices[0];
  const int param_dim  = grid.count() > 0 ? grid.candidate( 0 ).information.param_dim() : 0;
  // Carathéodory: at most p(p+1)/2 points needed for a D-optimal mixture;
  // +1 leaves room for one Fedorov insert before prune enforces the bound.
  const int max_support = PackedFim::packed_size( param_dim ) + 1;

  if ( start_size <= 0 || param_dim <= 0 )
    stop( "Fedorov-Wynn: invalid initial support size." );
  if ( initial_indices.size() < start_size + 1 )
    stop( "Fedorov-Wynn: initial protocol indices do not match support size." );
  if ( initial_weights.size() < start_size )
    stop( "Fedorov-Wynn: initial weights shorter than support size." );

  WeightedDesign design( param_dim, max_support, delta );
  for ( int slot = 0; slot < start_size; ++slot ) {
    const int candidate = initial_indices[ slot + 1 ] - 1;
    if ( candidate < 0 || candidate >= grid.count() )
      stop( "Fedorov-Wynn: candidate index out of range." );
    if ( design.includes( candidate ) )
      stop( "Fedorov-Wynn: duplicate initial candidate index." );
    if ( !design.append( candidate, initial_weights[ slot ], grid.candidate( candidate ) ) )
      stop( "Fedorov-Wynn: initial support exceeds capacity." );
  }
  return design;
}

}  // namespace fw

// --- Rcpp entry: unpack protocols, run solver, pack optimal design back to R ---
// [[Rcpp::export]]
Rcpp::List FedorovWynnAlgorithm_Rcpp(
    Rcpp::List    protocols,
    IntegerVector ndimen,
    IntegerVector nbprot,
    IntegerVector numprot,
    NumericVector freq,
    IntegerVector nbdata,
    NumericVector vectps,
    NumericVector fisher,
    IntegerVector error,
    IntegerVector protdep,
    NumericVector freqdep,
    bool          show_process = false,
    double        delta = 1e-4 )
{
  (void) nbprot;
  (void) error;

  if ( protocols.size() < 6 )
    stop( "Fedorov-Wynn: protocols list is incomplete." );
  if ( !R_finite( delta ) || delta < 0.0 )
    stop( "Fedorov-Wynn: delta must be a finite non-negative number." );

  // --- Unpack problem dimensions and candidate matrices from R list ---
  const IntegerVector n_candidates      = as<IntegerVector>( protocols[0] );
  const IntegerVector n_times_vec       = as<IntegerVector>( protocols[1] );
  const IntegerVector param_dimensions  = as<IntegerVector>( protocols[2] );
  const IntegerVector total_cost        = as<IntegerVector>( ndimen );
  const NumericMatrix sample_times      = as<NumericMatrix>( protocols[4] );
  const NumericMatrix information_rows  = as<NumericMatrix>( protocols[5] );

  const int candidate_count = n_candidates[0];
  const int param_dim       = param_dimensions[0];
  const int n_times         = n_times_vec[0];
  const int total_subjects  = total_cost[2];

  if ( candidate_count <= 0 || param_dim <= 0 || n_times <= 0 || total_subjects <= 0 )
    stop( "Fedorov-Wynn: non-positive problem dimensions." );

  // Uniform positive scale on every protocol FIM. Candidate rows already include
  // Evaluation arm-size scaling; this factor does not change the optimal simplex
  // (directional gains are invariant). Historical kernel used N / n_times.
  const double information_scale =
    static_cast<double>( total_subjects ) / static_cast<double>( n_times );

  // --- Build grid + seed design, then run outer Fedorov–Wynn loop ---
  // Pack lower-triangle rows into symmetric FIMs and apply information_scale.
  fw::ProtocolGrid grid = fw::ProtocolGrid::from_matrices(
    sample_times, information_rows, param_dim, n_times, information_scale
  );
  fw::WeightedDesign design = fw::initial_design( grid, protdep, freqdep, delta );

  fw::FedorovWynnSolver solver(
    std::move( grid ), std::move( design ), show_process
  );
  // solve(): cycles of Fedorov add, Wynn, prune, then optimality check.
  const fw::SolveStatus status = solver.solve();

  // --- Map SolveStatus to R-facing label / convergence flags ---
  const char* status_label = "success";
  switch ( status ) {
    case fw::SolveStatus::success:
      status_label = "success";
      break;
    case fw::SolveStatus::incomplete:
      status_label = "incomplete";
      break;
    case fw::SolveStatus::support_overflow:
      status_label = "support_overflow";
      break;
    case fw::SolveStatus::augmentation_failed:
      status_label = "augmentation_failed";
      break;
    case fw::SolveStatus::weight_optimization_failed:
      status_label = "weight_optimization_failed";
      break;
    case fw::SolveStatus::singular_fim:
      status_label = "singular_fim";
      break;
  }

  const bool converged = ( status == fw::SolveStatus::success );
  const bool usable =
    status == fw::SolveStatus::success || status == fw::SolveStatus::incomplete;

  // Row count matches active support (unique candidates); sized after solve.
  NumericMatrix optimal_sampling_times( 0, n_times );

  // --- Write optimal support into R inout vectors (1-based protocol indices) ---
  if ( usable ) {
    const fw::WeightedDesign& result = solver.design();
    const int active_size = result.active_size();
    const int packed_size = fw::PackedFim::packed_size( param_dim );
    if ( active_size > candidate_count )
      stop( "Fedorov-Wynn: active support exceeds candidate grid size." );
    if ( numprot.size() < active_size || freq.size() < active_size ||
         nbdata.size() < active_size || freqdep.size() < active_size )
      stop( "Fedorov-Wynn: output buffers shorter than active support." );
    if ( fisher.size() < packed_size )
      stop( "Fedorov-Wynn: fisher buffer shorter than packed FIM." );
    if ( vectps.size() < n_times )
      stop( "Fedorov-Wynn: vectps buffer shorter than n_times." );

    optimal_sampling_times = NumericMatrix( active_size, n_times );

    double weight_sum = 0.0;

    for ( int i = 0; i < active_size; ++i ) {
      const fw::ActiveProtocol& item = result.active_protocol( i );
      numprot[ i ] = item.candidate_index + 1;
      nbdata[ i ]  = static_cast<int>( item.sample_times.size() );
      // freq temporarily holds unnormalised mass (weight × information_scale).
      freq[ i ]    = item.weight * information_scale;
      weight_sum  += freq[ i ];
      freqdep[ i ] = item.weight;
    }

    if ( !R_finite( weight_sum ) || weight_sum <= 0.0 )
      stop( "Fedorov-Wynn: non-positive total weight; cannot renormalise frequencies." );

    // Renormalise frequencies to Σ freq = 1; copy sampling times row-wise.
    for ( int i = 0; i < active_size; ++i ) {
      freq[ i ] /= weight_sum;
      const fw::ActiveProtocol& item = result.active_protocol( i );
      for ( int j = 0; j < static_cast<int>( item.sample_times.size() ); ++j ) {
        vectps[ j ] = item.sample_times[ static_cast<size_t>( j ) ];
        optimal_sampling_times( i, j ) = item.sample_times[ static_cast<size_t>( j ) ];
      }
    }

    for ( int k = 0; k < packed_size; ++k )
      fisher[ k ] = result.mixture_fim()( k );
  }

  return List::create(
    Named( "freq" )                   = freq,
    Named( "optimal_sampling_times" ) = optimal_sampling_times,
    Named( "fisher" )                 = fisher,
    Named( "numprot" )                = numprot,
    Named( "status" )                 = status_label,
    Named( "converged" )              = converged
  );
}
