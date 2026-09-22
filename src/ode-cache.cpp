// ODE simulation-time grid cache and small hashing helpers.
//
// Builds the unique sorted time grid (sampling times ∪ event times) used by
// deSolve, with an optional process-wide LRU cache keyed by a string from R.
// Also exports cache_hash / env_cache_key wrappers used by R performance caches.
// [[Rcpp::depends(Rcpp)]]
#include <pfim/pfim-cache.hpp>
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// Process-wide LRU map: cache_key -> sorted unique simulation times.
struct OdeSimTimesCache {
  std::unordered_map<std::string, std::vector<double>> data_;
  std::list<std::string> lru_;

  // Move `key` to the most-recently-used end of the LRU list.
  void touch( const std::string& key ) {
    for ( auto it = lru_.begin(); it != lru_.end(); ++it ) {
      if ( *it == key ) {
        lru_.erase( it );
        break;
      }
    }
    lru_.push_back( key );
  }

  // Evict oldest entries until size <= max_entries (no-op if max_entries <= 0).
  void trim( int max_entries ) {
    if ( max_entries <= 0 )
      return;
    while ( static_cast<int>( lru_.size() ) > max_entries ) {
      const std::string drop = lru_.front();
      data_.erase( drop );
      lru_.pop_front();
    }
  }

  const std::vector<double>* get( const std::string& key ) {
    const auto it = data_.find( key );
    if ( it == data_.end() )
      return nullptr;
    touch( key );
    return &it->second;
  }

  void set( const std::string& key, std::vector<double> value, int max_entries ) {
    data_[ key ] = std::move( value );
    touch( key );
    trim( max_entries );
  }

  void clear() {
    data_.clear();
    lru_.clear();
  }

  int size() const { return static_cast<int>( data_.size() ); }
};

OdeSimTimesCache& ode_sim_times_cache() {
  static OdeSimTimesCache cache;
  return cache;
}

}  // namespace

// Hash character parts joined by a unit-separator (R: pfimCacheHash_Rcpp).
// [[Rcpp::export(name = "pfimCacheHash_Rcpp")]]
std::string pfimCacheHash_Rcpp( const Rcpp::CharacterVector& parts ) {
  return pfim::cache_hash_parts( parts );
}

// Shorten long environment binding names (Windows limit); pass-through if short.
// [[Rcpp::export(name = "pfimEnvCacheKey_Rcpp")]]
std::string pfimEnvCacheKey_Rcpp( const std::string& x ) {
  return pfim::env_cache_key( x );
}

// Return (and optionally cache) the unique sorted ODE simulation time grid.
// [[Rcpp::export(name = "pfimOdeSimTimesCached_Rcpp")]]
Rcpp::NumericVector pfimOdeSimTimesCached_Rcpp( const std::string& cache_key,
                                                const Rcpp::NumericVector& raw_samplings,
                                                const Rcpp::Nullable<Rcpp::NumericVector>& event_times,
                                                bool enabled,
                                                int max_entries = -1 ) {
  auto& cache = ode_sim_times_cache();
  if ( enabled ) {
    if ( const std::vector<double>* hit = cache.get( cache_key ) )
      return Rcpp::NumericVector( hit->begin(), hit->end() );
  }

  const std::vector<double> times = pfim::build_ode_sim_times( raw_samplings, event_times );
  if ( enabled )
    cache.set( cache_key, times, max_entries );

  return Rcpp::NumericVector( times.begin(), times.end() );
}

// [[Rcpp::export(name = "pfimOdeSimTimesCacheClear_Rcpp")]]
void pfimOdeSimTimesCacheClear_Rcpp() {
  ode_sim_times_cache().clear();
}

// [[Rcpp::export(name = "pfimOdeSimTimesCacheSize_Rcpp")]]
int pfimOdeSimTimesCacheSize_Rcpp() {
  return ode_sim_times_cache().size();
}
