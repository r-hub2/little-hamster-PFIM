// Cache hashing and ODE simulation-time grid helpers (header-only).
//
// Used by ode-cache.cpp and R performance caches: join/hash string parts,
// shorten long environment keys, and build unique sorted sim times from
// sampling times plus optional event times.
//
// Pitfall: near-duplicate times closer than 1e-8 are dropped — keep this
// tolerance aligned with any R-side grid builders that feed the same ODE path.
#pragma once
#include <Rcpp.h>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace pfim {

// Hex digest of std::hash over an arbitrary string buffer.
inline std::string cache_hash_string( const std::string& buf ) {
  const size_t h = std::hash<std::string>{}( buf );
  std::ostringstream oss;
  oss << std::hex << h;
  return oss.str();
}

// Join non-NA character parts with 0x1F, then hash (stable across call order).
inline std::string cache_hash_parts( const Rcpp::CharacterVector& parts ) {
  std::string buf;
  for ( int i = 0; i < parts.size(); ++i ) {
    if ( Rcpp::CharacterVector::is_na( parts[ i ] ) )
      continue;
    if ( !buf.empty() )
      buf += '\x1f';
    buf += Rcpp::as<std::string>( parts[ i ] );
  }
  return cache_hash_string( buf );
}

// Pass through short keys; hash long ones so Windows env binding names stay valid.
inline std::string env_cache_key( const std::string& x ) {
  if ( x.size() <= 500 )
    return x;
  return "h" + cache_hash_string( x );
}

// Unique sorted simulation times = samples ∪ events, dropping near-duplicates.
inline std::vector<double> build_ode_sim_times( const Rcpp::NumericVector& raw_samplings,
                                                const Rcpp::Nullable<Rcpp::NumericVector>& event_times ) {
  std::vector<double> times( raw_samplings.begin(), raw_samplings.end() );
  if ( event_times.isNotNull() ) {
    const Rcpp::NumericVector et = event_times.get();
    times.insert( times.end(), et.begin(), et.end() );
  }
  if ( times.empty() )
    return times;

  std::sort( times.begin(), times.end() );
  times.erase( std::unique( times.begin(), times.end() ), times.end() );

  // Drop points closer than 1e-8 to the previous kept time (solver stability).
  std::vector<double> out;
  out.reserve( times.size() );
  out.push_back( times.front() );
  for ( size_t i = 1; i < times.size(); ++i ) {
    if ( times[ i ] - out.back() > 1e-8 )
      out.push_back( times[ i ] );
  }
  return out;
}

}  // namespace pfim
