#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <stdint.h>
#include <cstdint>
#include <stdio.h>
#include <float.h>
#include <chrono>
#include <limits>
#include <cfloat>
#include <algorithm>    




using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]


double my_zeroin(double, double, std::function <double (double)>, double, long);

struct internal_state 
  {
  uint64_t state[4];
  };

class xoshiro256_plusplus {
 
  public:
  using result_type	= std::uint64_t;
  
  xoshiro256_plusplus () 
  {
   
  }
  void init_by_splitmix64()
  {
    s[0] = splitmix64();
    s[1] = splitmix64();
    s[2] = splitmix64();
    s[3] = splitmix64();
  }
  
  void set_state(internal_state u) 
  {
    s[0] = u.state[0]; s[1] = u.state[1]; s[2] = u.state[2]; s[3] = u.state[3];
  }
  
  internal_state view() 
  {
    internal_state u;
    u.state[0] =  s[0];  u.state[1] =  s[1];  u.state[2] =  s[2];  u.state[3] =  s[3];
    return(u);
  }
  std::uint64_t operator () ()
  {
    std::uint64_t u = next();
    return(u);
  }
  inline constexpr std::uint64_t min() 
  {
  return std::numeric_limits<std::uint64_t>::lowest ();
  }
  inline constexpr std::uint64_t max() 
  {
    return  std::numeric_limits<std::uint64_t>::max ();

  }
  std::uint64_t next() 
  {
    const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
    
    const std::uint64_t t = s[1] << 17;
    
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    
    s[2] ^= t;
    
    s[3] = rotl(s[3], 45);
    
    return result;
  }
  void jump(void) // change the state by a jump
  {
    static const std::uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
    
    std::uint64_t s0 = 0;
    std::uint64_t s1 = 0;
    std::uint64_t s2 = 0;
    std::uint64_t s3 = 0;
    for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
      for(int b = 0; b < 64; b++) {
        if (JUMP[i] & UINT64_C(1) << b) {
          s0 ^= s[0];
          s1 ^= s[1];
          s2 ^= s[2];
          s3 ^= s[3];
        }
        next();	
      }
      
    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
  }
  
private:
  std::uint64_t s[4]; // state vector of the algorithm
  static inline std::uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }
  std::uint64_t splitmix64()
  {
    std::random_device rd;
    std::uint64_t x = (std::uint64_t) rd();
    std::uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
  }
};