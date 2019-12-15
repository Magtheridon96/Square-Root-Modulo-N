#ifndef _SQRT_MOD_N_NUMBER_THEORETIC_H_
#define _SQRT_MOD_N_NUMBER_THEORETIC_H_

#include <algorithm>
#include <vector>
#include <tuple>
#include "utils.h"

bool is_power_of_two(int64_t n);
int64_t gcd(int64_t a, int64_t b);

// Computes a^b. Assumes b is non-negative.
int64_t ipow(int64_t a, int64_t b);

// Computes a*b (mod m).
int64_t mod_mul(int64_t a, int64_t b, int64_t m);
// Computes a^b (mod m)
int64_t mod_pow(int64_t a, int64_t b, int64_t m);

// Outputs a tuple (g, s, t) where g is the GCD of a and b,
// and s, t are such that g = a*s + b*t.
std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b);

// Computes a^(-1) (mod n)
int64_t mod_inv(int64_t a, int64_t n);

// Computes the (a p) Legendre Symbol.
int64_t legendre_symbol(int64_t a, int64_t p);
// Determines if a is a quadratic residue modulo p
int64_t is_quadratic_residue(int64_t a, int64_t p);

// Computes the square root of n modulo p.
// You can read about Tonelli-Shanks on Wikipedia.
std::vector<int64_t> tonelli_shanks(int64_t n, int64_t p);

// Chinese Remainder Theorem on a sequence (a_i, p_i).
// This will compute x such that x = a_i (mod p_i) for all i.
int64_t chinese_remainder_theorem(const std::vector<std::pair<int64_t, int64_t>>& pr);

// Determines if a 64-bit integer n is prime, fast.
// This can use Miller-Rabin for example (the implementation probably will).
bool is_prime(int64_t n);

// Finds the smallest prime factor of n, reasonably fast.
int64_t smallest_prime_factor(int64_t n);

// Generates the prime factorization of n reasonably fast.
// The output is a sequence of pairs (p_i, e_i) where p_i is a prime,
// and e_i is an exponent.
std::vector<std::pair<int64_t, int64_t>> factorize(int64_t n);

#endif
