#ifndef _SQRT_MOD_N_SQRT_MOD_N_H_
#define _SQRT_MOD_N_SQRT_MOD_N_H_

#include <set>
#include <cstdint>

// Computes the square roots of a (mod n).
// This internally is expected to use some kind of memoization to speed things up.
std::set<int64_t> sqrt_mod_n(int64_t a, int64_t n);

#endif
