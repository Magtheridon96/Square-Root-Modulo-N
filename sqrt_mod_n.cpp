#include "sqrt_mod_n.h"

#include <unordered_map>
#include <algorithm>
#include <vector>
#include "utils.h"
#include "number_theoretic.h"

// Computes the square roots of a modulo n.
// This will also memoize results.
std::set<int64_t> sqrt_mod_n(int64_t a, int64_t n) {
  static std::unordered_map<std::pair<int64_t, int64_t>, std::set<int64_t>> memo;
  if (n < 1) return {};
  if (n == 1) return {0ll};

  // Reduce a mod n first to reduce memory consumption in the dictionary,
  // and potentially allow for more usage of previously computed results.
  a %= n;
  const auto pr = std::make_pair(a, n);
  const auto it = memo.find(pr);
  if (it != memo.end()) return it->second;

  if (gcd(a, n) == 1) {
    if (is_power_of_two(n)) {
      // Powers of two are a special case.
      // For these, there are only solutions when a = 1 (mod 8).
      if (a % std::min(n, int64_t{8}) == 1) {
        int64_t k{0};
        int64_t t{n};
        while (t > 1) {
          t >>= 1;
          k += 1;
        }
        if (k == 1) return memo[pr] = {int64_t{1}};
        if (k == 2) return memo[pr] = {int64_t{1}, int64_t{3}};
        if (k == 3) return memo[pr] = {int64_t{1}, int64_t{3}, int64_t{5}, int64_t{7}};
        // Small optimization for the case of a == 1.
        if (a == 1)
          return memo[pr] = {int64_t{1}, (n >> 1) - 1ll, (n >> 1) + 1ll, n - 1ll};

        std::set<int64_t> roots;
        for (const auto x : sqrt_mod_n(a, n >> 1)) {
          const int64_t i = ((x*x - a) >> (k - 1)) & 1;
          const int64_t r = x + (i << (k - 2));
          roots.insert(r);
          roots.insert(n - r);
        }
        return memo[pr] = roots;
      }
    } else if (is_prime(n)) {
      // Primes are handled by the standard Tonelli-Shanks algorithm.
      std::set<int64_t> roots;
      for (const auto& r : tonelli_shanks(a, n)) {
        roots.insert(r);
      }
      return memo[pr] = roots;
    } else {
      // Composite numbers are not too hard. We decompose the number into
      // prime powers, solve the problem for the primes, lift the results
      // for each prime to results of prime powers using Hensel's Lifting
      // Lemma, and then we combine the results using C.R.T.
      const auto pe = factorize(n);

      // In the case of n being just an odd prime power.
      if (pe.size() == 1) {
        int64_t p{pe[0].first};
        int64_t k{pe[0].second};

        // Since n = p^k, we have to solve the equation x^2 = a (mod p), then
        // use Hensel's Lifting Lemma.
        auto roots = tonelli_shanks(a, p);
        int64_t pk{p};
        int64_t pi{p*p};
        for (int i{2}; i <= k; ++i) {
          const int64_t x{roots[0]};
          const int64_t y{mod_mul(mod_inv(2, pk), mod_inv(x, pk), pk)};
          roots[0] = (pi + x - mod_mul(mod_mul(x, x, pi) - a + pi, y, pi)) % pi;
          roots[1] = pi - roots[0];
          pk *= p;
          pi *= p;
        }
        return memo[pr] = {roots[0], roots[1]};
      } else {
        // Construct solutions for prime powers.
        std::vector<std::vector<std::pair<int64_t, int64_t>>> solutions(pe.size());
        for (int i{0}; i < pe.size(); ++i) {
          const auto m = ipow(pe[i].first, pe[i].second);
          const auto r = sqrt_mod_n(a, m);
          solutions[i].reserve(r.size());
          for (auto&& r0 : r) {
            solutions[i].push_back(std::make_pair(r0, m));
          }
        }
        // Construct all the possible square roots using
        // the Chinese Remainder Theorem.
        const auto cp = cartesian_product(solutions);
        std::set<int64_t> roots;
        for (auto&& p : cp) {
          roots.insert(chinese_remainder_theorem(p));
        }
        return memo[pr] = roots;
      }
    }
  }
  // No solutions.
  return {};
}
