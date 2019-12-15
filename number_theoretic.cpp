#include "number_theoretic.h"

// Determines if n is a power of 2
bool is_power_of_two(int64_t n) {
  return (n & (n - 1)) == 0;
}

// Computes a^b in log(b) time.
// Assumes b is non-negative.
int64_t ipow(int64_t a, int64_t b) {
  if (b == 0) return 1;
  else if (b & 1) return a*ipow(a, b - 1);
  else return ipow(a*a, b >> 1);
}

// Computes GCD(a, b)
int64_t gcd(int64_t a, int64_t b) {
  while (b) {
    auto t = b;
    b = a % b;
    a = t;
  }
  return a;
}

// Computes a*b (mod m).
int64_t mod_mul(int64_t a, int64_t b, int64_t m) {
#ifdef SQRT_MOD_N_SAFE_AND_SLOW_MUL
  // This will be needed if you want to do
  // computations with moduli greater than 2^30 - 1.
  // This is safer, and will depend on you defining
  // the appropriate flag.
  int64_t result = 0;
  a %= m;
  b %= m;
  while (b) {
    if (b & 1) {
      result = (result + a) % m;
    }
    a = (a + a) % m;
    b >>= 1;
  }
  return result;
#else
  // A lot faster, good when m is < roughly 2*30, and assuming
  // that the inputs a, b are <= m.
  return a * b % m;
#endif
}

// Computes a^b (mod m)
int64_t mod_pow(int64_t a, int64_t b, int64_t m) {
  // It should be noted that the multiplications that happen inside
  // the loop are unsafe when m is larger than roughly 2^30, and in
  // these cases, it is recommended to either use larger integer types
  // or to use mod_mul instead.
  int64_t result = 1;
  a %= m;
  while (b) {
    if (b & 1) {
      result = mod_mul(result, a, m);
    }
    a = mod_mul(a, a, m);
    b >>= 1;
  }
  return result;
}

// Extended GCD Algorithm
std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b) {
  if (a == 0) return std::make_tuple(b, 0ll, 1ll);
  auto t = extended_gcd(b % a, a);
  auto g = std::get<0>(t);
  auto y = std::get<1>(t);
  auto x = std::get<2>(t);
  return std::make_tuple(g, x - (b/a)*y, y);
}

// Computes a^(-1) (mod n)
int64_t mod_inv(int64_t a, int64_t n) {
  auto result = extended_gcd(a, n);
  if (std::get<0>(result) == 1)
    return (n + std::get<1>(result)) % n;
  throw std::runtime_error("mod_inv args not coprime.");
}

// Computes the Legendre Symbol
int64_t legendre_symbol(int64_t a, int64_t p) {
  return mod_pow(a, (p - 1)/2, p);
}

// Determines if a is a quadratic residue modulo p
int64_t is_quadratic_residue(int64_t a, int64_t p) {
  return legendre_symbol(a, p) == 1;
}

// Computes the square root of n modulo p
// using the Tonelli-Shanks algorithm.
std::vector<int64_t> tonelli_shanks(int64_t n, int64_t p) {
  if (!is_quadratic_residue(n, p)) return {};

  int64_t q{p - 1};
  int64_t s{0};
  while (~q & 1) {
  	q >>= 1;
  	s += 1;
  }

  // p = 3 (mod 4)
  // Hence, the solutions are trivial.
  if (s == 1) {
  	const auto x{mod_pow(n, (p + 1)/4, p)};
  	return {x, p - x};
  }

  // Select a quadratic non-residue (mod p)
  // This runs in expected logarithmic time
  // given Lagrange's theorem on the number of
  // quadratic residues modulo p.
  int64_t z{0};
  for (int64_t k{1}; k < p; ++k) {
    if (!is_quadratic_residue(k, p)) {
      z = k;
      break;
    }
  }

  auto c{mod_pow(z, q, p)};
  auto r{mod_pow(n, (q + 1)/2, p)};
  auto t{mod_pow(n, q, p)};
  auto m{s};
  while (t != 1) {
    auto i{1};
    auto x{mod_mul(t, t, p)};
    while (x != 1) {
      x = mod_mul(x, x, p);
      i += 1;
    }
    const auto b{mod_pow(c, (1ll << (m - i - 1)), p)};
    r = mod_mul(r, b, p);
    c = mod_mul(b, b, p);
    t = mod_mul(t, c, p);
    m = i;
  }
  return {r, p - r};
}

// Chinese Remainder Theorem on a sequence (ai, pi).
// This will compute x such that x = ai (mod pi) for all i
int64_t chinese_remainder_theorem(const std::vector<std::pair<int64_t, int64_t>>& pr) {
  int64_t m{1};
  for (int i{0}; i < pr.size(); ++i) {
    m *= pr[i].second;
  }

  int64_t x{0};
  for (int i{0}; i < pr.size(); ++i) {
    const auto a{pr[i].first};
    const auto pk{pr[i].second};
    const auto y0{m / pk};
    const auto y1{mod_inv(y0, pk)};
    x += mod_mul(mod_mul(a, y0, m), y1, m);
    if (x >= m) x -= m;
  }
  return x;
}

// Determines if an integer is prime using
// the Miller-Rabin primality test. This is
// the deterministic variant constructed
// for integers that can fit in 64 bits.
bool is_prime(int64_t n) {
  // This set of witnesses as it is now, is totally prime.
  // There is a better and smaller set of witnesses, and I believe it
  // has either 7 or 8 elements, which I do not believe are all prime,
  // though I'm not certain about this.
  static std::vector<std::pair<int64_t, bool>> witnesses = {
    { 2, true}, { 3, true}, { 5, true},
    { 7, true}, {11, true}, {13, true},
    {17, true}, {19, true}, {23, true},
    {29, true}, {31, true}, {37, true}
  };

  if (n <= 3) return n >= 2;
  if (~n & 1) return false;

  for (auto const& a_pair : witnesses) {
    if (n == a_pair.first) return a_pair.second;
  }

  int64_t s{0};
  int64_t d{n - 1};
  while (~d & 1) {
    d >>= 1;
    ++s;
  }

  for (auto const& a_pair : witnesses) {
    const auto a = a_pair.first;
    if (mod_pow(a, d, n) != 1) {
      bool ok{true};
      auto x = mod_pow(a, d, n);
      for (int64_t r{0}; r < s; ++r) {
        if (x == n - 1) {
          ok = false;
          break;
        }
        x = mod_mul(x, x, n);
      }
      if (ok) return false;
    }
  }
  return true;
}

// Finds the smallest prime factor of n
int64_t smallest_prime_factor(int64_t n) {
  static const int64_t sieve_limit = 20000000;
  static std::vector<int64_t> smf;
  static std::vector<int64_t> primes;

  // Perform the sieve on the first
  // call to this function.
  if (smf.size() == 0) {
    smf.resize(sieve_limit + 1);
    std::fill(smf.begin(), smf.end(), sieve_limit);
    std::vector<bool> composite(sieve_limit + 1, false);
    composite[0] = composite[1] = true;
    smf[0] = smf[1] = 1;
    for (int64_t i{2}; i <= sieve_limit; ++i) {
      if (!composite[i]) {
        primes.push_back(i);
        smf[i] = i;
        for (int64_t j{2*i}; j <= sieve_limit; j += i) {
          composite[j] = true;
          smf[j] = std::min(smf[j], i);
        }
      }
    }
  }

  if (n <= sieve_limit) return smf[n];

  // Try to quickly terminate
  // in case n happens to be prime.
  if (is_prime(n)) return n;

  // Try small primes.
  for (auto p : primes) {
    if (n % p == 0) return p;
  }

  // In this case, n is a composite number divisible by some prime greater
  // than the sieve limit. It will have one prime factor less than or equal
  // to its square root, and we'll find it by slow trial division. This can
  // be replaced by a careful implementation of Pollard's Rho algorithm if
  // needed.
  int64_t p = sieve_limit + ((sieve_limit & 1)? 0 : 1);
  while (n % p) p += 2;
  return p;
}

// Generates the prime factorization of n
std::vector<std::pair<int64_t, int64_t>> factorize(int64_t n) {
  std::vector<std::pair<int64_t, int64_t>> result;
  result.reserve(9);
  while (n > 1) {
    int64_t p{smallest_prime_factor(n)};
    int64_t e{0};
    while (n % p == 0) {
      n /= p;
      ++e;
    }
    result.push_back(std::make_pair(p, e));
  }
  return result;
}
