#include "sqrt_mod_n.h"
#include "utils.h"
#include <iostream>

int main() {
  // Example usage.
  std::cout << sqrt_mod_n(1, 5777) << '\n';
  std::cout << sqrt_mod_n(1, 19937) << '\n';
  std::cout << sqrt_mod_n(1, 9001) << '\n';

  // Expected to be empty.
  std::cout << sqrt_mod_n(2, 3) << '\n';
  return 0;
}
