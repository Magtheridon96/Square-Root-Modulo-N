#include "utils.h"

// Outputs a pair of whatever to a stream.
template<typename T, typename U>
std::ostream& operator<<(std::ostream& stream, const std::pair<T, U>& p) {
  std::cout << "(" << p.first << ", " << p.second << ")";
}

// Outputs a vector to a stream.
template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& container) {
  if (container.size() == 0) return stream << "{}";
  auto it = container.begin();
  stream << '{' << *it++;
  while (it != container.end()) {
    stream << ", " << *it++;
  }
  return stream << '}';
}

// Outputs a set to a stream.
template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::set<T>& container) {
  if (container.size() == 0) return stream << "{}";
  auto it = container.begin();
  stream << '{' << *it++;
  while (it != container.end()) {
    stream << ", " << *it++;
  }
  return stream << '}';
}

// Cartesian Product helper function
template<typename T>
void cartesian_product_helper(
  std::vector<std::vector<T>> const& v,
  std::vector<std::vector<T>>& result,
  std::vector<T>& path,
  int i)
{
  if (i == v.size()) {
    result.push_back(path);
  } else {
    for (int j = 0; j < v[i].size(); ++j) {
      path.push_back(v[i][j]);
      cartesian_product_helper(v, result, path, i + 1);
      path.pop_back();
    }
  }
}

// Generate the cartesian product of a bunch of vectors.
template<typename T>
std::vector<std::vector<T>> cartesian_product(std::vector<std::vector<T>> const& v) {
  std::vector<std::vector<T>> result;
  std::vector<T> path;
  result.reserve(50);
  cartesian_product_helper(v, result, path, 0);
  return result;
}

namespace std {
  template<>
  class hash<std::pair<int64_t, int64_t>> {
  public:
    std::size_t operator()(const std::pair<int64_t, int64_t>& t) const {
      static auto h = std::hash<int64_t>();
      return h(t.first) ^ h(t.second);
    }
  };
}
