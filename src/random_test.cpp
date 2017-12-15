#include <Omega_h_random.hpp>

int main() {
  constexpr int nbuckets = 10;
  int nsamples = 10 * 1000;
  RandomKey key = {199, 808};
  RandomState state;
  for (int i = 0; i < nsamples; ++i) {
    auto val = get_next_random_real(key, i, state);
    OMEGA_H_CHECK(0.0 <= val);
    OMEGA_H_CHECK(val <= 1.0);
  }
}
