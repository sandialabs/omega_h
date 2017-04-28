#include <Omega_h_align.hpp>


namespace Omega_h {

template <Int deg, typename T>
static Read<T> align_ev2v_deg(Read<T> ev2v, Read<I8> codes) {
  OMEGA_H_CHECK(ev2v.size() == codes.size() * deg);
  auto ne = codes.size();
  Write<T> ev2v_w(ev2v.size());
  auto f = OMEGA_H_LAMBDA(LO e) {
    align_adj<deg>(codes[e], ev2v, e * deg, ev2v_w, e * deg);
  };
  parallel_for(ne, f);
  return ev2v_w;
}

template <typename T>
Read<T> align_ev2v(Int deg, Read<T> ev2v, Read<I8> codes) {
  if (deg == 3) {
    return align_ev2v_deg<3>(ev2v, codes);
  } else {
    OMEGA_H_CHECK(deg == 2);
    return align_ev2v_deg<2>(ev2v, codes);
  }
}

#define INST(T)                                                                \
  template Read<T> align_ev2v(Int deg, Read<T> ev2v, Read<I8> codes);
INST(LO)
INST(GO)
#undef INST

}  // end namespace Omega_h
