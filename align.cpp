template <Int deg, typename T>
Read<T> align_ev2v(Read<T> ev2v, Read<I8> codes) {
  CHECK(ev2v.size() == codes.size() * deg);
  auto ne = codes.size();
  Write<T> ev2v_w(ev2v.size());
  auto f = LAMBDA(LO e) {
    align_adj<deg>(codes[e], &ev2v[e * deg], &ev2v_w[e * deg]);
  };
  parallel_for(ne, f);
  return ev2v_w;
}

#define INST(deg,T) \
template Read<T> align_ev2v<deg>(Read<T> ev2v, Read<I8> codes);
INST(2,LO)
INST(3,LO)
INST(2,GO)
INST(3,GO)
#undef INST
