struct Adj {
  Adj();
  Adj(LOs ab2b_);
  Adj(LOs ab2b_, Read<I8> codes_);
  Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_);
  LOs a2ab;
  LOs ab2b;
  Read<I8> codes;
};

INLINE I8 make_code(bool is_flipped, I8 rotation, I8 which_down) {
  return static_cast<I8>((which_down << 3) | (rotation << 1) | is_flipped);
}

INLINE bool code_is_flipped(I8 code) {
  return code & 1;
}

INLINE I8 code_rotation(I8 code) {
  return (code >> 1) & 3;
}

INLINE I8 code_which_down(I8 code) {
  return (code >> 3);
}

Adj invert(Adj down, I8 nlows_per_high, LO nlows,
    Read<GO> high_globals, map::InvertMethod method);
