struct Adj {
  Adj();
  Adj(LOs ab2b_);
  Adj(LOs ab2b_, Read<I8> codes_);
  Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_);
  LOs a2ab;
  LOs ab2b;
  Read<I8> codes;
};

Adj invert(Adj down, I8 nlows_per_high, LO nlows,
    Read<GO> high_globals, map::InvertMethod method = map::BY_ATOMICS);
