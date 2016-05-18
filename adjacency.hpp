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

/* given entity uses and unique entities,
   both defined by vertex lists, match
   uses to unique entities and derive their
   respective alignment codes.

   even though this is a downward adjacency, we'll
   define the code as describing how to transform
   the boundary entity into the entity use,
   since typically data is being pulled into an element
   from its boundary

   this function requires that all entities have
   at least one corresponding use ! */
template <Int deg>
void reflect_down(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);
