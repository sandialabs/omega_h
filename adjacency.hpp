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
    Read<GO> high_globals);

/* given the vertex lists for high entities,
   create vertex lists for all uses of low
   entities by high entities */
LOs form_uses(LOs hv2v, I8 high_dim, I8 low_dim);

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
void find_matches_by_sorting(LOs euv2v, LOs ev2v,
    LOs& eu2e, Read<I8>& eu2e_codes_);

Adj reflect_down_by_sorting(LOs hv2v, LOs lv2v, I8 high_dim, I8 low_dim);

LOs find_unique_by_sorting(LOs hv2v, I8 high_dim, I8 low_dim);

/* for each entity (or entity use), sort its vertex list.
   express the sorting transformation as
   an alignment code, output those too */
template <Int deg>
void make_canonical(LOs ev2v,
    LOs& canon_, Read<I8>& codes_);
template <Int deg>
Read<I8> find_jumps(LOs canon, LOs e_sorted2e);

template <Int deg>
void find_matches_by_upward(LOs av2v, LOs bv2v, Adj v2b,
    LOs& a2b, Read<I8>& codes);

Adj reflect_down_by_upward(LOs hv2v, LOs lv2v, Adj v2l,
    I8 high_dim, I8 low_dim);

namespace adj {
enum ReflectMethod {
  BY_SORTING,
  BY_UPWARD
};
}
Adj reflect_down(LOs hv2v, LOs lv2v, LO nv,
    I8 high_dim, I8 low_dim, adj::ReflectMethod method);
