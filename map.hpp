template <typename T, Int width = 1>
Read<T> permute(LOs out2in, Read<T> in);

LOs compound_maps(LOs a2b, LOs b2c);

LOs invert_permutation(LOs a2b);

LOs invert_funnel(LOs ab2a, LO na);

namespace map {

enum InvertMethod {
  BY_SORTING,
  BY_ATOMICS
};

void invert_by_sorting(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a);
void invert_by_atomics(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a);
void invert(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a, InvertMethod method);

}

LOs order_by_globals(
    LOs a2ab, LOs ab2b, Read<GO> b_global);
