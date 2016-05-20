template <typename T, Int width = 1>
Read<T> unmap(LOs a2b, Read<T> b_data);

LOs compound_maps(LOs a2b, LOs b2c);

LOs invert_permutation(LOs a2b);

LOs collect_marked(Read<I8> marks);

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

LOs get_degrees(LOs offsets);
