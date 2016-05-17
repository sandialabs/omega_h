template <typename T>
Read<T> permute(LOs out2in, Read<T> in);

LOs invert_funnel(LOs ab2a, LO na);

void invert_map_by_sorting(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a);
void invert_map_by_atomics(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a);
LOs order_by_globals(
    LOs a2ab, LOs ab2b, Read<GO> b_global);
