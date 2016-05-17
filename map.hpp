template <typename T>
Read<T> permute(LOs out2in, Read<T> in);

LOs invert_funnel(LOs ab2a, LO na);

void invert_map_by_sort(LOs a2b, LO nb,
    LOs& b2ba, LOs& ba2a);
