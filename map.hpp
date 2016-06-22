template <typename T>
void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width);

template <typename T>
Read<T> map_onto(Read<T> a_data, LOs a2b, LO nb, T init_val, Int width);

template <typename T>
Read<T> unmap(LOs a2b, Read<T> b_data, Int width);

template <typename T>
Read<T> expand(Read<T> a_data, LOs a2b, Int width);

LOs multiply_fans(LOs a2b, LOs a2c);

LOs compound_maps(LOs a2b, LOs b2c);

LOs invert_permutation(LOs a2b);

Read<I8> invert_marks(Read<I8> marks);

LOs collect_marked(Read<I8> marks);

Read<I8> mark_image(LOs a2b, LO nb);

LOs invert_injective_map(LOs a2b, LO nb);

LOs invert_funnel(LOs ab2a, LO na);

namespace map {

enum InvertMethod {
  BY_SORTING,
  BY_ATOMICS
};

Graph invert_by_sorting(LOs a2b, LO nb);
Graph invert_by_atomics(LOs a2b, LO nb);
Graph invert(LOs a2b, LO nb, InvertMethod method);

}

LOs get_degrees(LOs offsets);

LOs invert_fan(LOs a2b);

template <typename T>
Read<T> fan_sum(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_max(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_min(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_reduce(LOs a2b, Read<T> b_data, Int width, osh_op op);
