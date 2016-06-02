template <typename T>
struct ExclScan : public SumFunctor<LO> {
  typedef LO value_type;
  Read<T> in_;
  Write<LO> out_;
  ExclScan(Read<T> in, Write<LO> out):in_(in),out_(out) {}
  OSH_INLINE void operator()(Int i, value_type& update, bool final_pass) const {
    update += in_[i];
    if (final_pass)
      out_[i + 1] = update;
  }
};

template <typename T>
LOs offset_scan(Read<T> a) {
  Write<LO> out(a.size() + 1);
  out.set(0, 0);
  parallel_scan(a.size(), ExclScan<T>(a, out));
  return out;
}

template LOs offset_scan(Read<I8> a);
template LOs offset_scan(Read<I32> a);

struct FillRight : public MaxFunctor<LO> {
  typedef LO value_type;
  Write<LO> a_;
  FillRight(Write<LO> a):a_(a) {}
  OSH_INLINE void operator()(LO i, value_type& update, bool final_pass) const {
    if (a_[i] > update)
      update = a_[i];
    if (final_pass && (a_[i] == -1))
      a_[i] = update;
  }
};

void fill_right(Write<LO> a)
{
  parallel_scan(a.size(), FillRight(a));
}
