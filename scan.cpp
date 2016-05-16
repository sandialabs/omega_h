template <typename TO, typename TI>
struct ExclScan : public SumFunctor<TO> {
  typedef TO value_type;
  Read<TI> in_;
  Write<TO> out_;
  ExclScan(Read<TI> in, Write<TO> out):in_(in),out_(out) {}
  INLINE void operator()(Int i, value_type& update, bool final_pass) const {
    update += in_[i];
    if (final_pass)
      out_[i + 1] = update;
  }
};


template <typename TO, typename TI>
Read<TO> offset_scan(Read<TI> a) {
  Write<TO> out(a.size() + 1);
  parallel_scan(a.size(), ExclScan<TO,TI>(a, out));
  return out;
}

template Read<I32> offset_scan(Read<I8> a);
template Read<I32> offset_scan(Read<I32> a);
