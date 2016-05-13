template <typename T>
void parallel_for(UInt n, T const& f) {
  for (UInt i = 0; i < n; ++i)
    f(i);
}
