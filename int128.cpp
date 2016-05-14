void Int128::print(std::ostream& o) const {
  std::ios::fmtflags f(o.flags());
  o << std::hex << "int128(" << high << ',' << low << ')';
  o.flags(f);
}

std::ostream& operator<<(std::ostream& o, Int128 const& x) {
  x.print(o);
  return o;
}
