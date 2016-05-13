void fail(std::string const& why) {
  std::cerr << why;
  int* p = 0x0;
  std::cerr << *p;
  abort();
}

void fail_if(bool cond, std::string const& why) {
  if (cond)
    fail(why);
}
