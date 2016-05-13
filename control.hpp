void fail(std::string const& why) __attribute__((noreturn));
void fail(std::string const& why) {
  std::cerr << why;
//int* p = 0x0;
//std::cerr << *p;
  abort();
}

void fail_if(bool cond, std::string const& why);
void fail_if(bool cond, std::string const& why) {
  if (cond)
    fail(why);
}

#define CHECK(cond) fail_if(!(cond), \
    std::string("assertion ") + #cond + " failed at " \
    __FILE__ " +" + std::to_string(__LINE__) + "\n")
