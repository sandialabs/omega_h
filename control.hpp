void init(int& argc, char**& argv);

void fini();

void fail(std::string const& why) __attribute__((noreturn));

void fail_if(bool cond, std::string const& why);

#ifdef __CUDA_ARCH__
#define CHECK(cond) assert(cond)
#else
#define CHECK(cond) fail_if(!(cond), \
    std::string("assertion ") + #cond + " failed at " \
    __FILE__ " +" + std::to_string(__LINE__) + "\n")
#endif
