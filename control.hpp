void init(int& argc, char**& argv);

void fini();

void fail(char const* format, ...) __attribute__((noreturn));

#ifdef __CUDA_ARCH__
#define CHECK(cond) assert(cond)
#else
#define CHECK(cond) ((cond) ? ((void)0) : \
  fail("assertion %s failed at %s +%d\n",#cond,__FILE__,__LINE__))
#endif
