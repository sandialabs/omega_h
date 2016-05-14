#ifdef USE_CUDA
#include <thrust/complex.h>
typedef thrust::complex<Real> Complex;
#else
#include <complex>
typedef std::complex<Real> Complex;
#endif
