#include "Omega_h_future.hpp"

namespace Omega_h {

#if defined(OMEGA_H_USE_MPI) && defined(OMEGA_H_USE_CUDA) &&                   \
    !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
template <typename T>
Future<T>::Future(HostRead<T> sendbuf, HostWrite<T> recvbuf,
    const requests_type&& requests, const callback_type callback)
    : sendbuf_(sendbuf),
      recvbuf_(recvbuf),
      callback_(callback),
      requests_(requests),
      status_(Status::waiting) {}
#else

template <typename T>
Future<T>::Future(Read<T> sendbuf, Write<T> recvbuf,
    const requests_type&& requests, const callback_type callback)
    : sendbuf_(sendbuf),
      recvbuf_(recvbuf),
      callback_(callback),
      requests_(requests),
      status_(Status::waiting) {}

/// no-op: no comm, no callback
template <typename T>
Future<T>::Future(Write<T> recvbuf) : Future({} /* sendbuf */, recvbuf) {
  status_ = Status::completed;
}

#endif

/// no-op: no comm, no callback
template <typename T>
Future<T>::Future(Read<T> buf)
    : Future({} /* sendbuf */, {} /* recvbuf */, {} /* request */,
          [buf](recvbuf_type) { return buf; }) {
  status_ = Status::completed;
}

template <typename T>
bool Future<T>::completed() {
#ifdef OMEGA_H_USE_MPI
  if (status_ != Status::waiting) {
    return true;
  }
  int flag = 0;
  OMEGA_H_CHECK(
      MPI_SUCCESS == MPI_Testall(static_cast<int>(requests_.size()),
                         requests_.data(), &flag, MPI_STATUSES_IGNORE));
  if (flag == 0) {
    status_ = Status::completed;
  }
  return flag == 0;
#else   // !OMEGA_H_USE_MPI
  return true;
#endif  // OMEGA_H_USE_MPI
}

template <typename T>
Read<T> Future<T>::get() {
  if (status_ == Status::completed) {
    status_ = Status::consumed;
    return callback_(recvbuf_);
#ifdef OMEGA_H_USE_MPI
  } else if (status_ == Status::waiting) {
    OMEGA_H_CHECK(MPI_SUCCESS == MPI_Waitall(static_cast<int>(requests_.size()),
                                     requests_.data(), MPI_STATUS_IGNORE));
    status_ = Status::consumed;
    return callback_(recvbuf_);
#endif  // !OMEGA_H_USE_MPI
  } else {
    fail("Can not ask the result more than once.");
  }
}

#if defined(OMEGA_H_USE_MPI) && defined(OMEGA_H_USE_CUDA) &&                   \
    !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#define OMEGA_H_INST(T)                                                        \
  template Future<T>::Future(HostRead<T> sendbuf, HostWrite<T> recvbuf,        \
      const requests_type&&, const callback_type callback);
#else
#define OMEGA_H_INST(T)                                                        \
  template Future<T>::Future(Read<T> sendbuf, Write<T> recvbuf,                \
      const requests_type&& requests, const callback_type callback);           \
  template Future<T>::Future(Write<T> recvbuf);
#endif
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

#define OMEGA_H_INST(T)                                                        \
  template Future<T>::Future(Read<T> buf);                                     \
  template bool Future<T>::completed();                                        \
  template Read<T> Future<T>::get();
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

}  // namespace Omega_h
