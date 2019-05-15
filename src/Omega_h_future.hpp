#ifndef OMEGA_H_Future_HPP
#define OMEGA_H_Future_HPP

#include <functional>
#include <vector>

#include <Omega_h_mpi.h>
#include "Omega_h_array.hpp"
#include "Omega_h_fail.hpp"

namespace Omega_h {

/**
 * \brief Abstraction for asynchronous communication.
 *
 * Hold both receive/send buffers and post-communication callbacks
 */
template <typename T>
class Future {
 public:
#if defined(OMEGA_H_USE_MPI) && defined(OMEGA_H_USE_CUDA) &&                   \
    !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  using sendbuf_type = HostRead<T>;
  using recvbuf_type = HostWrite<T>;
#else
  using sendbuf_type = Read<T>;
  using recvbuf_type = Write<T>;
#endif

  /// post-processing callback for asynchronous communication
  using callback_type = std::function<Read<T>(recvbuf_type)>;
  /// optional additional processing callbacks (see \a add_callback method)
  using additional_callback_type = std::function<Read<T>(Read<T>)>;

#ifdef OMEGA_H_USE_MPI
  using requests_type = std::vector<MPI_Request>;
#else
  using requests_type = std::vector<int>;
#endif  // OMEGA_H_USE_MPI

#if defined(OMEGA_H_USE_MPI) && defined(OMEGA_H_USE_CUDA) &&                   \
    !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  Future(HostRead<T> sendbuf, HostWrite<T> recvbuf,
      const requests_type&& requests, const callback_type callback);
#else
  Future(Read<T> sendbuf, Write<T> recvbuf, const requests_type&& requests = {},
      const callback_type callback = [](recvbuf_type buf) -> Read<T> {
        return buf;
      });

  /// no-op: no comm, no callback
  explicit Future(Write<T> recvbuf);

#endif
  /// no-op: no comm, no callback
  explicit Future(Read<T> buf);

  /// register a additional post-communication callback
  void add_callback(const additional_callback_type& callback) {
    auto cur_cb = callback_;
    callback_ = [=](recvbuf_type recvbuf) -> Read<T> {
      return callback(cur_cb(recvbuf));
    };
  }

  /// \return true if asynchronous operation completed, false otherwise
  bool completed();

  /// wait for asynchronous operation completion, and return the result.
  /// This method can only be called once.
  Read<T> get();

 private:
  /// Object state
  enum class Status {
    waiting,    /// communication in progress
    completed,  /// communication is over, waiting for user to retrieve the
                /// result
    consumed    /// result already given to the user
  };

  sendbuf_type sendbuf_;
  recvbuf_type recvbuf_;
  callback_type callback_;
  requests_type requests_;
  Status status_;
};

#if defined(OMEGA_H_USE_MPI) && defined(OMEGA_H_USE_CUDA) &&                   \
    !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Future<T>::Future(HostRead<T> sendbuf, HostWrite<T> recvbuf, \
      const requests_type&&, const callback_type callback);
#else
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Future<T>::Future(Read<T> sendbuf, Write<T> recvbuf,         \
      const requests_type&& requests, const callback_type callback);           \
  extern template Future<T>::Future(Write<T> recvbuf);
#endif
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Future<T>::Future(Read<T> buf);                              \
  extern template bool Future<T>::completed();                                 \
  extern template Read<T> Future<T>::get();

OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // namespace Omega_h

#endif  // OMEGA_H_Future_HPP
