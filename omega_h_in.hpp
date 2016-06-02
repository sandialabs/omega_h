#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#cmakedefine OSH_USE_MPI
#cmakedefine OSH_USE_KOKKOS
#cmakedefine OSH_USE_OPENMP
#cmakedefine OSH_USE_CUDA
#cmakedefine OSH_USE_ZLIB
#cmakedefine OSH_CHECK_BOUNDS

#include <cassert>
#include <memory>
#include <string>

#ifdef OSH_USE_MPI

/* on BlueGene/Q the default install
 * defines MPICH2_CONST to empty if we
 * don't define it first, causing tons
 * of compile errors.
 *
 * in addition, the mpicxx.h header is full
 * of "const MPICH2_CONST", probably as a workaround
 * for the first mistake above.
 * as a result, properly defining MPICH2_CONST to const
 * also causes compile errors.
 * luckily, we can avoid including mpicxx.h with
 * MPICH_SKIP_MPICXX.
 */
#ifdef __bgq__
#define MPICH2_CONST const
#define MPICH_SKIP_MPICXX
#endif

/* have Clang diagnostics ignore everything
 * inside mpi.h
 */
#ifdef __clang__
#pragma clang system_header
#endif
#include <mpi.h>

#endif //OSH_USE_MPI

#ifdef OSH_USE_KOKKOS
#ifdef __clang__
#pragma clang system_header
#endif
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <Kokkos_Core.hpp>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#endif //OSH_USE_KOKKOS

#ifdef OSH_USE_KOKKOS
#define OSH_INLINE KOKKOS_INLINE_FUNCTION
#define OSH_LAMBDA KOKKOS_LAMBDA
#else
#define OSH_INLINE inline
#define OSH_LAMBDA [=]
#endif //OSH_USE_KOKKOS

#ifdef OSH_USE_CUDA
#define OSH_DEVICE __device__ inline
#else
#define OSH_DEVICE inline
#endif //OSH_USE_CUDA

//namespace osh will start here

void init(int& argc, char**& argv);

void fini();

void fail(char const* format, ...) __attribute__((noreturn));

#ifdef __CUDA_ARCH__
#define OSH_CHECK(cond) assert(cond)
#else
#define OSH_CHECK(cond) ((cond) ? ((void)0) : \
  fail("assertion %s failed at %s +%d\n",#cond,__FILE__,__LINE__))
#endif

typedef std::int8_t  I8;
typedef std::int16_t I16;
typedef std::int32_t I32;
typedef std::int64_t I64;
typedef I32          Int;
typedef I32          LO;
typedef I64          GO;
typedef double       Real;

template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef OSH_USE_KOKKOS
  Kokkos::View<T*> view_;
#else
  std::shared_ptr<T> ptr_;
  LO size_;
#endif
public:
  Write();
  Write(LO size);
  Write(LO size, T value);
  Write(LO size, T offset, T stride);
  Write(HostWrite<T> host_write);
  LO size() const;
  OSH_INLINE T& operator[](LO i) const {
#ifdef OSH_CHECK_BOUNDS
    OSH_CHECK(0 <= i);
    OSH_CHECK(i < size());
#endif
#ifdef OSH_USE_KOKKOS
    return view_(i);
#else
    return ptr_.get()[i];
#endif
  }
  T* data() const;
#ifdef OSH_USE_KOKKOS
  Kokkos::View<T*> view() const;
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
  bool exists() const;
};

template <typename T>
class Read {
  Write<T> write_;
public:
  Read();
  Read(Write<T> write);
  Read(LO size, T value);
  Read(LO size, T offset, T stride);
  Read(std::initializer_list<T> l);
  LO size() const;
  OSH_INLINE T const& operator[](LO i) const {
    return write_[i];
  }
  T const* data() const;
#ifdef OSH_USE_KOKKOS
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T last() const;
  bool exists() const;
};

class LOs : public Read<LO> {
public:
  LOs();
  LOs(Read<LO> base);
  LOs(Write<LO> write);
  LOs(LO size, LO value);
  LOs(LO size, LO offset, LO stride);
  LOs(std::initializer_list<LO> l);
};

class Reals : public Read<Real> {
public:
  Reals();
  Reals(Read<Real> base);
  Reals(Write<Real> write);
  Reals(LO size, Real value);
  Reals(std::initializer_list<Real> l);
};

template <typename T>
class HostRead {
  Read<T> read_;
#ifdef OSH_USE_KOKKOS
  typename Kokkos::View<const T*>::HostMirror mirror_;
#endif
public:
  HostRead();
  HostRead(Read<T> read);
  LO size() const;
  inline T const& operator[](LO i) const {
#ifdef OSH_USE_KOKKOS
#ifdef OSH_CHECK_BOUNDS
    OSH_CHECK(0 <= i);
    OSH_CHECK(i < size());
#endif
    return mirror_(i);
#else
    return read_[i];
#endif
  }
  T const* data() const;
  T last() const;
};

enum Xfer {
  OSH_DONT_TRANSFER,
  OSH_INHERIT,
  OSH_LINEAR_INTERP,
  OSH_POINTWISE,
  OSH_CONSERVE,
  OSH_GLOBAL,
  OSH_LENGTH,
  OSH_QUALITY,
  OSH_METRIC
};

enum TagType {
  OSH_I8  = 0,
  OSH_I32 = 2,
  OSH_I64 = 3,
  OSH_F64 = 5,
};

class TagBase {
  public:
    TagBase(std::string const& name, Int ncomps, Xfer xfer);
    virtual ~TagBase();
    std::string const& name() const;
    Int ncomps() const;
    Xfer xfer() const;
    virtual TagType type() const = 0;
  private:
    std::string name_;
    Int ncomps_;
    Xfer xfer_;
};

template <typename T>
class Tag : public TagBase {
  public:
    Tag(std::string const& name, Int ncomps, Xfer xfer);
    Read<T> array() const;
    void set_array(Read<T> array);
    virtual TagType type() const override;
  private:
    Read<T> array_;
};

struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, LOs idxs_):
    ranks(ranks_),idxs(idxs_) {
  }
  Read<I32> ranks;
  LOs idxs;
};

struct Int128
{
  std::int64_t high;
  std::uint64_t low;
  OSH_INLINE Int128();
  OSH_INLINE Int128(std::int64_t h, std::uint64_t l);
  OSH_INLINE Int128(std::int64_t value);
  OSH_INLINE void operator=(Int128 const& rhs) volatile;
  OSH_INLINE Int128(Int128 const& rhs);
  OSH_INLINE Int128(const volatile Int128& rhs);
  double to_double(double unit) const;
  void print(std::ostream& o) const;
  static OSH_INLINE Int128 from_double(double value, double unit);
};

enum ReduceOp {
  MIN,
  MAX,
  SUM
};

class Comm;

typedef std::shared_ptr<Comm> CommPtr;

class Comm {
#ifdef OSH_USE_MPI
  MPI_Comm impl_;
#endif
  Read<I32> srcs_;
  HostRead<I32> host_srcs_;
  Read<I32> dsts_;
  HostRead<I32> host_dsts_;
public:
  Comm();
#ifdef OSH_USE_MPI
  Comm(MPI_Comm impl);
#else
  Comm(bool sends_to_self);
#endif
  ~Comm();
  static CommPtr world();
  static CommPtr self();
  I32 rank() const;
  I32 size() const;
  CommPtr dup() const;
  CommPtr split(I32 color, I32 key) const;
  CommPtr graph(Read<I32> dsts) const;
  CommPtr graph_adjacent(Read<I32> srcs, Read<I32> dsts) const;
  CommPtr graph_inverse() const;
  Read<I32> sources() const;
  Read<I32> destinations() const;
  template <typename T>
  T allreduce(T x, ReduceOp op) const;
  Int128 add_int128(Int128 x) const;
  template <typename T>
  T exscan(T x, ReduceOp op) const;
  template <typename T>
  void bcast(T& x) const;
  void bcast_string(std::string& s) const;
  template <typename T>
  Read<T> allgather(T x) const;
  template <typename T>
  Read<T> alltoall(Read<T> x) const;
  template <typename T>
  Read<T> alltoallv(Read<T> sendbuf,
      Read<LO> sendcounts, Read<LO> sdispls,
      Read<LO> recvcounts, Read<LO> rdispls) const;
  void barrier() const;
};

class Dist {
  CommPtr parent_comm_;
  LOs roots2items_[2];
  LOs items2content_[2];
  LOs msgs2content_[2];
  CommPtr comm_[2];
public:
  Dist();
  Dist(CommPtr comm, Remotes fitems2rroots, LO nrroots);
  void set_parent_comm(CommPtr parent_comm);
  void set_dest_ranks(Read<I32> items2ranks);
  void set_dest_idxs(LOs fitems2rroots, LO nrroots);
  void set_roots2items(LOs froots2fitems);
  Dist invert() const;
  template <typename T>
  Read<T> exch(Read<T> data, Int width) const;
  template <typename T>
  Read<T> exch_sum(Read<T> data, Int width) const;
  CommPtr parent_comm() const;
  CommPtr comm() const;
  LOs content2msgs() const;
  LOs items2msgs() const;
  LOs roots2items() const;
  Read<I32> msgs2ranks() const;
  Read<I32> items2ranks() const;
  LOs items2dest_idxs() const;
  Remotes items2dests() const;
  LO nitems() const;
  LO nroots() const;
  LO ndests() const;
  void change_comm(CommPtr new_comm);
  Remotes exch(Remotes data, Int width) const;
private:
  enum { F, R };
};

namespace inertia {
struct Rib;
}

enum Partition {
  ELEMENT_BASED,
  GHOSTED,
  VERTEX_BASED
};

//namespace osh will end here

#endif
