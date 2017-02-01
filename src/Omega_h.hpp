#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#include <cassert>
#include <iosfwd>
#include <map>
#include <new>
#include <string>
#include <type_traits>
#include <vector>

#include <Omega_h_c.h>
#include <Omega_h_array.hpp>

namespace Omega_h {

template <typename T>
Read<T> permute(Read<T> a_data, LOs a2b, Int width);

class TagBase {
 public:
  TagBase(std::string const& name, Int ncomps, Int xfer, Int outflags);
  virtual ~TagBase();
  std::string const& name() const;
  Int ncomps() const;
  Int xfer() const;
  Int outflags() const;
  virtual Omega_h_Type type() const = 0;

 private:
  std::string name_;
  Int ncomps_;
  Int xfer_;
  Int outflags_;
};

template <typename T>
class Tag : public TagBase {
 public:
  Tag(std::string const& name, Int ncomps, Int xfer, Int outflags);
  Read<T> array() const;
  void set_array(Read<T> array);
  virtual Omega_h_Type type() const override;

 private:
  Read<T> array_;
};

struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, LOs idxs_) : ranks(ranks_), idxs(idxs_) {}
  Read<I32> ranks;
  LOs idxs;
};

struct Int128 {
  std::int64_t high;
  std::uint64_t low;
  OMEGA_H_INLINE Int128();
  OMEGA_H_INLINE Int128(std::int64_t h, std::uint64_t l);
  OMEGA_H_INLINE Int128(std::int64_t value);
  OMEGA_H_INLINE void operator=(Int128 const& rhs) volatile;
  OMEGA_H_INLINE Int128(Int128 const& rhs);
  OMEGA_H_INLINE Int128(const volatile Int128& rhs);
  double to_double(double unit) const;
  void print(std::ostream& o) const;
  static OMEGA_H_INLINE Int128 from_double(double value, double unit);
};

class Library;
class Comm;

typedef std::shared_ptr<Comm> CommPtr;

class Comm {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl_;
#endif
  Library* library_;
  Read<I32> srcs_;
  Read<I32> dsts_;
  HostRead<I32> host_srcs_;
  HostRead<I32> host_dsts_;
  LO self_src_;
  LO self_dst_;

 public:
  Comm();
#ifdef OMEGA_H_USE_MPI
  Comm(Library* library, MPI_Comm impl);
#else
  Comm(Library* library, bool is_graph, bool sends_to_self);
#endif
  ~Comm();
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
  T allreduce(T x, Omega_h_Op op) const;
  bool reduce_or(bool x) const;
  bool reduce_and(bool x) const;
  Int128 add_int128(Int128 x) const;
  template <typename T>
  T exscan(T x, Omega_h_Op op) const;
  template <typename T>
  void bcast(T& x) const;
  void bcast_string(std::string& s) const;
  template <typename T>
  Read<T> allgather(T x) const;
  template <typename T>
  Read<T> alltoall(Read<T> x) const;
  template <typename T>
  Read<T> alltoallv(Read<T> sendbuf, Read<LO> sendcounts, Read<LO> sdispls,
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
  Dist(Dist const& other);
  Dist& operator=(Dist const& other);
  Dist(CommPtr comm, Remotes fitems2rroots, LO nrroots);
  void set_parent_comm(CommPtr parent_comm);
  void set_dest_ranks(Read<I32> items2ranks);
  void set_dest_idxs(LOs fitems2rroots, LO nrroots);
  void set_roots2items(LOs froots2fitems);
  Dist invert() const;
  template <typename T>
  Read<T> exch(Read<T> data, Int width) const;
  template <typename T>
  Read<T> exch_reduce(Read<T> data, Int width, Omega_h_Op op) const;
  CommPtr parent_comm() const;
  CommPtr comm() const;
  LOs msgs2content() const;
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
  LO nsrcs() const;
  void change_comm(CommPtr new_comm);
  Remotes exch(Remotes data, Int width) const;

 private:
  void copy(Dist const& other);
  enum { F, R };
};

struct Graph {
  Graph() {}
  explicit Graph(LOs ab2b_) : ab2b(ab2b_) {}
  Graph(LOs a2ab_, LOs ab2b_) : a2ab(a2ab_), ab2b(ab2b_) {}
  LOs a2ab;
  LOs ab2b;
  LO nnodes() const;
  LO nedges() const;
};

enum { DIMS = OMEGA_H_DIMS };

enum {
  VERT = OMEGA_H_VERT,
  EDGE = OMEGA_H_EDGE,
  TRI = OMEGA_H_TRI,
  TET = OMEGA_H_TET
};

struct Adj : public Graph {
  Adj() {}
  explicit Adj(LOs ab2b_) : Graph(ab2b_) {}
  Adj(LOs ab2b_, Read<I8> codes_) : Graph(ab2b_), codes(codes_) {}
  Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_)
      : Graph(a2ab_, ab2b_), codes(codes_) {}
  Adj(LOs a2ab_, LOs ab2b_) : Graph(a2ab_, ab2b_) {}
  Adj(Graph g) : Graph(g) {}
  Read<I8> codes;
};

void find_matches(
    Int dim, LOs av2v, LOs bv2v, Adj v2b, LOs* a2b_out, Read<I8>* codes_out);

class Library {
 public:
  Library(Library const&);
  inline Library(int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi = MPI_COMM_WORLD
#endif
      ) {
    initialize(OMEGA_H_VERSION, argc, argv
#ifdef OMEGA_H_USE_MPI
        ,
        comm_mpi
#endif
        );
  }
  ~Library();
  CommPtr world();
  CommPtr self();
  void add_to_timer(std::string const& name, double nsecs);
  LO self_send_threshold() const;

 private:
  void initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi
#endif
      );
  CommPtr world_;
  CommPtr self_;
#ifdef OMEGA_H_USE_MPI
  bool we_called_mpi_init;
#endif
#ifdef OMEGA_H_USE_KOKKOS
  bool we_called_kokkos_init;
#endif
  bool should_time_;
  std::map<std::string, double> timers;
  LO self_send_threshold_;
};

namespace inertia {
struct Rib;
}

class Mesh {
 public:
  Mesh(Library* library);
  Library* library() const;
  void set_comm(CommPtr const& comm);
  void set_dim(Int dim);
  void set_verts(LO nverts);
  void set_ents(Int dim, Adj down);
  void keep_canonical_globals(bool yn);
  CommPtr comm() const;
  Omega_h_Parting parting() const;
  inline Int dim() const {
    OMEGA_H_CHECK(0 <= dim_ && dim_ <= 3);
    return dim_;
  }
  LO nents(Int dim) const;
  LO nelems() const;
  LO ntets() const;
  LO ntris() const;
  LO nedges() const;
  LO nverts() const;
  GO nglobal_ents(Int dim);
  template <typename T>
  void add_tag(
      Int dim, std::string const& name, Int ncomps, Int xfer, Int outflags);
  template <typename T>
  void add_tag(Int dim, std::string const& name, Int ncomps, Int xfer,
      Int outflags, Read<T> array, bool internal = false);
  template <typename T>
  void set_tag(
      Int dim, std::string const& name, Read<T> array, bool internal = false);
  TagBase const* get_tagbase(Int dim, std::string const& name) const;
  template <typename T>
  Tag<T> const* get_tag(Int dim, std::string const& name) const;
  template <typename T>
  Read<T> get_array(Int dim, std::string const& name) const;
  void remove_tag(Int dim, std::string const& name);
  bool has_tag(Int dim, std::string const& name) const;
  Int ntags(Int dim) const;
  TagBase const* get_tag(Int dim, Int i) const;
  bool has_ents(Int dim) const;
  bool has_adj(Int from, Int to) const;
  Adj get_adj(Int from, Int to) const;
  Adj ask_down(Int from, Int to);
  LOs ask_verts_of(Int dim);
  LOs ask_elem_verts();
  Adj ask_up(Int from, Int to);
  Graph ask_star(Int dim);
  Graph ask_dual();

 public:
  typedef std::shared_ptr<TagBase> TagPtr;
  typedef std::shared_ptr<Adj> AdjPtr;
  typedef std::shared_ptr<Dist> DistPtr;
  typedef std::shared_ptr<inertia::Rib> RibPtr;

 private:
  typedef std::vector<TagPtr> TagVector;
  typedef TagVector::iterator TagIter;
  typedef TagVector::const_iterator TagCIter;
  TagIter tag_iter(Int dim, std::string const& name);
  TagCIter tag_iter(Int dim, std::string const& name) const;
  void check_dim(Int dim) const;
  void check_dim2(Int dim) const;
  void add_adj(Int from, Int to, Adj adj);
  Adj derive_adj(Int from, Int to);
  Adj ask_adj(Int from, Int to);
  void react_to_set_tag(Int dim, std::string const& name);
  Int dim_;
  CommPtr comm_;
  Int parting_;
  Int nghost_layers_;
  LO nents_[DIMS];
  TagVector tags_[DIMS];
  AdjPtr adjs_[DIMS][DIMS];
  Remotes owners_[DIMS];
  DistPtr dists_[DIMS];
  RibPtr rib_hints_;
  bool keeps_canonical_globals_;
  Library* library_;

 public:
  void add_coords(Reals array);
  Reals coords() const;
  void set_coords(Reals const& array);
  Read<GO> ask_globals(Int dim);
  void reset_globals();
  Reals ask_lengths();
  Reals ask_qualities();
  void set_owners(Int dim, Remotes owners);
  Remotes ask_owners(Int dim);
  Read<I8> owned(Int dim);
  Dist ask_dist(Int dim);
  Int nghost_layers() const;
  void set_parting(Omega_h_Parting parting, Int nlayers, bool verbose);
  void set_parting(Omega_h_Parting parting, bool verbose = false);
  void migrate(Remotes new_elems2old_owners, bool verbose = false);
  void reorder();
  void balance(bool predictive = false);
  Graph ask_graph(Int from, Int to);
  template <typename T>
  Read<T> sync_array(Int ent_dim, Read<T> a, Int width);
  template <typename T>
  Read<T> sync_subset_array(
      Int ent_dim, Read<T> a_data, LOs a2e, T default_val, Int width);
  template <typename T>
  Read<T> reduce_array(Int ent_dim, Read<T> a, Int width, Omega_h_Op op);
  template <typename T>
  Read<T> owned_array(Int ent_dim, Read<T> a, Int width);
  void sync_tag(Int dim, std::string const& name);
  void reduce_tag(Int dim, std::string const& name, Omega_h_Op op);
  bool operator==(Mesh& other);
  Real min_quality();
  Real max_length();
  bool could_be_shared(Int ent_dim) const;
  bool owners_have_all_upward(Int ent_dim) const;
  Mesh copy_meta() const;
  bool keeps_canonical_globals() const;
  RibPtr rib_hints() const;
  void set_rib_hints(RibPtr hints);
  Real imbalance(Int ent_dim = -1) const;
};

#ifdef OMEGA_H_USE_LIBMESHB
namespace meshb {
void read(Mesh* mesh, const char* filepath);
void write(Mesh* mesh, const char* filepath, int version);
}
#endif

namespace gmsh {
void read(std::istream& stream, Mesh* mesh);
void read(std::string const& filename, Mesh* mesh);
}

namespace vtk {
void write_vtu(std::ostream& stream, Mesh* mesh, Int cell_dim);
void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim);
void write_parallel(std::string const& path, Mesh* mesh, Int cell_dim);
class Writer {
  Mesh* mesh_;
  std::string root_path_;
  Int cell_dim_;
  Int step_;
  std::streampos pvd_pos_;

 public:
  Writer();
  Writer(Writer const&);
  Writer& operator=(Writer const&);
  ~Writer();
  Writer(Mesh* mesh, std::string const& root_path, Int cell_dim);
  void write(Real time);
  void write();
};
class FullWriter {
  std::vector<Writer> writers_;

 public:
  FullWriter();
  FullWriter(FullWriter const&);
  FullWriter& operator=(FullWriter const&);
  ~FullWriter();
  FullWriter(Mesh* mesh, std::string const& root_path);
  void write(Real time);
  void write();
};
}  // end namespace vtk

enum Verbosity { SILENT, EACH_ADAPT, EACH_REBUILD, EXTRA_STATS };

#ifdef OMEGA_H_USE_EGADS
struct Egads;
#endif

struct AdaptOpts {
  AdaptOpts(Int dim);     // sets defaults
  AdaptOpts(Mesh* mesh);  // calls above
  Real min_length_desired;
  Real max_length_desired;
  Real max_length_allowed;
  Real min_quality_allowed;
  Real min_quality_desired;
  Int nsliver_layers;
  Verbosity verbosity;
  Real length_histogram_min;
  Real length_histogram_max;
#ifdef OMEGA_H_USE_EGADS
  Egads* egads_model;
  bool should_smooth_snap;
  Real snap_smooth_tolerance;
#endif
  Int max_motion_steps;
  Real motion_step_size;
  bool should_refine;
  bool should_coarsen;
  bool should_swap;
  bool should_coarsen_slivers;
  bool should_move_for_quality;
};

/* returns false if the mesh was not modified. */
bool adapt(Mesh* mesh, AdaptOpts const& opts);

bool print_adapt_status(Mesh* mesh, AdaptOpts const& opts);
void print_adapt_histograms(Mesh* mesh, AdaptOpts const& opts);

namespace binary {
void write(std::string const& path, Mesh* mesh);
I32 read(std::string const& path, CommPtr comm, Mesh* mesh);
I32 read_nparts(std::string const& path, CommPtr comm);
I32 read_version(std::string const& path, CommPtr comm);
void read_in_comm(
    std::string const& path, CommPtr comm, Mesh* mesh, I32 version);
}

void build_from_elems2verts(
    Mesh* mesh, CommPtr comm, Int edim, LOs ev2v, Read<GO> vert_globals);
void build_from_elems2verts(Mesh* mesh, Int edim, LOs ev2v, LO nverts);
void build_from_elems_and_coords(Mesh* mesh, Int edim, LOs ev2v, Reals coords);
void build_box(Mesh* mesh, Real x, Real y, Real z, LO nx, LO ny, LO nz);

void classify_by_angles(Mesh* mesh, Real sharp_angle);

Adj reflect_down(LOs hv2v, LOs lv2v, Adj v2l, Int high_dim, Int low_dim);

Remotes owners_from_globals(
    CommPtr comm, Read<GO> globals, Read<I32> own_ranks);

Real repro_sum(Reals a);
Real repro_sum(CommPtr comm, Reals a);
void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]);
Real repro_sum_owned(Mesh* mesh, Int dim, Reals a);

OMEGA_H_INLINE bool code_is_flipped(I8 code) { return code & 1; }

OMEGA_H_INLINE Int code_rotation(I8 code) { return (code >> 1) & 3; }

OMEGA_H_INLINE Int code_which_down(I8 code) { return (code >> 3); }

Read<I8> mark_class_closure(
    Mesh* mesh, Int ent_dim, Int class_dim, I32 class_id);
Read<I8> mark_class_closures(Mesh* mesh, Int ent_dim,
    std::vector<Int> class_dims, std::vector<I32> class_ids);
void fix_momentum_velocity_verts(
    Mesh* mesh, Int class_dim, I32 class_id, Int comp);

template <typename T>
Read<I8> each_eq_to(Read<T> a, T b);
template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a);
template <typename T>
Read<T> min_each(Read<T> a, Read<T> b);
LOs collect_marked(Read<I8> marks);

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts);
bool approach_size_field(Mesh* mesh, AdaptOpts const& opts);

Reals find_implied_size(Mesh* mesh);
Reals find_implied_metric(Mesh* mesh);
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
Reals limit_size_field_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol = 1e-3);
Reals expected_elems_per_elem_iso(Mesh* mesh, Reals v2h);
Reals expected_elems_per_elem_metric(Mesh* mesh, Reals v2m);
Real size_scalar_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems);
Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m);
Reals smooth_isos_once(Mesh* mesh, Reals v2h);
Reals get_curvature_isos(Mesh* mesh, Real segment_angle, Real max_size);

Reals recover_hessians(Mesh* mesh, Reals vert_values);
Reals metric_from_hessians(Int dim, Reals hessians, Real eps, Real hmax);
Reals metric_for_nelems_from_hessians(
    Mesh* mesh, Real target_nelems, Real tolerance, Reals hessians, Real hmax);

template <typename T, Int n>
class Few {
  using UninitT = typename std::aligned_storage<sizeof(T), alignof(T)>::type;
  UninitT array_[n];

 public:
  enum { size = n };
  OMEGA_H_INLINE T* data() { return reinterpret_cast<T*>(array_); }
  OMEGA_H_INLINE T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
  OMEGA_H_INLINE T volatile* data() volatile {
    return reinterpret_cast<T volatile*>(array_);
  }
  OMEGA_H_INLINE T const volatile* data() const volatile {
    return reinterpret_cast<T const volatile*>(array_);
  }
  OMEGA_H_INLINE T& operator[](Int i) { return data()[i]; }
  OMEGA_H_INLINE T const& operator[](Int i) const { return data()[i]; }
  OMEGA_H_INLINE T volatile& operator[](Int i) volatile { return data()[i]; }
  OMEGA_H_INLINE T const volatile& operator[](Int i) const volatile {
    return data()[i];
  }
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (data() + (i++)) T(*it);
    }
  }
  OMEGA_H_INLINE Few() {
    for (Int i = 0; i < n; ++i) new (data() + i) T();
  }
  OMEGA_H_INLINE ~Few() {
    for (Int i = 0; i < n; ++i) (data()[i]).~T();
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  OMEGA_H_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
};

template <typename T>
OMEGA_H_INLINE T max2(T a, T b) {
  return (b > a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE T min2(T a, T b) {
  return (b < a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE void swap2(T& a, T& b) {
  T c = a;
  a = b;
  b = c;
}

bool ends_with(std::string const& s, std::string const& suffix);

/* begin explicit instantiation declarations */
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Read<T> permute(Read<T> a_data, LOs a2b, Int width);         \
  extern template Read<I8> each_eq_to(Read<T> a, T b);                         \
  extern template Read<T> multiply_each_by(T factor, Read<T> x);               \
  extern template Read<T> min_each(Read<T> a, Read<T> b);                      \
  extern template T Comm::allreduce(T x, Omega_h_Op op) const;                 \
  extern template T Comm::exscan(T x, Omega_h_Op op) const;                    \
  extern template void Comm::bcast(T& x) const;                                \
  extern template Read<T> Comm::allgather(T x) const;                          \
  extern template Read<T> Comm::alltoall(Read<T> x) const;                     \
  extern template Read<T> Comm::alltoallv(Read<T> sendbuf,                     \
      Read<LO> sendcounts, Read<LO> sdispls, Read<LO> recvcounts,              \
      Read<LO> rdispls) const;                                                 \
  extern template Read<T> Dist::exch(Read<T> data, Int width) const;           \
  extern template Read<T> Dist::exch_reduce<T>(                                \
      Read<T> data, Int width, Omega_h_Op op) const;                           \
  extern template Tag<T> const* Mesh::get_tag<T>(                              \
      Int dim, std::string const& name) const;                                 \
  extern template Read<T> Mesh::get_array<T>(Int dim, std::string const& name) \
      const;                                                                   \
  extern template void Mesh::add_tag<T>(                                       \
      Int dim, std::string const& name, Int ncomps, Int xfer, Int outflags);   \
  extern template void Mesh::add_tag<T>(Int dim, std::string const& name,      \
      Int ncomps, Int xfer, Int outflags, Read<T> array, bool internal);       \
  extern template void Mesh::set_tag(                                          \
      Int dim, std::string const& name, Read<T> array, bool internal);         \
  extern template Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width); \
  extern template Read<T> Mesh::owned_array(                                   \
      Int ent_dim, Read<T> a, Int width);                                      \
  extern template Read<T> Mesh::sync_subset_array(                             \
      Int ent_dim, Read<T> a_data, LOs a2e, T default_val, Int width);         \
  extern template Read<T> Mesh::reduce_array(                                  \
      Int ent_dim, Read<T> a, Int width, Omega_h_Op op);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL
/* end explicit instantiation declarations */

}  // end namespace Omega_h

#endif
