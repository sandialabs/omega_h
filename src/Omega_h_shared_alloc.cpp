#include <Omega_h_fail.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_malloc.hpp>
#include <Omega_h_profile.hpp>
#include <Omega_h_shared_alloc.hpp>
#include <sstream>

namespace Omega_h {

OMEGA_H_DLL bool entering_parallel = false;
Allocs* global_allocs = nullptr;

void start_tracking_allocations() {
  OMEGA_H_CHECK(global_allocs == nullptr);
  global_allocs = new Allocs();
  global_allocs->first = nullptr;
  global_allocs->last = nullptr;
  global_allocs->total_bytes = 0;
  global_allocs->high_water_bytes = 0;
}

void stop_tracking_allocations(Library* lib) {
  OMEGA_H_CHECK(global_allocs != nullptr);
  auto comm = lib->world();
  auto mpi_high_water_bytes =
      comm->allreduce(I64(global_allocs->high_water_bytes), OMEGA_H_MAX);
  auto candidate_rank =
      (I64(global_allocs->high_water_bytes) == mpi_high_water_bytes)
          ? comm->rank()
          : comm->size();
  auto winner_rank = comm->allreduce(candidate_rank, OMEGA_H_MIN);
  if (comm->rank() == winner_rank) {
    std::stringstream ss;
    ss << "MPI rank " << comm->rank()
       << " had the highest memory high water mark.\n";
    ss << "The high water memory usage was " << global_allocs->high_water_bytes
       << " bytes, allocated as follows:\n";
    for (auto& rec : global_allocs->high_water_records) {
      ss << rec.name << ": " << rec.bytes << " bytes.\n";
    }
    auto s = ss.str();
    std::printf("%s\n", s.c_str());
  }
  delete global_allocs;
  global_allocs = nullptr;
}

Alloc::Alloc(std::size_t size_in, std::string const& name_in)
    : size(size_in), name(name_in) {
  init();
}

Alloc::Alloc(std::size_t size_in, std::string&& name_in)
    : size(size_in), name(name_in) {
  init();
}

OMEGA_H_DLL Alloc::~Alloc() {
  ::Omega_h::maybe_pooled_device_free(ptr, size);
  auto ga = global_allocs;
  if (ga) {
    if (next == nullptr) {
      ga->last = prev;
    } else {
      next->prev = prev;
    }
    if (prev == nullptr) {
      ga->first = next;
    } else {
      prev->next = next;
    }
    ga->total_bytes -= size;
  }
}

void Alloc::init() {
  ptr = ::Omega_h::maybe_pooled_device_malloc(size);
  use_count = 1;
  auto ga = global_allocs;
  if (size && (ptr == nullptr)) {
    std::stringstream ss;
    ss << "Failed to allocate " << size << " bytes for " << name << '\n';
    if (ga) {
      ss << "at the time, " << ga->total_bytes
         << " total bytes were allocated as follows:\n";
      for (auto a = ga->first; a; a = a->next) {
        ss << a->name << ": " << a->size << " bytes\n";
      }
    }
    auto s = ss.str();
    Omega_h_fail("%s\n", s.c_str());
  }
  if (ga) {
    auto old_last = ga->last;
    this->prev = old_last;
    this->next = nullptr;
    if (old_last) {
      old_last->next = this;
    } else {
      ga->first = this;
      ga->last = this;
    }
    ga->total_bytes += size;
    if (ga->total_bytes > ga->high_water_bytes) {
      Omega_h::ScopedTimer high_water_timer("high water update");
      ga->high_water_bytes = ga->total_bytes;
      ga->high_water_records.clear();
      for (auto a = ga->first; a; a = a->next) {
        ga->high_water_records.push_back({a->name, a->size});
      }
    }
  }
}

SharedAlloc::SharedAlloc(std::size_t size_in, std::string const& name_in) {
  alloc = new Alloc(size_in, name_in);
  direct_ptr = alloc->ptr;
}

SharedAlloc::SharedAlloc(std::size_t size_in, std::string&& name_in) {
  alloc = new Alloc(size_in, name_in);
  direct_ptr = alloc->ptr;
}

SharedAlloc::SharedAlloc(std::size_t size_in) : SharedAlloc(size_in, "") {}

SharedAlloc SharedAlloc::identity(std::size_t size_in) {
  SharedAlloc out;
  out.direct_ptr = nullptr;
  out.alloc = reinterpret_cast<Alloc*>(
      (static_cast<std::uintptr_t>(size_in) << 3) & IS_IDENTITY);
  return out;
}

}  // namespace Omega_h
