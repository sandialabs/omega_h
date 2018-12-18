#ifndef OMEGA_H_MALLOC_HPP
#define OMEGA_H_MALLOC_HPP

#include <cstddef>

namespace Omega_h {

void* device_malloc(std::size_t size);
void device_free(void* ptr, std::size_t size);
void* host_malloc(std::size_t size);
void host_free(void* ptr, std::size_t size);

void enable_pooling();
void disable_pooling();

void* maybe_pooled_device_malloc(std::size_t size);
void maybe_pooled_device_free(void* ptr, std::size_t size);
void* maybe_pooled_host_malloc(std::size_t size);
void maybe_pooled_host_free(void* ptr, std::size_t size);
}  // namespace Omega_h

#endif
