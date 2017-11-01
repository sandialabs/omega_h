#ifndef CONTROL_HPP
#define CONTROL_HPP

#include <string>

namespace Omega_h {
extern bool should_log_memory;
extern char* max_memory_stacktrace;
void print_stacktrace(std::ostream& out, int max_frames);
}  // namespace Omega_h

#endif
