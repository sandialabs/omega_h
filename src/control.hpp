#ifndef CONTROL_HPP
#define CONTROL_HPP

#include <string>

namespace Omega_h {
extern bool should_log_memory;
extern char* max_memory_stacktrace;
void add_to_global_timer(std::string const& name, double nsecs);
}

#endif
