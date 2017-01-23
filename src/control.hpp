#ifndef CONTROL_HPP
#define CONTROL_HPP

#include <string>

namespace Omega_h {
extern bool should_log_memory;
extern char* max_memory_stacktrace;
void print_stacktrace(std::ostream& out, int max_frames);
void add_to_global_timer(std::string const& name, double nsecs);
}

#endif
