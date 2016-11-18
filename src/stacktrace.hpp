#ifndef STACKTRACE_HPP
#define STACKTRACE_HPP

#include <iosfwd>

namespace Omega_h {

void print_stacktrace(std::ostream& out, int max_frames);

}

#endif
