// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

/* Dan Ibanez: found this code for stacktrace printing,
   cleaned it up for compilation as strict C++11 */

#ifndef STACKTRACE_HPP
#define STACKTRACE_HPP

#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <cxxabi.h>

/** Print a demangled stack backtrace of the caller function to FILE* out. */
static inline void print_stacktrace(FILE *out = stderr, int max_frames = 63)
{
  fprintf(out, "stack trace:\n");
  // storage array for stack trace address data
  void** addrlist = new void*[max_frames];
  // retrieve current stack addresses
  int addrlen = backtrace(addrlist, max_frames);
  if (addrlen == 0) {
    fprintf(out, "  <empty, possibly corrupt>\n");
    return;
  } else {
    fprintf(out, "addrlen %d\n", addrlen);
  }
  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char** symbollist = backtrace_symbols(addrlist, addrlen);
  delete [] addrlist;
  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = new char[funcnamesize];
  // iterate over the returned symbol lines. skip the first, it is the
  // address of this function.
  for (int i = 1; i < addrlen; i++)
  {
    // TODO: this parsing code assumes certain GNU/Linux formatting.
    // try adding OSX code when time allows
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;
    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = symbollist[i]; *p; ++p)
    {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset) {
        end_offset = p;
        break;
      }
    }
    if (begin_name && begin_offset && end_offset
        && begin_name < begin_offset)
    {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';
      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():
      int status;
      char* ret = abi::__cxa_demangle(begin_name,
          funcname, &funcnamesize, &status);
      if (status == 0) {
        funcname = ret; // use possibly realloc()-ed string
        fprintf(out, "  %s : %s+%s\n",
            symbollist[i], funcname, begin_offset);
      }
      else {
        // demangling failed. Output function name as a C function with
        // no arguments.
        fprintf(out, "  %s : %s()+%s\n",
            symbollist[i], begin_name, begin_offset);
      }
    }
    else
    {
      // couldn't parse the line? print the whole line.
      fprintf(out, " (noparse)  %s\n", symbollist[i]);
    }
  }
  delete [] funcname;
  free(symbollist);
}

#endif // STACKTRACE_H
