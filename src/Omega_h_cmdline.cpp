#include "Omega_h_cmdline.hpp"
#include "internal.hpp"

#include <iostream>

namespace Omega_h {

template <typename T>
struct ItemTraits;

template <>
struct ItemTraits<std::string> {
  static char const* name() { return "std::string"; }
};

template <>
struct ItemTraits<double> {
  static char const* name() { return "double"; }
};

template <>
struct ItemTraits<int> {
  static char const* name() { return "int"; }
};

static void shift(int* p_argc, char** argv, int i) {
  --(*p_argc);
  for (auto j = i; j < *p_argc; ++j) argv[j] = argv[j + 1];
}

static bool parse_arg(char const* arg, std::string* p_value) {
  *p_value = arg;
  return true;
}

static bool parse_arg(char const* arg, int* p_value) {
  *p_value = atoi(arg);
  return true;
}

static bool parse_arg(char const* arg, double* p_value) {
  *p_value = atof(arg);
  return true;
}

CmdLineItem::~CmdLineItem() {}

CmdLineItem::CmdLineItem(std::string const& name)
    : name_(name), parsed_(false) {}

bool CmdLineItem::parse(CommPtr comm, int* p_argc, char** argv, int i) {
  if (parsed_) {
    if (!comm->rank())
      std::cout << "argument <" << name() << " being parsed twice!\n";
    return false;
  }
  auto ok = parse_impl(comm, p_argc, argv, i);
  if (ok) parsed_ = true;
  return ok;
}

std::string const& CmdLineItem::name() const { return name_; }

bool CmdLineItem::parsed() const { return parsed_; }

template <typename T>
CmdLineArg<T>::CmdLineArg(std::string const& name, T const& defval)
    : CmdLineItem(name), value_(defval) {}

template <typename T>
CmdLineArg<T>::~CmdLineArg() {}

template <typename T>
bool CmdLineArg<T>::parse_impl(CommPtr comm, int* p_argc, char** argv, int i) {
  if (!parse_arg(argv[i], &value_)) {
    if (!comm->rank()) {
      std::cout << "could not parse \"" << argv[i] << "\" as type "
                << ItemTraits<T>::name() << '\n';
    }
    return false;
  }
  shift(p_argc, argv, i);
  return true;
}

template <typename T>
T CmdLineArg<T>::get() const {
  return value_;
}

CmdLineFlag::CmdLineFlag(std::string const& name, std::string const& desc)
    : CmdLineItem(name), desc_(desc) {}

bool CmdLineFlag::parse_impl(CommPtr comm, int* p_argc, char** argv, int i) {
  shift(p_argc, argv, i);
  if (((*p_argc) - i) < int(args_.size())) {
    if (!comm->rank()) {
      std::cout << "flag " << name() << " takes " << args_.size()
                << " arguments\n";
    }
  }
  for (auto const& arg : args_) {
    if (!arg->parse(comm, p_argc, argv, i)) {
      if (!comm->rank()) {
        std::cout << "could not parse argument <" << arg->name() << "> of flag "
                  << name() << '\n';
      }
      return false;
    }
  }
  return true;
}

template <typename T>
void CmdLineFlag::add_arg(std::string const& name, T const& defval) {
  args_.push_back(
      std::unique_ptr<CmdLineArg<T>>(new CmdLineArg<T>(name, defval)));
}

std::string const& CmdLineFlag::desc() const { return desc_; }

CmdLineItem* CmdLineFlag::arg(std::size_t i) { return args_.at(i).get(); }

CmdLineItem* CmdLineFlag::arg(std::string const& arg_name) {
  for (auto const& arg : args_)
    if (arg->name() == arg_name) return arg.get();
  NORETURN(nullptr);
}

std::size_t CmdLineFlag::nargs() const { return args_.size(); }

CmdLine::CmdLine() : nargs_parsed_(0) {}

bool CmdLine::parse(CommPtr comm, int* p_argc, char** argv) {
  for (int i = 1; i < *p_argc;) {
    bool parsed_by_flag = false;
    for (auto const& flag : flags_) {
      if (flag->name() == argv[i]) {
        if (!flag->parse(comm, p_argc, argv, i)) {
          if (!comm->rank()) {
            std::cout << "failed to parse flag " << flag->name() << '\n';
          }
          return false;
        }
        parsed_by_flag = true;
      }
    }
    if (parsed_by_flag) {
      continue;
    } else if (std::string("-help") == argv[i]) {
      return false;
    } else if (std::string("--help") == argv[i]) {
      return false;
    } else if (args_.size() > nargs_parsed_) {
      if (args_[nargs_parsed_]->parse(comm, p_argc, argv, i)) {
        ++nargs_parsed_;
      }
    } else {
      ++i;
    }
  }
  if (nargs_parsed_ < args_.size()) {
    if (!comm->rank()) {
      for (std::size_t i = nargs_parsed_; i < args_.size(); ++i) {
        std::cout << "missing required argument <" << args_[i]->name() << ">\n";
      }
    }
    return false;
  }
  return true;
}

CmdLineFlag& CmdLine::add_flag(
    std::string const& name, std::string const& desc) {
  flags_.push_back(std::unique_ptr<CmdLineFlag>(new CmdLineFlag(name, desc)));
  return *(flags_.back());
}

template <typename T>
void CmdLine::add_arg(std::string const& name, T const& defval) {
  args_.push_back(
      std::unique_ptr<CmdLineArg<T>>(new CmdLineArg<T>(name, defval)));
}

bool CmdLine::parsed(std::string const& flag_name) const {
  for (auto const& flag : flags_)
    if (flag->name() == flag_name) return flag->parsed();
  return false;
}

template <typename T>
static T get(CmdLineItem const* p) {
  auto p2 = dynamic_cast<CmdLineArg<T> const*>(p);
  CHECK(p2);
  return p2->get();
}

template <typename T>
T CmdLine::get(std::string const& flag_name, std::string const& arg_name) const {
  for (auto const& flag : flags_)
    if (flag->name() == flag_name) return Omega_h::get<T>(flag->arg(arg_name));
  NORETURN(T());
}

bool CmdLine::parsed(std::string const& flag_name, std::size_t i) const {
  for (auto const& flag : flags_)
    if (flag->name() == flag_name) return flag->arg(i)->parsed();
  return false;
}

template <typename T>
T CmdLine::get(std::string const& arg_name) const {
  for (auto const& arg : args_)
    if (arg->name() == arg_name) return Omega_h::get<T>(arg.get());
  NORETURN(T());
}

bool CmdLine::parsed(std::size_t i) const { return args_.at(i)->parsed(); }

bool CmdLine::check_empty(CommPtr comm, int argc, char** argv) {
  if (argc == 1) return true;
  if (!comm->rank()) {
    for (int i = 1; i < argc; ++i) {
      std::cout << "unknown argument \"" << argv[i] << "\"\n";
    }
  }
  return false;
}

void CmdLine::show_help(CommPtr comm, char** argv) const {
  if (comm->rank()) return;
  std::cout << "usage: " << argv[0];
  if (!flags_.empty()) std::cout << " [options]";
  for (auto const& arg : args_) std::cout << " " << arg->name();
  std::cout << '\n';
  if (!flags_.empty()) {
    std::cout << "options:\n";
    for (auto const& flag : flags_) {
      std::cout << "  " << flag->name();
      for (std::size_t i = 0; i < flag->nargs(); ++i) {
        std::cout << " " << flag->arg(i)->name();
      }
      std::cout << '\n';
      std::cout << "    " << flag->desc() << '\n';
    }
  }
}

#define OMEGA_H_EXPL_INST(T)                                                   \
  template class CmdLineArg<T>;                                                \
  template void CmdLineFlag::add_arg<T>(                                       \
      std::string const& name, T const& defval);                               \
  template void CmdLine::add_arg<T>(std::string const& name, T const& defval); \
  template T CmdLine::get<T>(                                                  \
      std::string const& flag_name, std::string const& arg_name) const;              \
  template T CmdLine::get<T>(std::string const& arg_name) const;
OMEGA_H_EXPL_INST(int)
OMEGA_H_EXPL_INST(double)
OMEGA_H_EXPL_INST(std::string)
}
