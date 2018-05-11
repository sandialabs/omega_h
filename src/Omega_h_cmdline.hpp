#ifndef OMEGA_H_CMDLINE_HPP
#define OMEGA_H_CMDLINE_HPP

#include <memory>
#include <string>
#include <vector>

#include <Omega_h_comm.hpp>

namespace Omega_h {

class CmdLineItem {
 public:
  CmdLineItem(std::string const& name_in);
  virtual ~CmdLineItem();
  bool parse(int* p_argc, char** argv, int i, bool should_print);
  virtual bool parse_impl(
      int* p_argc, char** argv, int i, bool should_print) = 0;
  std::string const& name() const;
  bool parsed() const;

 protected:
  std::string name_;
  bool parsed_;
};

template <typename T>
class CmdLineArg : public CmdLineItem {
 public:
  CmdLineArg(std::string const& name_in, T const& defval);
  ~CmdLineArg() override;
  bool parse_impl(int* p_argc, char** argv, int i, bool should_print) override;
  T get() const;

 private:
  T value_;
};

class CmdLineFlag : public CmdLineItem {
 public:
  CmdLineFlag(std::string const& name_in, std::string const& desc_in);
  virtual bool parse_impl(
      int* p_argc, char** argv, int i, bool should_print) override;
  template <typename T>
  void add_arg(std::string const& arg_name, T const& defval = T());
  std::string const& desc() const;
  CmdLineItem* arg(std::size_t i);
  CmdLineItem* arg(std::string const& arg_name);
  std::size_t nargs() const;

 private:
  CmdLineFlag(CmdLineFlag const& other) = delete;
  std::string desc_;
  std::vector<std::unique_ptr<CmdLineItem>> args_;
};

class CmdLine {
 public:
  CmdLine();
  OMEGA_H_NODISCARD
  bool parse_final(CommPtr comm, int* p_argc, char** argv);
  OMEGA_H_NODISCARD
  bool parse(CommPtr comm, int* p_argc, char** argv);
  CmdLineFlag& add_flag(std::string const& name, std::string const& desc);
  template <typename T>
  void add_arg(std::string const& name, T const& defval = T());
  bool parsed(std::string const& flag_name) const;
  template <typename T>
  T get(std::string const& flag_name, std::string const& arg_name) const;
  bool parsed(std::string const& flag_name, std::size_t i) const;
  template <typename T>
  T get(std::string const& arg_name) const;
  bool parsed(std::size_t i) const;
  static bool check_empty(CommPtr comm, int argc, char** argv);
  static bool check_empty(int argc, char** argv, bool should_print);
  void show_help(CommPtr comm, char** argv) const;
  void show_help(char** argv) const;

 private:
  bool parse(int* p_argc, char** argv, bool should_print);
  std::vector<std::unique_ptr<CmdLineItem>> args_;
  std::vector<std::unique_ptr<CmdLineFlag>> flags_;
  std::size_t nargs_parsed_;
  bool parsed_help_;
};

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template class CmdLineArg<T>;                                         \
  extern template void CmdLineFlag::add_arg<T>(                                \
      std::string const& name, T const& defval);                               \
  extern template void CmdLine::add_arg<T>(                                    \
      std::string const& name, T const& defval);                               \
  extern template T CmdLine::get<T>(                                           \
      std::string const& flag_name, std::string const& arg_name) const;        \
  extern template T CmdLine::get<T>(std::string const& arg_name) const;
OMEGA_H_EXPL_INST_DECL(int)
OMEGA_H_EXPL_INST_DECL(double)
OMEGA_H_EXPL_INST_DECL(std::string)
#undef OMEGA_H_EXPL_INST_DECL

}  // end namespace Omega_h

#endif
