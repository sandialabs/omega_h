#ifndef OMEGA_H_CMDLINE_HPP
#define OMEGA_H_CMDLINE_HPP

#include "Omega_h.hpp"

namespace Omega_h {

class CmdLineItem {
  public:
    CmdLineItem(std::string const& name);
    virtual ~CmdLineItem();
    bool parse(CommPtr comm, int* p_argc, char** argv, int i);
    virtual bool parse_impl(CommPtr comm, int* p_argc, char** argv, int i) = 0;
    std::string const& name() const;
    bool parsed() const;
  protected:
    std::string name_;
    bool parsed_;
};

template <typename T>
class CmdLineArg : public CmdLineItem {
  public:
    CmdLineArg(std::string const& name, T const& defval);
    virtual ~CmdLineArg();
    virtual bool parse_impl(CommPtr comm, int* p_argc, char** argv, int i);
    T get() const;
  private:
    T value_;
};

class CmdLineFlag : public CmdLineItem {
  public:
    CmdLineFlag(std::string const& name, std::string const& desc);
    virtual bool parse_impl(CommPtr comm, int* p_argc, char** argv, int i);
    template <typename T>
    void add_arg(std::string const& name, T const& defval = T());
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
    bool parse(CommPtr comm, int* p_argc, char** argv);
    CmdLineFlag& add_flag(std::string const& name, std::string const& desc);
    template <typename T>
    void add_arg(std::string const& name, T const& defval = T());
    bool parsed(std::string const& flag_name);
    template <typename T>
    T get(std::string const& flag_name, std::string const& arg_name);
    bool parsed(std::string const& flag_name, std::size_t i);
    template <typename T>
    T get(std::string const& arg_name);
    bool parsed(std::size_t i);
    static bool check_empty(CommPtr comm, int argc, char** argv);
    void show_help(CommPtr comm, char** argv);
  private:
    std::vector<std::unique_ptr<CmdLineItem>> args_;
    std::vector<std::unique_ptr<CmdLineFlag>> flags_;
    std::size_t nargs_parsed_;
};

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template class CmdLineArg<T>;                                         \
  extern template \
  void CmdLineFlag::add_arg<T>(std::string const& name, T const& defval); \
  extern template \
  void CmdLine::add_arg<T>(std::string const& name, T const& defval); \
  extern template \
  T CmdLine::get<T>(std::string const& flag_name, std::string const& arg_name); \
  extern template \
  T CmdLine::get<T>(std::string const& arg_name);
OMEGA_H_EXPL_INST_DECL(int)
OMEGA_H_EXPL_INST_DECL(double)
OMEGA_H_EXPL_INST_DECL(std::string)

}

#endif
