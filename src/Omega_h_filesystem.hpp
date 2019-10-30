#ifndef OMEGA_H_FILESYSTEM_HPP
#define OMEGA_H_FILESYSTEM_HPP

#include <cstdint>
#include <iosfwd>
#include <string>
#include <system_error>

// our own tiny subset of std::filesystem while we wait for C++17

namespace Omega_h {

namespace filesystem {

class filesystem_error : public std::system_error {
 public:
  filesystem_error(filesystem_error const&) = default;
  virtual ~filesystem_error() override;
  filesystem_error(int ev, const char* what_arg);
};

class path {
 public:
  using value_type = char;
#ifdef _MSC_VER
  static constexpr value_type preferred_separator = '\\';
#else
  static constexpr value_type preferred_separator = '/';
#endif
  using string_type = std::basic_string<value_type>;
  path() = default;
  path(value_type const* source);
  path(string_type const& source);
  value_type const* c_str() const noexcept;
  const string_type& native() const noexcept;
  std::string string() const;
  path filename() const;
  path extension() const;
  path stem() const;
  path parent_path() const;
  path& operator/=(path const&);
  path& operator/=(std::string const&);
  path& operator/=(char const*);
  path& operator+=(path const&);
  path& operator+=(std::string const&);
  path& operator+=(char const*);
  path& operator+=(char);

 public:
  string_type impl;
};

path operator/(path const& a, path const& b);
std::ostream& operator<<(std::ostream& os, path const& p);

bool create_directory(path const& p);
path current_path();
bool remove(path const& p);
std::uintmax_t remove_all(path const& p);
bool exists(path const& p);

enum class file_type {
  none,
  not_found,
  regular,
  directory,
  symlink,
  block,
  character,
  fifo,
  socket,
  unknown
};

class file_status {
 public:
  file_status(file_type const& type_in);
  file_type type() const noexcept;

 private:
  file_type type_variable;
};

file_status status(path const& p);

class directory_entry {
 public:
  directory_entry() = default;
  directory_entry(filesystem::path const& p);
  const filesystem::path& path() const noexcept { return path_variable; }
  bool is_regular_file() const;
  bool is_directory() const;
  bool is_symlink() const;

 private:
  filesystem::path path_variable;
};

struct IteratorImpl;

class directory_iterator {
 public:
  directory_iterator();
  directory_iterator(directory_iterator&& other);
  directory_iterator(path const& p);
  ~directory_iterator();
  directory_iterator& operator++();
  const directory_entry& operator*() const;
  const directory_entry* operator->() const;
  bool operator==(directory_iterator const& other) const;
  // hassle to implement, not needed
  directory_iterator(directory_iterator const& other) = delete;

 private:
  IteratorImpl* impl;
  directory_entry entry;
};

inline bool operator!=(
    directory_iterator const& a, directory_iterator const& b) {
  return !(a == b);
}

}  // namespace filesystem

}  // namespace Omega_h

#endif
