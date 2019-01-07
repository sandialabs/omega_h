#include <Omega_h_filesystem.hpp>
#include <cerrno>
#include <cstdio>
#include <vector>
#include <cstring>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

namespace Omega_h {

namespace filesystem {

bool create_directory(path const& p, std::error_code& ec) {
  ::mode_t const mode = S_IRWXU | S_IRWXG | S_IRWXO;
  ::errno = 0;
  int err = ::mkdir(path, mode);
  if (err != 0) {
    if (::errno != EEXIST) {
      throw filesystem_error(::errno, "Omega_h::filesystem::create_directory");
    }
    return false;
  }
  return true;
}

path current_path(std::error_code& ec) {
  char buf[1024];
  ::errno = 0;
  char* ret = ::getcwd(buf, sizeof(buf));
  if (ret == nullptr) {
    throw filesystem_error(::errno, "Omega_h::filesystem::current_path");
  }
  return ret;
}

bool remove(path const& p) {
  ::errno = 0;
  int err = ::remove(p.impl.c_str());
  if (err != 0) {
    sys_str = std::strerror(::errno);
    throw filesystem_error(::errno, "Omega_h::filesystem::remove");
  }
}

std::uintmax_t remove_all(path const& p) {
  recursive_directory_iterator it(p);
}

bool exists(path const& p) {
  return ::access(p.impl.c_str(), F_OK) != -1;
}

struct IteratorImpl {
  IteratorImpl():stream(nullptr),entry(nullptr) {}
  IteratorImpl(path const& p):root(p),entry(nullptr) {
    ::errno = 0;
    stream = ::opendir(p.impl.c_str());
    if (stream == nullptr) {
      throw filesystem_error(::errno, "Omega_h::filesystem::directory_iterator");
    }
    increment();
  }
  ~IteratorImpl() {
    close();
  }
  void close() {
    if (stream == nullptr) return;
    ::errno = 0;
    int err = ::closedir(stream);
    stream = nullptr;
    if (err != 0) {
      throw filesystem_error(::errno, "Omega_h::filesystem::directory_iterator");
    }
  }
  void increment() {
    ::errno = 0;
    entry = ::readdir(stream);
    if (entry == nullptr) {
      if (::errno != 0) {
        throw filesystem_error(::errno, "Omega_h::filesystem::directory_iterator");
      }
      // safely reached the end of the directory
      close();
    } else {
      // we just extracted a good entry from the stream
      // skip dot and dot-dot, max 3-call recursion
      if (0 == strcmp(entry->d_name, ".")) increment();
      else if (0 == strcmp(entry->d_name, "..")) increment();
    }
  }
  directory_entry deref() {
    return directory_entry(root / entry->d_name);
  }
  bool is_end() { return entry == nullptr; }
  bool equal(IteratorImpl const& other) {
    if (entry == nullptr && other.entry == nullptr) return true;
    if (root.impl != other.root.impl) return false;
    return 0 == strcmp(entry->d_name, other.entry->d_name);
  }
  path root;
  DIR* stream;
  ::dirent* entry;
};

file_status status(path const& p) {
  ::errno = 0;
  ::stat buf;
  int ret = ::stat(p.c_str(), &buf);
  if (ret != 0) {
    throw filesystem_error(::errno, "Omega_h::filesystem::status");
  }
  file_type type;
  if (buf.st_mode & S_IFREG) type = file_type::regular;
  else if (buf.st_mode & S_IFDIR) type = file_type::directory;
  else if (buf.st_mode & S_IFLNK) type = file_type::symlink;
  else if (buf.st_mode & S_IFBLK) type = file_type::block;
  else if (buf.st_mode & S_IFCHR) type = file_type::character;
  else if (buf.st_mode & S_IFIFO) type = file_type::fifo;
  return file_status(type);
}

// end of OS-specific stuff

path::path(char const* source):impl(source) {
}

path::path(std::string const& source):impl(source) {
}

value_type const* path::c_str() const noexcept {
  return native().c_str();
}

const string_type& path::native() const noexcept {
  return impl;
}

filesystem_error::filesystem_error(int ev, const char* what_arg)
  : std::system_error(ev, std::system_category(), what_arg)
{
}

file_status::file_status(file_type const& type_in)
  :type_variable(type_in)
{
}

file_type file_status::type() const noexcept {
  return type_variable;
}

directory_entry::directory_entry(path const& p):path_variable(p) {}

const filesystem::path& directory_entry::path() const noexcept {
  return path_variable;
}

bool directory_entry::is_regular_file() const noexcept {
  return status(path_variable).type() == file_type::regular;
}

bool directory_entry::is_directory() const noexcept {
  return status(path_variable).type() == file_type::directory;
}

bool directory_entry::is_symlink() const noexcept {
  return status(path_variable).type() == file_type::symlink;
}

directory_iterator::directory_iterator()
 : impl(new IteratorImpl())
{
}

directory_iterator::directory_iterator(directory_iterator&& other) {
  delete impl;
  impl = other.impl;
  other.impl = new IteratorImpl();
}

directory_iterator::directory_iterator(path const& p)
 : impl(new IteratorImpl(p))
 , entry(impl->deref())
{
}

directory_iterator::~directory_iterator() {
  delete impl;
  impl = nullptr;
}

directory_iterator& directory_iterator::operator++() {
  impl->increment();
  if (!impl->is_ent()) entry = impl->deref();
  return *this;
}

const directory_entry& operator*() const {
  return entry;
}

const directory_entry* operator-> const {
  return &entry;
}

bool operator!=(directory_iterator const& other) {
  return !impl->equal(other.impl);
}

}

}
