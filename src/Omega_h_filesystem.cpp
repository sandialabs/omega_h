#include <Omega_h_fail.hpp>
#include <Omega_h_filesystem.hpp>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace Omega_h {

namespace filesystem {

// begin OS-specific stuff

#ifdef _MSC_VER

bool create_directory(path const& p) {
  BOOL success = ::CreateDirectoryA(p.impl.c_str(), nullptr);
  if (!success) {
    if (GetLastError() != ERROR_ALREADY_EXISTS) {
      throw filesystem_error(GetLastError(), "create_directory");
    }
    return false;
  }
  return true;
}

path current_path() {
  char buf[1024];
  static_assert(std::is_same<TCHAR, char>::value, "expecting TCHAR to be char");
  int ret = ::GetCurrentDirectoryA(sizeof(buf), buf);
  if (ret == 0) {
    throw filesystem_error(GetLastError(), "current_path");
  }
  return std::string(buf);
}

bool remove(path const& p) {
  BOOL success = ::RemoveDirectoryA(p.impl.c_str());
  if (!success) {
    throw filesystem_error(GetLastError(), "remove");
  }
  return true;
}

bool exists(path const& p) {
  DWORD dwAttrib = ::GetFileAttributesA(p.impl.c_str());
  return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
          !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}

struct IteratorImpl {
  IteratorImpl() : stream(nullptr), entry_is_good(false) {}
  IteratorImpl(path const& p) : root(p), entry_is_good(false) {
    std::string search_string = p.impl;
    search_string += "\\*";
    stream = ::FindFirstFileA(search_string.c_str(), &entry);
    if (stream == INVALID_HANDLE_VALUE) {
      throw filesystem_error(GetLastError(), "directory_iterator");
    }
    entry_is_good = true;
    if ((0 == strcmp(entry.cFileName, ".")) || (0 == strcmp(entry.cFileName, ".."))) {
      increment();
    }
  }
  ~IteratorImpl() { close(); }
  void close() {
    if (stream == nullptr) return;
    BOOL success = ::FindClose(stream);
    stream = nullptr;
    if (!success) {
      throw filesystem_error(GetLastError(), "directory_iterator");
    }
  }
  void increment() {
    entry_is_good = ::FindNextFileA(stream, &entry);
    if (!entry_is_good) {
      if (GetLastError() != ERROR_NO_MORE_FILES) {
        throw filesystem_error(GetLastError(), "directory_iterator");
      }
      // safely reached the end of the directory
      close();
    } else {
      if ((0 == strcmp(entry.cFileName, ".")) ||
          (0 == strcmp(entry.cFileName, ".."))) {
        increment();
      }
    }
  }
  directory_entry deref() {
    OMEGA_H_CHECK(entry_is_good);
    return directory_entry(root / entry.cFileName);
  }
  bool is_end() { return !entry_is_good; }
  bool equal(IteratorImpl const& other) const {
    if (!entry_is_good && !other.entry_is_good) return true;
    if (root.impl != other.root.impl) return false;
    return 0 == strcmp(entry.cFileName, other.entry.cFileName);
  }
  path root;
  HANDLE stream;
  WIN32_FIND_DATAA entry;
  BOOL entry_is_good;
};

file_status status(path const& p) {
  DWORD dwAttrib = ::GetFileAttributesA(p.impl.c_str());
  if (dwAttrib == INVALID_FILE_ATTRIBUTES) {
    throw filesystem_error(GetLastError(), "status");
  }
  file_type type;
  if (dwAttrib & FILE_ATTRIBUTE_DIRECTORY)
    type = file_type::directory;
  else
    type = file_type::regular;
  return file_status(type);
}

#else

bool create_directory(path const& p) {
  ::mode_t const mode = S_IRWXU | S_IRWXG | S_IRWXO;
  errno = 0;
  int err = ::mkdir(p.c_str(), mode);
  if (err != 0) {
    if (errno != EEXIST) {
      throw filesystem_error(errno, "create_directory");
    }
    return false;
  }
  return true;
}

path current_path() {
  char buf[1024];
  errno = 0;
  char* ret = ::getcwd(buf, sizeof(buf));
  if (ret == nullptr) {
    throw filesystem_error(errno, "current_path");
  }
  return ret;
}

bool remove(path const& p) {
  errno = 0;
  int err = ::remove(p.impl.c_str());
  if (err != 0) {
    throw filesystem_error(errno, "remove");
  }
  return true;
}

bool exists(path const& p) { return ::access(p.impl.c_str(), F_OK) != -1; }

struct IteratorImpl {
  IteratorImpl() : stream(nullptr), entry(nullptr) {}
  IteratorImpl(path const& p) : root(p), entry(nullptr) {
    errno = 0;
    stream = ::opendir(p.impl.c_str());
    if (stream == nullptr) {
      throw filesystem_error(errno, "directory_iterator");
    }
    increment();
  }
  ~IteratorImpl() { close(); }
  void close() {
    if (stream == nullptr) return;
    errno = 0;
    int err = ::closedir(stream);
    stream = nullptr;
    if (err != 0) {
      throw filesystem_error(errno, "directory_iterator");
    }
  }
  void increment() {
    errno = 0;
    entry = ::readdir(stream);
    if (entry == nullptr) {
      if (errno != 0) {
        throw filesystem_error(errno, "directory_iterator");
      }
      // safely reached the end of the directory
      close();
    } else {
      // we just extracted a good entry from the stream
      // skip dot and dot-dot, max 3-call recursion
      if (0 == strcmp(entry->d_name, "."))
        increment();
      else if (0 == strcmp(entry->d_name, ".."))
        increment();
    }
  }
  directory_entry deref() {
    OMEGA_H_CHECK(entry != nullptr);
    return directory_entry(root / entry->d_name);
  }
  bool is_end() { return entry == nullptr; }
  bool equal(IteratorImpl const& other) const {
    if (entry == nullptr && other.entry == nullptr) return true;
    if (root.impl != other.root.impl) return false;
    return 0 == strcmp(entry->d_name, other.entry->d_name);
  }
  path root;
  DIR* stream;
  ::dirent* entry;
};

file_status status(path const& p) {
  errno = 0;
  struct ::stat buf;
  int ret = ::stat(p.c_str(), &buf);
  if (ret != 0) {
    throw filesystem_error(errno, "status");
  }
  file_type type;
  if (buf.st_mode & S_IFREG)
    type = file_type::regular;
  else if (buf.st_mode & S_IFDIR)
    type = file_type::directory;
  else if (buf.st_mode & S_IFLNK)
    type = file_type::symlink;
  else if (buf.st_mode & S_IFBLK)
    type = file_type::block;
  else if (buf.st_mode & S_IFCHR)
    type = file_type::character;
  else if (buf.st_mode & S_IFIFO)
    type = file_type::fifo;
  else
    throw filesystem_error(errno, "type_status");
  return file_status(type);
}

#endif

// end of OS-specific stuff

std::uintmax_t remove_all(path const& p) {
  directory_iterator end;
  std::uintmax_t count = 0;
  while (true) {
    directory_iterator first(p);
    if (first == end) break;
    if (first->is_directory()) {
      count += remove_all(first->path());
    } else {
      remove(first->path());
      ++count;
    }
  }
  remove(p);
  ++count;
  return count;
}

filesystem_error::~filesystem_error() {}

path::path(char const* source) : impl(source) {}

path::path(std::string const& source) : impl(source) {}

path::value_type const* path::c_str() const noexcept {
  return native().c_str();
}

const path::string_type& path::native() const noexcept { return impl; }

std::string path::string() const { return impl; }

path path::filename() const {
  auto const last_sep_pos = impl.find_last_of(preferred_separator);
  if (last_sep_pos == std::string::npos) return impl;
  return impl.substr(last_sep_pos + 1, std::string::npos);
}

path path::parent_path() const {
  auto const last_sep_pos = impl.find_last_of(preferred_separator);
  if (last_sep_pos == std::string::npos) return impl;
  return impl.substr(0, last_sep_pos);
}

path path::extension() const {
  auto const filename_str = filename().native();
  auto const last_dot_pos = filename_str.find_last_of('.');
  // If the pathname is either . or .., or if filename() does not contain the .
  // character, then empty path is returned.
  if (last_dot_pos == std::string::npos) return path();
  if (filename_str == "." || filename_str == "..") return path();
  // If the first character in the filename is a period, that period is ignored
  // (a filename like ".profile" is not treated as an extension)
  if (last_dot_pos == 0) return path();
  // If the filename() component of the generic-format path contains a period
  // (.), and is not one of the special filesystem elements dot or dot-dot, then
  // the extension is the substring beginning at the rightmost period (including
  // the period) and until the end of the pathname.
  return filename_str.substr(last_dot_pos, std::string::npos);
}

// see extension
path path::stem() const {
  auto const filename_str = filename().native();
  auto const last_dot_pos = filename_str.find_last_of('.');
  if (last_dot_pos == std::string::npos) return filename_str;
  if (filename_str == "." || filename_str == "..") return filename_str;
  if (last_dot_pos == 0) return filename_str;
  return filename_str.substr(0, last_dot_pos);
}

path& path::operator/=(path const& p) {
  impl.push_back(preferred_separator);
  impl += p.impl;
  return *this;
}

path& path::operator/=(std::string const& s) {
  impl.push_back(preferred_separator);
  impl += s;
  return *this;
}

path& path::operator/=(const char* s) {
  impl.push_back(preferred_separator);
  impl += s;
  return *this;
}

path& path::operator+=(path const& p) {
  impl += p.impl;
  return *this;
}

path& path::operator+=(std::string const& s) {
  impl += s;
  return *this;
}

path& path::operator+=(const char* s) {
  impl += s;
  return *this;
}

path& path::operator+=(char c) {
  impl.push_back(c);
  return *this;
}

path operator/(path const& a, path const& b) {
  path c = a;
  c /= b;
  return c;
}

std::ostream& operator<<(std::ostream& os, path const& p) {
  os << p.c_str();
  return os;
}

filesystem_error::filesystem_error(int ev, const char* what_arg)
    : std::system_error(ev, std::system_category(), what_arg) {}

file_status::file_status(file_type const& type_in) : type_variable(type_in) {}

file_type file_status::type() const noexcept { return type_variable; }

directory_entry::directory_entry(class path const& p) : path_variable(p) {}

bool directory_entry::is_regular_file() const {
  return status(path_variable).type() == file_type::regular;
}

bool directory_entry::is_directory() const {
  return status(path_variable).type() == file_type::directory;
}

bool directory_entry::is_symlink() const {
  return status(path_variable).type() == file_type::symlink;
}

directory_iterator::directory_iterator() : impl(new IteratorImpl()) {}

directory_iterator::directory_iterator(directory_iterator&& other) {
  delete impl;
  impl = other.impl;
  other.impl = new IteratorImpl();
}

directory_iterator::directory_iterator(path const& p)
    : impl(new IteratorImpl(p)) {
  if (!impl->is_end()) entry = impl->deref();
}

directory_iterator::~directory_iterator() {
  delete impl;
  impl = nullptr;
}

directory_iterator& directory_iterator::operator++() {
  impl->increment();
  if (!impl->is_end()) entry = impl->deref();
  return *this;
}

const directory_entry& directory_iterator::operator*() const { return entry; }

const directory_entry* directory_iterator::operator->() const { return &entry; }

bool directory_iterator::operator==(directory_iterator const& other) const {
  return impl->equal(*other.impl);
}

}  // namespace filesystem

}  // namespace Omega_h
