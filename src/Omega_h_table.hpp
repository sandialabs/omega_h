#ifndef OMEGA_H_TABLE_HPP
#define OMEGA_H_TABLE_HPP

#include <Omega_h_std_vector.hpp>

namespace Omega_h {

/* pretty simple 2D array */
template <typename T>
struct Table {
  std::vector<T> data;
  int ncols;
  using Ref = typename std::vector<T>::reference;
  using ConstRef = typename std::vector<T>::const_reference;
  Table() = default;
  Table(int ncols_init, int nrows_reserve) : ncols(ncols_init) {
    OMEGA_H_CHECK(0 <= ncols_init);
    reserve(data, ncols * nrows_reserve);
  }
};

template <typename T>
int get_nrows(Table<T> const& t) {
  OMEGA_H_CHECK(t.ncols > 0);
  OMEGA_H_CHECK(size(t.data) % t.ncols == 0);
  return size(t.data) / t.ncols;
}

template <typename T>
int get_ncols(Table<T> const& t) {
  return t.ncols;
}

template <typename T>
void resize(Table<T>& t, int new_nrows, int new_ncols) {
  OMEGA_H_CHECK(new_ncols == t.ncols);  // pretty specialized right now
  Omega_h::resize(t.data, new_nrows * t.ncols);
}

template <typename T>
typename Table<T>::Ref at(Table<T>& t, int row, int col) {
  OMEGA_H_CHECK(0 <= col);
  OMEGA_H_CHECK(col < t.ncols);
  OMEGA_H_CHECK(0 <= row);
  OMEGA_H_CHECK(row < get_nrows(t));
  return Omega_h::at(t.data, row * t.ncols + col);
}

template <typename T>
typename Table<T>::ConstRef at(Table<T> const& t, int row, int col) {
  OMEGA_H_CHECK(0 <= col);
  OMEGA_H_CHECK(col < t.ncols);
  OMEGA_H_CHECK(0 <= row);
  OMEGA_H_CHECK(row < get_nrows(t));
  return Omega_h::at(t.data, row * t.ncols + col);
}

}  // namespace Omega_h

#endif
