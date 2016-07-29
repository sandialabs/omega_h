#ifndef LOCAL_SORT_HPP
#define LOCAL_SORT_HPP

#include "algorithm.hpp"
#include "internal.hpp"

namespace osh {

template <typename Item>
struct SortableArray {
  typedef Item value_type;
  typedef Item key_type;
  Item* ptr;
  SortableArray(Item array[]):ptr(array) {}
  key_type key(Int i) const { return ptr[i]; }
  value_type const& value(Int i) const { return ptr[i]; }
  void set(Int i, value_type const& x) const { ptr[i] = x; }
};

template <typename Sortable>
INLINE void copy_back(Sortable const& sortable, Int begin, Int end,
    typename Sortable::value_type scratch[]) {
  for (auto k = begin; k < end; ++k)
    sortable.set(k, scratch[k]);
}

template <typename Sortable>
INLINE void top_down_merge(Sortable const& sortable, Int begin, Int middle, Int end,
    typename Sortable::value_type scratch[]) {
  auto i = begin;
  auto j = middle;
  for (auto k = begin; k < end; ++k) {
    if (i < middle && (j >= end || (sortable.key(i) <= sortable.key(j))))
      scratch[k] = sortable.value(i++);
    else
      scratch[k] = sortable.value(j++);
  }
}

template <typename Sortable>
INLINE void top_down_split_merge(Sortable const& sortable, Int begin, Int end,
    typename Sortable::value_type scratch[]) {
  if (end - begin < 2) return;
  auto middle = (end + begin) / 2;
  top_down_split_merge(sortable, begin, middle, scratch);
  top_down_split_merge(sortable, middle, end, scratch);
  top_down_merge(sortable, begin, middle, end, scratch);
  copy_back(sortable, begin, end, scratch);
}

template <Int cap, typename Sortable>
INLINE void top_down_merge_sort(Sortable const& sortable, Int n) {
  typename Sortable::value_type scratch[cap];
  top_down_split_merge(sortable, 0, n, scratch);
}

template <Int cap, typename Item>
INLINE void top_down_merge_sort(Item array[], Int n) {
  top_down_merge_sort<cap>(SortableArray<Item>(array), n);
}

template <typename Sortable>
INLINE void selection_sort(Sortable const& sortable, Int n) {
  for (Int i = 0; i < n; ++i) {
    auto k = i;
    auto i_key = sortable.key(i);
    for (Int j = i + 1; j < n; ++j) {
      if (sortable.key(j) < i_key) k = j;
    }
    auto tmp = sortable.value(i);
    sortable.set(i, sortable.value(k));
    sortable.set(k, tmp);
  }
}

template <typename Item>
INLINE void selection_sort(Item array[], Int n) {
  selection_sort(SortableArray<Item>(array), n);
}

INLINE Int heapsort_parent(Int i) { return (i - 1) / 2; }
INLINE Int heapsort_left_child(Int i) { return 2 * i + 1; }
INLINE Int heapsort_right_child(Int i) { return 2 * i + 2; }

}

#endif
