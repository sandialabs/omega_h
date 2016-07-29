#ifndef LOCAL_SORT_HPP
#define LOCAL_SORT_HPP

#include "algorithm.hpp"
#include "internal.hpp"

namespace osh {

template <typename Item>
struct SortableArray {
  typedef Item value_type;
  Item* ptr;
  SortableArray(Item array[]):ptr(array) {}
  bool is_less(Int i, Int j) const { return ptr[i] < ptr[j]; }
  value_type const& get(Int i) const { return ptr[i]; }
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
    if (i < middle && (j >= end || (!sortable.is_less(j, i))))
      scratch[k] = sortable.get(i++);
    else
      scratch[k] = sortable.get(j++);
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
    Int k = i;
    for (Int j = i + 1; j < n; ++j) {
      if (sortable.is_less(j, k)) k = j;
    }
    auto tmp = sortable.get(i);
    sortable.set(i, sortable.get(k));
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
