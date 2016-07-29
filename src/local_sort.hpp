#ifndef LOCAL_SORT_HPP
#define LOCAL_SORT_HPP

#include "internal.hpp"

namespace osh {

template <typename Item, Int cap>
struct LocalMergeSort {
  Int n;
  Item array[cap];
  Item scratch[cap];
  INLINE void run() {
    top_down_split_merge(0, n);
  }
  INLINE void top_down_split_merge(Int begin, Int end) {
    if (end - begin < 2) return;
    auto middle = (end + begin) / 2;
    top_down_split_merge(begin, middle);
    top_down_split_merge(middle, end);
    top_down_merge(begin, middle, end);
    copy_back(begin, end);
  }
  INLINE void top_down_merge(Int begin, Int middle, Int end) {
    auto i = begin;
    auto j = middle;
    for (auto k = begin; k < end; ++k) {
      if (i < middle && (j >= end || (!(array[j] < array[i]))))
        scratch[k] = array[i++];
      else
        scratch[k] = array[j++];
    }
  }
  INLINE void copy_back(Int begin, Int end) {
    for (auto k = begin; k < end; ++k)
      array[k] = scratch[k];
  }
};

}

#endif
