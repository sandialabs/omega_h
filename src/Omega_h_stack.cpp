#include <Omega_h_stack.hpp>
#include <iostream>

namespace Omega_h {
namespace perf {

History* global_singleton_history = nullptr;

History::History():current_frame(invalid),last_root(invalid) {}

std::size_t History::first(std::size_t parent) const {
  if (parent != invalid) return frames[parent].first_child;
  if (!frames.empty()) return 0;
  return invalid;
}

std::size_t History::next(std::size_t sibling) const {
  return frames[sibling].next_sibling;
}

std::size_t History::parent(std::size_t child) const {
  return frames[child].parent;
}

std::size_t History::pre_order_next(std::size_t frame) const {
  auto first2 = first(frame);
  if (first2 != invalid) {
    return first2;
  }
  auto next2 = next(frame);
  if (next2 != invalid) {
    return next2;
  }
  for (; frame != invalid; frame = parent(frame)) {
    auto next3 = next(frame);
    if (next3 != invalid) {
      return next3;
    }
  }
  return invalid;
}

double History::time(std::size_t frame) {
  return frames[frame].total_runtime;
}

double History::calls(std::size_t frame) {
  return frames[frame].number_of_calls;
}

struct PreOrderIterator {
  using reference = std::size_t;
  PreOrderIterator& operator++() {
    frame = history.pre_order_next(frame);
    return *this;
  }
  reference operator*() {
    return frame;
  }
  PreOrderIterator(History const& history_in, std::size_t i_in)
    :history(history_in)
    ,frame(i_in)
  {
  }
  bool operator!=(PreOrderIterator const& other) {
    return frame != other.frame;
  }
  History const& history;
  std::size_t frame;
};

static PreOrderIterator begin(History const& history) {
  return PreOrderIterator(history, history.frames.empty() ? invalid : 0);
}

static PreOrderIterator end(History const& history) {
  return PreOrderIterator(history, invalid);
}

static std::vector<std::size_t> compute_depths(History const& history) {
  std::vector<std::size_t> depths(history.frames.size());
  for (auto frame : history) {
    auto parent = history.parent(frame);
    if (parent != invalid) {
      depths[frame] = depths[parent] + 1;
    } else {
      depths[frame] = 0;
    }
  }
  return depths;
}

static void simple_print(History const& history, std::vector<std::size_t> const& depths) {
  for (auto frame : history) {
    std::size_t depth = depths[frame];
    for (std::size_t i = 0; i < depth; ++i) std::cout << "  ";
    std::cout << history.get_name(frame) << ' ' << history.frames[frame].total_runtime << '\n';
  }
}

void simple_print(History const& history) {
  auto depths = compute_depths(history);
  simple_print(history, depths);
}

}}
