#include <Omega_h_profile.hpp>
#include <algorithm>
#include <iostream>
#include <queue>

namespace Omega_h {
namespace profile {

OMEGA_H_DLL History* global_singleton_history = nullptr;

History::History() : current_frame(invalid), last_root(invalid) {}

std::size_t History::first(std::size_t parent_index) const {
  if (parent_index != invalid) return frames[parent_index].first_child;
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

double History::time(std::size_t frame) const {
  return frames[frame].total_runtime;
}

std::size_t History::calls(std::size_t frame) const {
  return frames[frame].number_of_calls;
}

struct PreOrderIterator {
  using reference = std::size_t;
  PreOrderIterator& operator++() {
    frame = history.pre_order_next(frame);
    return *this;
  }
  reference operator*() { return frame; }
  PreOrderIterator(History const& history_in, std::size_t i_in)
      : history(history_in), frame(i_in) {}
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

static void simple_print(
    History const& history, std::vector<std::size_t> const& depths) {
  for (auto frame : history) {
    std::size_t depth = depths[frame];
    for (std::size_t i = 0; i < depth; ++i) std::cout << "  ";
    std::cout << history.get_name(frame) << ' '
              << history.frames[frame].total_runtime << '\n';
  }
}

void simple_print(History const& history) {
  auto depths = compute_depths(history);
  simple_print(history, depths);
}

History invert(History const& h) {
  History invh;
  std::queue<std::size_t> q;
  for (std::size_t s = h.first(invalid); s != invalid; s = h.next(s)) {
    q.push(s);
  }
  while (!q.empty()) {
    auto node = q.front();
    q.pop();
    auto self_time = h.time(node);
    auto calls = h.calls(node);
    for (auto child = h.first(node); child != invalid; child = h.next(child)) {
      self_time -= h.time(child);
      q.push(child);
    }
    self_time = std::max(self_time,
        0.);  // floating-point may give negative epsilon instead of zero
    auto inv_node = invalid;
    for (; node != invalid; node = h.parent(node)) {
      auto name = h.get_name(node);
      inv_node = invh.find_or_create_child_of(inv_node, name);
      invh.frames[inv_node].total_runtime += self_time;
      invh.frames[inv_node].number_of_calls += calls;
    }
  }
  return invh;
}

static void print_time_sorted_recursive(History const& h, std::size_t frame,
    std::vector<std::size_t> const& depths) {
  std::vector<std::size_t> child_frames;
  for (std::size_t child = h.first(frame); child != invalid;
       child = h.next(child)) {
    child_frames.push_back(child);
  }
  std::stable_sort(begin(child_frames), end(child_frames),
      [&](std::size_t a, std::size_t b) { return h.time(a) > h.time(b); });
  for (auto child : child_frames) {
    std::size_t depth = depths[child];
    for (std::size_t i = 0; i < depth; ++i) std::cout << "|  ";
    std::cout << h.get_name(child) << ' ' << h.time(child) << ' '
              << h.calls(child) << '\n';
    print_time_sorted_recursive(h, child, depths);
  }
}

void print_time_sorted(History const& h) {
  auto depths = compute_depths(h);
  print_time_sorted_recursive(h, invalid, depths);
}

void print_top_down_and_bottom_up(History const& h) {
  std::cout << "\n";
  std::cout << "TOP-DOWN:\n";
  std::cout << "=========\n";
  print_time_sorted(h);
  auto h_inv = invert(h);
  std::cout << "\n";
  std::cout << "BOTTOM-UP:\n";
  std::cout << "==========\n";
  print_time_sorted(h_inv);
}

}  // namespace profile
}  // namespace Omega_h
