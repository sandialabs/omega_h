/* Dan Ibanez took the code from HP/SGI and modified it to:
 * 1. be user-space instead of standard-namespace
 * 2. Not have Unicode foreign language comments
 * 3. Not have allocators (just uses new and delete)
 */

/*
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 */

#ifndef OMEGA_H_RB_TREE_HPP
#define OMEGA_H_RB_TREE_HPP

/*

Red-black tree class, designed for use in implementing STL
associative containers (set, multiset, map, and multimap). The
insertion and deletion algorithms are based on those in Cormen,
Leiserson, and Rivest, Introduction to Algorithms (MIT Press, 1990),
except that

(1) the header cell is maintained with links not only to the root
but also to the leftmost node of the tree, to enable constant time
begin(), and to the rightmost node of the tree, to enable linear time
performance when used with the generic set algorithms (set_union,
etc.);

(2) when a node being deleted has two children its successor node is
relinked into its place, rather than copied, so that the only
iterators invalidated are those referring to the deleted node.

*/

#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <utility>

namespace Omega_h {

namespace nonstd {
template <typename T>
void swap(T& a, T& b) {
  T c = a;
  a = b;
  b = c;
}
}  // namespace nonstd

typedef bool Rb_tree_Color_type;
static constexpr Rb_tree_Color_type S_rb_tree_red = false;
static constexpr Rb_tree_Color_type S_rb_tree_black = true;

struct Rb_tree_node_base {
  typedef Rb_tree_Color_type Color_type;
  typedef Rb_tree_node_base* Base_ptr;

  Color_type M_color;
  Base_ptr M_parent;
  Base_ptr M_left;
  Base_ptr M_right;

  static Base_ptr S_minimum(Base_ptr x) {
    while (x->M_left != nullptr) x = x->M_left;
    return x;
  }

  static Base_ptr S_maximum(Base_ptr x) {
    while (x->M_right != nullptr) x = x->M_right;
    return x;
  }
};

template <class Value>
struct Rb_tree_node : public Rb_tree_node_base {
  typedef Rb_tree_node<Value>* Link_type;
  Value M_value_field;
};

struct Rb_tree_base_iterator {
  typedef Rb_tree_node_base::Base_ptr Base_ptr;
  typedef std::ptrdiff_t difference_type;
  using iterator_category = std::bidirectional_iterator_tag;
  Base_ptr M_node;

  void M_increment() {
    if (M_node->M_right != nullptr) {
      M_node = M_node->M_right;
      while (M_node->M_left != nullptr) M_node = M_node->M_left;
    } else {
      Base_ptr y = M_node->M_parent;
      while (M_node == y->M_right) {
        M_node = y;
        y = y->M_parent;
      }
      if (M_node->M_right != y) M_node = y;
    }
  }

  void M_decrement() {
    if (M_node->M_color == S_rb_tree_red &&
        M_node->M_parent->M_parent == M_node)
      M_node = M_node->M_right;
    else if (M_node->M_left != nullptr) {
      Base_ptr y = M_node->M_left;
      while (y->M_right != nullptr) y = y->M_right;
      M_node = y;
    } else {
      Base_ptr y = M_node->M_parent;
      while (M_node == y->M_left) {
        M_node = y;
        y = y->M_parent;
      }
      M_node = y;
    }
  }
};

template <class Value, class Ref, class Ptr>
struct Rb_tree_iterator : public Rb_tree_base_iterator {
  typedef Value value_type;
  typedef Ref reference;
  typedef Ptr pointer;
  typedef Rb_tree_iterator<Value, Value&, Value*> iterator;
  typedef Rb_tree_iterator<Value, const Value&, const Value*> const_iterator;
  typedef Rb_tree_iterator<Value, Ref, Ptr> Self;
  typedef Rb_tree_node<Value>* Link_type;

  Rb_tree_iterator() {}
  Rb_tree_iterator(Link_type x) { M_node = x; }
  Rb_tree_iterator(const Self&) = default;
  Rb_tree_iterator& operator=(Rb_tree_iterator const&) = default;

  reference operator*() const { return Link_type(M_node)->M_value_field; }
  pointer operator->() const { return &(operator*()); }

  Self& operator++() {
    M_increment();
    return *this;
  }
  Self operator++(int) {
    Self tmp = *this;
    M_increment();
    return tmp;
  }

  Self& operator--() {
    M_decrement();
    return *this;
  }
  Self operator--(int) {
    Self tmp = *this;
    M_decrement();
    return tmp;
  }
};

inline bool operator==(
    const Rb_tree_base_iterator& x, const Rb_tree_base_iterator& y) {
  return x.M_node == y.M_node;
}

inline bool operator!=(
    const Rb_tree_base_iterator& x, const Rb_tree_base_iterator& y) {
  return x.M_node != y.M_node;
}

inline void Rb_tree_rotate_left(
    Rb_tree_node_base* x, Rb_tree_node_base*& root) {
  Rb_tree_node_base* y = x->M_right;
  x->M_right = y->M_left;
  if (y->M_left != nullptr) y->M_left->M_parent = x;
  y->M_parent = x->M_parent;

  if (x == root)
    root = y;
  else if (x == x->M_parent->M_left)
    x->M_parent->M_left = y;
  else
    x->M_parent->M_right = y;
  y->M_left = x;
  x->M_parent = y;
}

inline void Rb_tree_rotate_right(
    Rb_tree_node_base* x, Rb_tree_node_base*& root) {
  Rb_tree_node_base* y = x->M_left;
  x->M_left = y->M_right;
  if (y->M_right != nullptr) y->M_right->M_parent = x;
  y->M_parent = x->M_parent;

  if (x == root)
    root = y;
  else if (x == x->M_parent->M_right)
    x->M_parent->M_right = y;
  else
    x->M_parent->M_left = y;
  y->M_right = x;
  x->M_parent = y;
}

inline void Rb_tree_rebalance(Rb_tree_node_base* x, Rb_tree_node_base*& root) {
  x->M_color = S_rb_tree_red;
  while (x != root && x->M_parent->M_color == S_rb_tree_red) {
    if (x->M_parent == x->M_parent->M_parent->M_left) {
      Rb_tree_node_base* y = x->M_parent->M_parent->M_right;
      if (y && y->M_color == S_rb_tree_red) {
        x->M_parent->M_color = S_rb_tree_black;
        y->M_color = S_rb_tree_black;
        x->M_parent->M_parent->M_color = S_rb_tree_red;
        x = x->M_parent->M_parent;
      } else {
        if (x == x->M_parent->M_right) {
          x = x->M_parent;
          Rb_tree_rotate_left(x, root);
        }
        x->M_parent->M_color = S_rb_tree_black;
        x->M_parent->M_parent->M_color = S_rb_tree_red;
        Rb_tree_rotate_right(x->M_parent->M_parent, root);
      }
    } else {
      Rb_tree_node_base* y = x->M_parent->M_parent->M_left;
      if (y && y->M_color == S_rb_tree_red) {
        x->M_parent->M_color = S_rb_tree_black;
        y->M_color = S_rb_tree_black;
        x->M_parent->M_parent->M_color = S_rb_tree_red;
        x = x->M_parent->M_parent;
      } else {
        if (x == x->M_parent->M_left) {
          x = x->M_parent;
          Rb_tree_rotate_right(x, root);
        }
        x->M_parent->M_color = S_rb_tree_black;
        x->M_parent->M_parent->M_color = S_rb_tree_red;
        Rb_tree_rotate_left(x->M_parent->M_parent, root);
      }
    }
  }
  root->M_color = S_rb_tree_black;
}

inline Rb_tree_node_base* Rb_tree_rebalance_for_erase(Rb_tree_node_base* z,
    Rb_tree_node_base*& root, Rb_tree_node_base*& leftmost,
    Rb_tree_node_base*& rightmost) {
  Rb_tree_node_base* y = z;
  Rb_tree_node_base* x = nullptr;
  Rb_tree_node_base* x_parent = nullptr;
  if (y->M_left == nullptr)        // z has at most one non-null child. y == z.
    x = y->M_right;                // x might be null.
  else if (y->M_right == nullptr)  // z has exactly one non-null child. y == z.
    x = y->M_left;                 // x is not null.
  else {                           // z has two non-null children.  Set y to
    y = y->M_right;                //   z's successor.  x might be null.
    while (y->M_left != nullptr) y = y->M_left;
    x = y->M_right;
  }
  if (y != z) {  // relink y in place of z.  y is z's successor
    z->M_left->M_parent = y;
    y->M_left = z->M_left;
    if (y != z->M_right) {
      x_parent = y->M_parent;
      if (x) x->M_parent = y->M_parent;
      y->M_parent->M_left = x;  // y must be a child of M_left
      y->M_right = z->M_right;
      z->M_right->M_parent = y;
    } else
      x_parent = y;
    if (root == z)
      root = y;
    else if (z->M_parent->M_left == z)
      z->M_parent->M_left = y;
    else
      z->M_parent->M_right = y;
    y->M_parent = z->M_parent;
    nonstd::swap(y->M_color, z->M_color);
    y = z;
    // y now points to node to be actually deleted
  } else {  // y == z
    x_parent = y->M_parent;
    if (x) x->M_parent = y->M_parent;
    if (root == z)
      root = x;
    else if (z->M_parent->M_left == z)
      z->M_parent->M_left = x;
    else
      z->M_parent->M_right = x;
    if (leftmost == z) {
      if (z->M_right == nullptr) {  // z->M_left must be null also
        leftmost = z->M_parent;
        // makes leftmost == M_header if z == root
      } else {
        leftmost = Rb_tree_node_base::S_minimum(x);
      }
    }
    if (rightmost == z) {
      if (z->M_left == nullptr) {  // z->M_right must be null also
        rightmost = z->M_parent;
        // makes rightmost == M_header if z == root
      } else {  // x == z->M_left
        rightmost = Rb_tree_node_base::S_maximum(x);
      }
    }
  }
  if (y->M_color != S_rb_tree_red) {
    while (x != root && (x == nullptr || x->M_color == S_rb_tree_black))
      if (x == x_parent->M_left) {
        Rb_tree_node_base* w = x_parent->M_right;
        if (w->M_color == S_rb_tree_red) {
          w->M_color = S_rb_tree_black;
          x_parent->M_color = S_rb_tree_red;
          Rb_tree_rotate_left(x_parent, root);
          w = x_parent->M_right;
        }
        if ((w->M_left == nullptr || w->M_left->M_color == S_rb_tree_black) &&
            (w->M_right == nullptr || w->M_right->M_color == S_rb_tree_black)) {
          w->M_color = S_rb_tree_red;
          x = x_parent;
          x_parent = x_parent->M_parent;
        } else {
          if (w->M_right == nullptr || w->M_right->M_color == S_rb_tree_black) {
            if (w->M_left) w->M_left->M_color = S_rb_tree_black;
            w->M_color = S_rb_tree_red;
            Rb_tree_rotate_right(w, root);
            w = x_parent->M_right;
          }
          w->M_color = x_parent->M_color;
          x_parent->M_color = S_rb_tree_black;
          if (w->M_right) w->M_right->M_color = S_rb_tree_black;
          Rb_tree_rotate_left(x_parent, root);
          break;
        }
      } else {  // same as above, with M_right <-> M_left.
        Rb_tree_node_base* w = x_parent->M_left;
        if (w->M_color == S_rb_tree_red) {
          w->M_color = S_rb_tree_black;
          x_parent->M_color = S_rb_tree_red;
          Rb_tree_rotate_right(x_parent, root);
          w = x_parent->M_left;
        }
        if ((w->M_right == nullptr || w->M_right->M_color == S_rb_tree_black) &&
            (w->M_left == nullptr || w->M_left->M_color == S_rb_tree_black)) {
          w->M_color = S_rb_tree_red;
          x = x_parent;
          x_parent = x_parent->M_parent;
        } else {
          if (w->M_left == nullptr || w->M_left->M_color == S_rb_tree_black) {
            if (w->M_right) w->M_right->M_color = S_rb_tree_black;
            w->M_color = S_rb_tree_red;
            Rb_tree_rotate_left(w, root);
            w = x_parent->M_left;
          }
          w->M_color = x_parent->M_color;
          x_parent->M_color = S_rb_tree_black;
          if (w->M_left) w->M_left->M_color = S_rb_tree_black;
          Rb_tree_rotate_right(x_parent, root);
          break;
        }
      }
    if (x) x->M_color = S_rb_tree_black;
  }
  return y;
}

// Base class to encapsulate the differences between old SGI-style
// allocators and standard-conforming allocators.  In order to avoid
// having an empty base class, we arbitrarily move one of rb_tree's
// data members into the base class.

template <class Tp>
struct Rb_tree_base {
  Rb_tree_base() : M_header(nullptr) { M_header = M_get_node(); }
  ~Rb_tree_base() { M_put_node(M_header); }

 protected:
  Rb_tree_node<Tp>* M_header;

  Rb_tree_node<Tp>* M_get_node() {
    return static_cast<Rb_tree_node<Tp>*>(
        std::malloc(sizeof(Rb_tree_node<Tp>)));
  }
  void M_put_node(Rb_tree_node<Tp>* p) { std::free(p); }
};

template <class Key, class Value, class KeyOfValue, class Compare>
class Rb_tree : protected Rb_tree_base<Value> {
  typedef Rb_tree_base<Value> Base;

 protected:
  typedef Rb_tree_node_base* Base_ptr;
  typedef Rb_tree_node<Value> Node_type;
  typedef Rb_tree_Color_type Color_type;

 public:
  typedef Key key_type;
  typedef Value value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef Node_type* Link_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

 protected:
  using Base::M_get_node;
  using Base::M_header;
  using Base::M_put_node;

 protected:
  Link_type M_create_node(const value_type& x) {
    Link_type tmp = M_get_node();
    try {
      new (&tmp->M_value_field) Value(x);
    } catch (...) {
      M_put_node(tmp);
      throw;
    }
    return tmp;
  }

  Link_type M_create_node(value_type&& x) {
    Link_type tmp = M_get_node();
    try {
      new (&tmp->M_value_field) Value(x);
    } catch (...) {
      M_put_node(tmp);
      throw;
    }
    return tmp;
  }

  Link_type M_clone_node(Link_type x) {
    Link_type tmp = M_create_node(x->M_value_field);
    tmp->M_color = x->M_color;
    tmp->M_left = nullptr;
    tmp->M_right = nullptr;
    return tmp;
  }

  void destroy_node(Link_type p) {
    p->M_value_field.~Value();
    M_put_node(p);
  }

 protected:
  size_type M_node_count;  // keeps track of size of tree
  Compare M_key_compare;

// although the following casts may technically break aliasing rules,
// they provide an excellent separation of the generic tree code (base)
// from the specific value type selected by the user
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#endif
  Link_type& M_root() const { return (Link_type&)(M_header->M_parent); }
  Link_type& M_leftmost() const { return (Link_type&)(M_header->M_left); }
  Link_type& M_rightmost() const { return (Link_type&)(M_header->M_right); }

  static Link_type& S_left(Link_type x) { return (Link_type&)(x->M_left); }
  static Link_type& S_right(Link_type x) { return (Link_type&)(x->M_right); }
  static Link_type& S_parent(Link_type x) { return (Link_type&)(x->M_parent); }
  static reference S_value(Link_type x) { return x->M_value_field; }
  static const Key& S_key(Link_type x) { return KeyOfValue()(S_value(x)); }
  static Color_type& S_color(Link_type x) { return (x->M_color); }

  static Link_type& S_left(Base_ptr x) { return (Link_type&)(x->M_left); }
  static Link_type& S_right(Base_ptr x) { return (Link_type&)(x->M_right); }
  static Link_type& S_parent(Base_ptr x) { return (Link_type&)(x->M_parent); }
  static reference S_value(Base_ptr x) { return Link_type(x)->M_value_field; }
  static const Key& S_key(Base_ptr x) {
    return KeyOfValue()(S_value(Link_type(x)));
  }
  static Color_type& S_color(Base_ptr x) { return Link_type(x)->M_color; }
#ifdef __clang__
#pragma clang diagnostic pop
#endif

  static Link_type S_minimum(Link_type x) {
    return static_cast<Link_type>(Rb_tree_node_base::S_minimum(x));
  }

  static Link_type S_maximum(Link_type x) {
    return static_cast<Link_type>(Rb_tree_node_base::S_maximum(x));
  }

 public:
  typedef Rb_tree_iterator<value_type, reference, pointer> iterator;
  typedef Rb_tree_iterator<value_type, const_reference, const_pointer>
      const_iterator;

 private:
  iterator M_insert(Base_ptr x, Base_ptr y, const value_type& v);
  iterator M_insert(Base_ptr x, Base_ptr y, value_type&& v);
  Link_type M_copy(Link_type x, Link_type p);
  void M_erase(Link_type x);

 public:
  // allocation/deallocation
  Rb_tree() : Base(), M_node_count(0), M_key_compare() { M_empty_initialize(); }

  Rb_tree(const Compare& comp) : Base(), M_node_count(0), M_key_compare(comp) {
    M_empty_initialize();
  }

  Rb_tree(const Rb_tree<Key, Value, KeyOfValue, Compare>& x)
      : Base(), M_node_count(0), M_key_compare(x.M_key_compare) {
    if (x.M_root() == nullptr)
      M_empty_initialize();
    else {
      S_color(M_header) = S_rb_tree_red;
      M_root() = M_copy(x.M_root(), M_header);
      M_leftmost() = S_minimum(M_root());
      M_rightmost() = S_maximum(M_root());
    }
    M_node_count = x.M_node_count;
  }
  ~Rb_tree() { clear(); }
  Rb_tree<Key, Value, KeyOfValue, Compare>& operator=(
      const Rb_tree<Key, Value, KeyOfValue, Compare>& x);

 private:
  void M_empty_initialize() {
    S_color(M_header) = S_rb_tree_red;  // used to distinguish header from
                                        // root, in iterator.operator++
    M_root() = nullptr;
    M_leftmost() = M_header;
    M_rightmost() = M_header;
  }

 public:
  // accessors:
  Compare key_comp() const { return M_key_compare; }
  iterator begin() { return M_leftmost(); }
  const_iterator begin() const { return M_leftmost(); }
  iterator end() { return M_header; }
  const_iterator end() const { return M_header; }
  bool empty() const { return M_node_count == 0; }
  size_type size() const { return M_node_count; }
  size_type max_size() const { return size_type(-1); }

  void swap(Rb_tree<Key, Value, KeyOfValue, Compare>& t) {
    nonstd::swap(M_header, t.M_header);
    nonstd::swap(M_node_count, t.M_node_count);
    nonstd::swap(M_key_compare, t.M_key_compare);
  }

 public:
  // insert/erase
  std::pair<iterator, bool> insert(const value_type& x);
  std::pair<iterator, bool> insert(value_type&& x);

  iterator insert(iterator position, const value_type& x);
  iterator insert(iterator position, value_type&& x);

  template <class InputIterator>
  void insert(InputIterator first, InputIterator last);

  void erase(iterator position);
  size_type erase(const key_type& x);
  void erase(iterator first, iterator last);
  void erase(const key_type* first, const key_type* last);
  void clear() {
    if (M_node_count != 0) {
      M_erase(M_root());
      M_leftmost() = M_header;
      M_root() = nullptr;
      M_rightmost() = M_header;
      M_node_count = 0;
    }
  }

 public:
  // set operations:
  iterator find(const key_type& x);
  const_iterator find(const key_type& x) const;
  size_type count(const key_type& x) const;
  iterator lower_bound(const key_type& x);
  const_iterator lower_bound(const key_type& x) const;
  iterator upper_bound(const key_type& x);
  const_iterator upper_bound(const key_type& x) const;
  std::pair<iterator, iterator> equal_range(const key_type& x);
  std::pair<const_iterator, const_iterator> equal_range(
      const key_type& x) const;
};

template <class Key, class Value, class KeyOfValue, class Compare>
inline bool operator==(const Rb_tree<Key, Value, KeyOfValue, Compare>& x,
    const Rb_tree<Key, Value, KeyOfValue, Compare>& y) {
  return x.size() == y.size() && equal(x.begin(), x.end(), y.begin());
}

template <class Key, class Value, class KeyOfValue, class Compare>
inline bool operator<(const Rb_tree<Key, Value, KeyOfValue, Compare>& x,
    const Rb_tree<Key, Value, KeyOfValue, Compare>& y) {
  return lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

template <class Key, class Value, class KeyOfValue, class Compare>
Rb_tree<Key, Value, KeyOfValue, Compare>&
Rb_tree<Key, Value, KeyOfValue, Compare>::operator=(
    const Rb_tree<Key, Value, KeyOfValue, Compare>& x) {
  if (this != &x) {  // Note that Key may be a constant type.
    clear();
    M_node_count = 0;
    M_key_compare = x.M_key_compare;
    if (x.M_root() == nullptr) {
      M_root() = nullptr;
      M_leftmost() = M_header;
      M_rightmost() = M_header;
    } else {
      M_root() = M_copy(x.M_root(), M_header);
      M_leftmost() = S_minimum(M_root());
      M_rightmost() = S_maximum(M_root());
      M_node_count = x.M_node_count;
    }
  }
  return *this;
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::M_insert(
    Base_ptr x_, Base_ptr y_, const Value& v) {
  Link_type x = Link_type(x_);
  Link_type y = Link_type(y_);
  Link_type z;

  if (y == M_header || x != nullptr ||
      M_key_compare(KeyOfValue()(v), S_key(y))) {
    z = M_create_node(v);
    S_left(y) = z;  // also makes M_leftmost() = z
                    //    when y == M_header
    if (y == M_header) {
      M_root() = z;
      M_rightmost() = z;
    } else if (y == M_leftmost())
      M_leftmost() = z;  // maintain M_leftmost() pointing to min node
  } else {
    z = M_create_node(v);
    S_right(y) = z;
    if (y == M_rightmost())
      M_rightmost() = z;  // maintain M_rightmost() pointing to max node
  }
  S_parent(z) = y;
  S_left(z) = nullptr;
  S_right(z) = nullptr;
  Rb_tree_rebalance(z, M_header->M_parent);
  ++M_node_count;
  return iterator(z);
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::M_insert(
    Base_ptr x_, Base_ptr y_, Value&& v) {
  Link_type x = Link_type(x_);
  Link_type y = Link_type(y_);
  Link_type z;

  if (y == M_header || x != nullptr ||
      M_key_compare(KeyOfValue()(v), S_key(y))) {
    z = M_create_node(std::move(v));
    S_left(y) = z;  // also makes M_leftmost() = z
                    //    when y == M_header
    if (y == M_header) {
      M_root() = z;
      M_rightmost() = z;
    } else if (y == M_leftmost())
      M_leftmost() = z;  // maintain M_leftmost() pointing to min node
  } else {
    z = M_create_node(std::move(v));
    S_right(y) = z;
    if (y == M_rightmost())
      M_rightmost() = z;  // maintain M_rightmost() pointing to max node
  }
  S_parent(z) = y;
  S_left(z) = nullptr;
  S_right(z) = nullptr;
  Rb_tree_rebalance(z, M_header->M_parent);
  ++M_node_count;
  return iterator(z);
}

template <class Key, class Value, class KeyOfValue, class Compare>
std::pair<typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator, bool>
Rb_tree<Key, Value, KeyOfValue, Compare>::insert(const Value& v) {
  Link_type y = M_header;
  Link_type x = M_root();
  bool comp = true;
  while (x != nullptr) {
    y = x;
    comp = M_key_compare(KeyOfValue()(v), S_key(x));
    x = comp ? S_left(x) : S_right(x);
  }
  iterator j = iterator(y);
  if (comp) {
    if (j == begin()) {
      return std::pair<iterator, bool>(M_insert(x, y, v), true);
    } else {
      --j;
    }
  }
  if (M_key_compare(S_key(j.M_node), KeyOfValue()(v))) {
    return std::pair<iterator, bool>(M_insert(x, y, v), true);
  }
  return std::pair<iterator, bool>(j, false);
}

template <class Key, class Value, class KeyOfValue, class Compare>
std::pair<typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator, bool>
Rb_tree<Key, Value, KeyOfValue, Compare>::insert(Value&& v) {
  Link_type y = M_header;
  Link_type x = M_root();
  bool comp = true;
  while (x != nullptr) {
    y = x;
    comp = M_key_compare(KeyOfValue()(v), S_key(x));
    x = comp ? S_left(x) : S_right(x);
  }
  iterator j = iterator(y);
  if (comp) {
    if (j == begin()) {
      return std::pair<iterator, bool>(M_insert(x, y, std::move(v)), true);
    } else {
      --j;
    }
  }
  if (M_key_compare(S_key(j.M_node), KeyOfValue()(v))) {
    return std::pair<iterator, bool>(M_insert(x, y, std::move(v)), true);
  }
  return std::pair<iterator, bool>(j, false);
}

template <class Key, class Val, class KeyOfValue, class Compare>
typename Rb_tree<Key, Val, KeyOfValue, Compare>::iterator
Rb_tree<Key, Val, KeyOfValue, Compare>::insert(
    iterator position, const Val& v) {
  if (position.M_node == M_header->M_left) {  // begin()
    if (size() > 0 && M_key_compare(KeyOfValue()(v), S_key(position.M_node)))
      return M_insert(position.M_node, position.M_node, v);
    // first argument just needs to be non-null
    else
      return insert(v).first;
  } else if (position.M_node == M_header) {  // end()
    if (M_key_compare(S_key(M_rightmost()), KeyOfValue()(v)))
      return M_insert(nullptr, M_rightmost(), v);
    else
      return insert(v).first;
  } else {
    iterator before = position;
    --before;
    if (M_key_compare(S_key(before.M_node), KeyOfValue()(v)) &&
        M_key_compare(KeyOfValue()(v), S_key(position.M_node))) {
      if (S_right(before.M_node) == nullptr)
        return M_insert(nullptr, before.M_node, v);
      else
        return M_insert(position.M_node, position.M_node, v);
      // first argument just needs to be non-null
    } else
      return insert(v).first;
  }
}

template <class Key, class Val, class KeyOfValue, class Compare>
typename Rb_tree<Key, Val, KeyOfValue, Compare>::iterator
Rb_tree<Key, Val, KeyOfValue, Compare>::insert(iterator position, Val&& v) {
  if (position.M_node == M_header->M_left) {  // begin()
    if (size() > 0 && M_key_compare(KeyOfValue()(v), S_key(position.M_node)))
      return M_insert(position.M_node, position.M_node, std::move(v));
    // first argument just needs to be non-null
    else
      return insert(std::move(v)).first;
  } else if (position.M_node == M_header) {  // end()
    if (M_key_compare(S_key(M_rightmost()), KeyOfValue()(v)))
      return M_insert(nullptr, M_rightmost(), std::move(v));
    else
      return insert(std::move(v)).first;
  } else {
    iterator before = position;
    --before;
    if (M_key_compare(S_key(before.M_node), KeyOfValue()(v)) &&
        M_key_compare(KeyOfValue()(v), S_key(position.M_node))) {
      if (S_right(before.M_node) == nullptr)
        return M_insert(nullptr, before.M_node, std::move(v));
      else
        return M_insert(position.M_node, position.M_node, std::move(v));
      // first argument just needs to be non-null
    } else
      return insert(std::move(v)).first;
  }
}

template <class Key, class Val, class KeyOfValue, class Cmp>
template <class II>
void Rb_tree<Key, Val, KeyOfValue, Cmp>::insert(II first, II last) {
  for (; first != last; ++first) insert(*first);
}

template <class Key, class Value, class KeyOfValue, class Compare>
inline void Rb_tree<Key, Value, KeyOfValue, Compare>::erase(iterator position) {
  Link_type y = Rb_tree_rebalance_for_erase(
      position.M_node, M_header->M_parent, M_header->M_left, M_header->M_right);
  destroy_node(y);
  --M_node_count;
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::size_type
Rb_tree<Key, Value, KeyOfValue, Compare>::erase(const Key& x) {
  std::pair<iterator, iterator> p = equal_range(x);
  size_type n = size_type(std::distance(p.first, p.second));
  erase(p.first, p.second);
  return n;
}

template <class Key, class Val, class KeyOfValue, class Compare>
typename Rb_tree<Key, Val, KeyOfValue, Compare>::Link_type
Rb_tree<Key, Val, KeyOfValue, Compare>::M_copy(Link_type x, Link_type p) {
  // structural copy.  x and p must be non-null.
  Link_type top = M_clone_node(x);
  top->M_parent = p;

  try {
    if (x->M_right) top->M_right = M_copy(S_right(x), top);
    p = top;
    x = S_left(x);

    while (x != nullptr) {
      Link_type y = M_clone_node(x);
      p->M_left = y;
      y->M_parent = p;
      if (x->M_right) y->M_right = M_copy(S_right(x), y);
      p = y;
      x = S_left(x);
    }
  } catch (...) {
    M_erase(top);
    throw;
  }

  return top;
}

template <class Key, class Value, class KeyOfValue, class Compare>
void Rb_tree<Key, Value, KeyOfValue, Compare>::M_erase(
    Link_type x) {  // erase without rebalancing
  while (x != nullptr) {
    M_erase(S_right(x));
    Link_type y = S_left(x);
    destroy_node(x);
    x = y;
  }
}

template <class Key, class Value, class KeyOfValue, class Compare>
void Rb_tree<Key, Value, KeyOfValue, Compare>::erase(
    iterator first, iterator last) {
  if (first == begin() && last == end())
    clear();
  else
    while (first != last) erase(first++);
}

template <class Key, class Value, class KeyOfValue, class Compare>
void Rb_tree<Key, Value, KeyOfValue, Compare>::erase(
    const Key* first, const Key* last) {
  while (first != last) erase(*first++);
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::find(const Key& k) {
  Link_type y = M_header;  // Last node which is not less than k.
  Link_type x = M_root();  // Current node.

  while (x != nullptr)
    if (!M_key_compare(S_key(x), k))
      y = x, x = S_left(x);
    else
      x = S_right(x);

  iterator j = iterator(y);
  return (j == end() || M_key_compare(k, S_key(j.M_node))) ? end() : j;
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::const_iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::find(const Key& k) const {
  Link_type y = M_header; /* Last node which is not less than k. */
  Link_type x = M_root(); /* Current node. */

  while (x != nullptr) {
    if (!M_key_compare(S_key(x), k))
      y = x, x = S_left(x);
    else
      x = S_right(x);
  }
  const_iterator j = const_iterator(y);
  return (j == end() || M_key_compare(k, S_key(j.M_node))) ? end() : j;
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::size_type
Rb_tree<Key, Value, KeyOfValue, Compare>::count(const Key& k) const {
  std::pair<const_iterator, const_iterator> p = equal_range(k);
  return size_type(std::distance(p.first, p.second));
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::lower_bound(const Key& k) {
  Link_type y = M_header; /* Last node which is not less than k. */
  Link_type x = M_root(); /* Current node. */

  while (x != nullptr)
    if (!M_key_compare(S_key(x), k))
      y = x, x = S_left(x);
    else
      x = S_right(x);

  return iterator(y);
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::const_iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::lower_bound(const Key& k) const {
  Link_type y = M_header; /* Last node which is not less than k. */
  Link_type x = M_root(); /* Current node. */

  while (x != nullptr)
    if (!M_key_compare(S_key(x), k))
      y = x, x = S_left(x);
    else
      x = S_right(x);

  return const_iterator(y);
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::upper_bound(const Key& k) {
  Link_type y = M_header; /* Last node which is greater than k. */
  Link_type x = M_root(); /* Current node. */

  while (x != nullptr)
    if (M_key_compare(k, S_key(x)))
      y = x, x = S_left(x);
    else
      x = S_right(x);

  return iterator(y);
}

template <class Key, class Value, class KeyOfValue, class Compare>
typename Rb_tree<Key, Value, KeyOfValue, Compare>::const_iterator
Rb_tree<Key, Value, KeyOfValue, Compare>::upper_bound(const Key& k) const {
  Link_type y = M_header; /* Last node which is greater than k. */
  Link_type x = M_root(); /* Current node. */

  while (x != nullptr)
    if (M_key_compare(k, S_key(x)))
      y = x, x = S_left(x);
    else
      x = S_right(x);

  return const_iterator(y);
}

template <class Key, class Value, class KeyOfValue, class Compare>
inline std::pair<typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator,
    typename Rb_tree<Key, Value, KeyOfValue, Compare>::iterator>
Rb_tree<Key, Value, KeyOfValue, Compare>::equal_range(const Key& k) {
  return std::pair<iterator, iterator>(lower_bound(k), upper_bound(k));
}

template <class Key, class Value, class KeyOfValue, class Compare>
inline std::pair<
    typename Rb_tree<Key, Value, KeyOfValue, Compare>::const_iterator,
    typename Rb_tree<Key, Value, KeyOfValue, Compare>::const_iterator>
Rb_tree<Key, Value, KeyOfValue, Compare>::equal_range(const Key& k) const {
  return std::pair<const_iterator, const_iterator>(
      lower_bound(k), upper_bound(k));
}

template <class Key, class Value, class KeyOfValue,
    class Compare = std::less<Key>>
struct rb_tree : public Rb_tree<Key, Value, KeyOfValue, Compare> {
  typedef Rb_tree<Key, Value, KeyOfValue, Compare> Base;

  rb_tree(const Compare& comp = Compare()) : Base(comp) {}

  rb_tree(rb_tree const&) = default;
  ~rb_tree() = default;
};

template <class T>
struct Rb_value_as_key {
  T const& operator()(T const& value) const { return value; }
};

template <class Key>
using set = rb_tree<Key, Key, Rb_value_as_key<Key>>;

template <class T1, class T2>
struct Rb_first_as_key {
  T1 const& operator()(std::pair<T1, T2> const& value) const {
    return value.first;
  }
};

template <class Key, class Value>
struct map
    : public rb_tree<Key, std::pair<Key, Value>, Rb_first_as_key<Key, Value>> {
  Value& operator[](const Key& key) {
    auto it = this->upper_bound(key);
    if (it == this->end() || (!(it->first == key))) {
      it = this->insert(it, std::make_pair(key, Value()));
    }
    return it->second;
  }
};

}  // namespace Omega_h

#endif
