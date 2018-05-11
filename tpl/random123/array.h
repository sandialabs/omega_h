/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef R123_ARRAY_H
#define R123_ARRAY_H

#include <random123/features/compilerfeatures.h>
#include <random123/features/sse.h>

#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <iterator>
#include <limits>
#include <iostream>

/** @defgroup arrayNxW The r123arrayNxW classes 

    Each of the r123arrayNxW is a fixed size array of N W-bit unsigned integers.
    It is functionally equivalent to the C++0x std::array<N, std::uintW_t>,
    but does not require C++0x features or libraries.

    In addition to meeting most of the requirements of a Container,
    it also has a member function, incr(), which increments the zero-th
    element and carrys overflows into higher indexed elements.  Thus,
    by using incr(), sequences of up to 2^(N*W) distinct values
    can be produced. 
*/

template <typename value_type>
inline R123_CUDA_DEVICE value_type assemble_from_u32(std::uint32_t *p32){
    value_type v=0;
    for(size_t i=0; i<(3+sizeof(value_type))/4; ++i)
        v |= (value_type(*p32++)) << (32*i);
    return v;
}

// Work-alike methods and typedefs modeled on std::array:
#define R123_CXXMETHODS(_N, W, T)                                            \
    typedef T value_type;                                               \
    typedef T* iterator;                                                \
    typedef const T* const_iterator;                                    \
    typedef value_type& reference;                                      \
    typedef const value_type& const_reference;                          \
    typedef size_t size_type;                                           \
    typedef ptrdiff_t difference_type;                                  \
    typedef T* pointer;                                                 \
    typedef const T* const_pointer;                                     \
    /* Boost.array has static_size.  C++11 specializes tuple_size */    \
    enum {static_size = _N};                                            \
    R123_CUDA_DEVICE reference operator[](size_type i){return v[i];}                     \
    R123_CUDA_DEVICE const_reference operator[](size_type i) const {return v[i];}        \
    R123_CUDA_DEVICE reference at(size_type i){ if(i >=  _N) R123_THROW(std::out_of_range("array index out of range")); return (*this)[i]; } \
    R123_CUDA_DEVICE const_reference at(size_type i) const { if(i >=  _N) R123_THROW(std::out_of_range("array index out of range")); return (*this)[i]; } \
    R123_CUDA_DEVICE size_type size() const { return  _N; }                              \
    R123_CUDA_DEVICE size_type max_size() const { return _N; }                           \
    R123_CUDA_DEVICE bool empty() const { return _N==0; }                                \
    R123_CUDA_DEVICE iterator begin() { return &v[0]; }                                  \
    R123_CUDA_DEVICE iterator end() { return &v[_N]; }                                   \
    R123_CUDA_DEVICE const_iterator begin() const { return &v[0]; }                      \
    R123_CUDA_DEVICE const_iterator end() const { return &v[_N]; }                       \
    R123_CUDA_DEVICE const_iterator cbegin() const { return &v[0]; }                     \
    R123_CUDA_DEVICE const_iterator cend() const { return &v[_N]; }                      \
    R123_CUDA_DEVICE pointer data(){ return &v[0]; }                                     \
    R123_CUDA_DEVICE const_pointer data() const{ return &v[0]; }                         \
    R123_CUDA_DEVICE reference front(){ return v[0]; }                                   \
    R123_CUDA_DEVICE const_reference front() const{ return v[0]; }                       \
    R123_CUDA_DEVICE reference back(){ return v[_N-1]; }                                 \
    R123_CUDA_DEVICE const_reference back() const{ return v[_N-1]; }                     \
    R123_CUDA_DEVICE bool operator==(const r123array##_N##x##W& rhs) const{ \
      /* CUDA3 does not have std::equal */ \
      for (size_t i = 0; i < _N; ++i) \
        if (v[i] != rhs.v[i]) return false; \
      return true; \
    } \
    R123_CUDA_DEVICE bool operator!=(const r123array##_N##x##W& rhs) const{ return !(*this == rhs); } \
    /* CUDA3 does not have std::fill_n */ \
    R123_CUDA_DEVICE void fill(const value_type& val){ for (size_t i = 0; i < _N; ++i) v[i] = val; } \
    R123_CUDA_DEVICE void swap(r123array##_N##x##W& rhs){ \
      /* CUDA3 does not have std::swap_ranges */ \
      for (size_t i = 0; i < _N; ++i) { \
          T tmp = v[i]; \
          v[i] = rhs.v[i]; \
          rhs.v[i] = tmp; \
      } \
    } \
    /* seed(SeedSeq) would be a constructor if having a constructor */  \
    /* didn't cause headaches with defaults */                          \
    template <typename SeedSeq>                                         \
    R123_CUDA_DEVICE static r123array##_N##x##W seed(SeedSeq &ss){      \
        r123array##_N##x##W ret;                                        \
        const size_t Ngen = _N*((3+sizeof(value_type))/4);              \
        std::uint32_t u32[Ngen];                                             \
        std::uint32_t *p32 = &u32[0];                                        \
        ss.generate(&u32[0], &u32[Ngen]);                               \
        for(size_t i=0; i<_N; ++i){                                     \
            ret.v[i] = assemble_from_u32<value_type>(p32);              \
            p32 += (3+sizeof(value_type))/4;                            \
        }                                                               \
        return ret;                                                     \
    }                                                                   \
protected:                                                              \
    
                                                                        
// There are several tricky considerations for the insertion and extraction
// operators:
// - we would like to be able to print r123array16x8 as a sequence of 16 integers,
//   not as 16 bytes.
// Solution: r123arrayinsertable and r123arrayextractable

template<typename T>
struct r123arrayinsertable{
    const T& v;
    r123arrayinsertable(const T& t_) : v(t_) {} 
    friend std::ostream& operator<<(std::ostream& os, const r123arrayinsertable<T>& t){
        return os << t.v;
    }
};

template<>
struct r123arrayinsertable<std::uint8_t>{
    const std::uint8_t& v;
    r123arrayinsertable(const std::uint8_t& t_) : v(t_) {} 
    friend std::ostream& operator<<(std::ostream& os, const r123arrayinsertable<std::uint8_t>& t){
        return os << int(t.v);
    }
};

template<typename T>
struct r123arrayextractable{
    T& v;
    r123arrayextractable(T& t_) : v(t_) {}
    friend std::istream& operator>>(std::istream& is, r123arrayextractable<T>& t){
        return is >> t.v;
    }
};

template<>
struct r123arrayextractable<std::uint8_t>{
    std::uint8_t& v;
    r123arrayextractable(std::uint8_t& t_) : v(t_) {} 
    friend std::istream& operator>>(std::istream& is, r123arrayextractable<std::uint8_t>& t){
        int i;
        is >>  i;
        t.v = std::uint8_t(i);
        return is;
    }
};

#define R123_CXXOVERLOADS(_N, W, T)                                          \
                                                                        \
inline std::ostream& operator<<(std::ostream& os, const r123array##_N##x##W& a){   \
    os << r123arrayinsertable<T>(a.v[0]);                                  \
    for(size_t i=1; i<_N; ++i)                                          \
        os << " " << r123arrayinsertable<T>(a.v[i]);                       \
    return os;                                                          \
}                                                                       \
                                                                        \
inline std::istream& operator>>(std::istream& is, r123array##_N##x##W& a){         \
    for(size_t i=0; i<_N; ++i){                                         \
        r123arrayextractable<T> x(a.v[i]);                                 \
        is >> x;                                                        \
    }                                                                   \
    return is;                                                          \
}                                                                       \
                                                                        \
namespace r123{                                                        \
 typedef r123array##_N##x##W Array##_N##x##W;                          \
}

/* R123_ARRAY_TPL expands to a declaration of struct r123arrayNxW.  

   In C, it's nothing more than a struct containing an array of N
   objects of type T.

   In C++ it's the same, but endowed with an assortment of member
   functions, typedefs and friends.  In C++, r123arrayNxW looks a lot
   like std::array<T,N>, has most of the capabilities of a container,
   and satisfies the requirements outlined in compat/Engine.hpp for
   counter and key types.  ArrayNxW, in the r123 namespace is
   a typedef equivalent to r123arrayNxW.
*/

#define R123_ARRAY_TPL(_N, W, T)                   \
    /** @ingroup arrayNxW */                        \
    /** @see arrayNxW */                            \
struct r123array##_N##x##W{                         \
 T v[_N];                                       \
 R123_CXXMETHODS(_N, W, T)                           \
};                                              \
                                                \
R123_CXXOVERLOADS(_N, W, T)

/** @endcond */

R123_ARRAY_TPL(1, 32, std::uint32_t)  /* r123array1x32 */
R123_ARRAY_TPL(2, 32, std::uint32_t)  /* r123array2x32 */
R123_ARRAY_TPL(4, 32, std::uint32_t)  /* r123array4x32 */
R123_ARRAY_TPL(8, 32, std::uint32_t)  /* r123array8x32 */

R123_ARRAY_TPL(1, 64, std::uint64_t)  /* r123array1x64 */
R123_ARRAY_TPL(2, 64, std::uint64_t)  /* r123array2x64 */
R123_ARRAY_TPL(4, 64, std::uint64_t)  /* r123array4x64 */

R123_ARRAY_TPL(16, 8, std::uint8_t)  /* r123array16x8 for ARSsw, AESsw */

/* In C++, it's natural to use sizeof(a::value_type), but in C it's
   pretty convoluted to figure out the width of the value_type of an
   r123arrayNxW:
*/
#define R123_W(a)   (8*sizeof(((a *)0)->v[0]))

/** @namespace r123
  Most of the Random123 C++ API is contained in the r123 namespace. 
*/

#endif
