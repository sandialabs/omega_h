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
#ifndef R123_SEE_H
#define R123_SEE_H

#if R123_USE_SSE

#if R123_USE_X86INTRIN_H
#include <x86intrin.h>
#endif
#if R123_USE_IA32INTRIN_H
#include <ia32intrin.h>
#endif
#if R123_USE_XMMINTRIN_H
#include <xmmintrin.h>
#endif
#if R123_USE_EMMINTRIN_H
#include <emmintrin.h>
#endif
#if R123_USE_SMMINTRIN_H
#include <smmintrin.h>
#endif
#if R123_USE_WMMINTRIN_H
#include <wmmintrin.h>
#endif
#if R123_USE_INTRIN_H
#include <intrin.h>
#endif
#ifdef __cplusplus
#include <iostream>
#include <limits>
#include <stdexcept>
#endif

#if R123_USE_ASM_GNU

/* bit25 of CX tells us whether AES is enabled. */
R123_STATIC_INLINE int haveAESNI(){
    unsigned int eax, ebx, ecx, edx;
    __asm__ __volatile__ ("cpuid": "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) :
                      "a" (1));
    return (ecx>>25) & 1;
}
#elif R123_USE_CPUID_MSVC
R123_STATIC_INLINE int haveAESNI(){
    int CPUInfo[4];
    __cpuid(CPUInfo, 1);
    return (CPUInfo[2]>>25)&1;
}
#else /* R123_USE_CPUID_??? */
#warning "No R123_USE_CPUID_XXX method chosen.  haveAESNI will always return false"
R123_STATIC_INLINE int haveAESNI(){
    return 0;
}
#endif /* R123_USE_ASM_GNU || R123_USE_CPUID_MSVC */

// There is a lot of annoying and inexplicable variation in the
// SSE intrinsics available in different compilation environments.
// The details seem to depend on the compiler, the version and
// the target architecture.  Rather than insisting on
// R123_USE_feature tests for each of these in each of the
// compilerfeatures.h files we just keep the complexity localized
// to here...
#if (defined(__ICC) && __ICC<1210) || (defined(_MSC_VER) && !defined(_WIN64))
/* Is there an intrinsic to assemble an __m128i from two 64-bit words? 
   If not, use the 4x32-bit intrisic instead.  N.B.  It looks like Intel
   added _mm_set_epi64x to icc version 12.1 in Jan 2012.
*/
R123_STATIC_INLINE __m128i _mm_set_epi64x(std::uint64_t v1, std::uint64_t v0){
    union{
        std::uint64_t u64;
        std::uint32_t u32[2];
    } u1, u0;
    u1.u64 = v1;
    u0.u64 = v0;
    return _mm_set_epi32(u1.u32[1], u1.u32[0], u0.u32[1], u0.u32[0]);
}
#endif
/* _mm_extract_lo64 abstracts the task of extracting the low 64-bit
   word from an __m128i.  The _mm_cvtsi128_si64 intrinsic does the job
   on 64-bit platforms.  Unfortunately, both MSVC and Open64 fail
   assertions in ut_M128.cpp and ut_carray.cpp when we use the
   _mm_cvtsi128_si64 intrinsic.  (See
   https://bugs.open64.net/show_bug.cgi?id=873 for the Open64 bug).
   On 32-bit platforms, there's no MOVQ, so there's no intrinsic.
   Finally, even if the intrinsic exists, it may be spelled with or
   without the 'x'.
*/
#if !defined(__x86_64__) || defined(_MSC_VER) || defined(__OPEN64__)
R123_STATIC_INLINE std::uint64_t _mm_extract_lo64(__m128i si){
    union{
        std::uint64_t u64[2];
        __m128i m;
    }u;
    _mm_store_si128(&u.m, si);
    return u.u64[0];
}
#elif defined(__llvm__) || defined(__ICC)
R123_STATIC_INLINE std::uint64_t _mm_extract_lo64(__m128i si){
    return std::uint64_t(_mm_cvtsi128_si64(si));
}
#else /* GNUC, others */
/* FWIW, gcc's emmintrin.h has had the 'x' spelling
   since at least gcc-3.4.4.  The no-'x' spelling showed up
   around 4.2. */
R123_STATIC_INLINE std::uint64_t _mm_extract_lo64(__m128i si){
    return std::uint64_t(_mm_cvtsi128_si64x(si));
}
#endif
#if defined(__GNUC__) && __GNUC__ < 4
/* the cast builtins showed up in gcc4. */
R123_STATIC_INLINE __m128 _mm_castsi128_ps(__m128i si){
    return (__m128)si;
}
#endif

#else /* !R123_USE_SSE */
R123_STATIC_INLINE int haveAESNI(){
    return 0;
}
#endif /* R123_USE_SSE */

#endif /* _Random123_sse_dot_h__ */
