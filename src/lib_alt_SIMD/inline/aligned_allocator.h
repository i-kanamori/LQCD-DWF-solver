/*!

        @file    aligned_allocator.h

        @brief   allocator with alignment

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================

#ifndef __ALINGED_ALLOCATOR_SIMD_H__
#define __ALINGED_ALLOCATOR_SIMD_H__

#include <stdlib.h>
#include <cstddef>

#ifdef HHHH_256_HHHH
#define ALIGNMENT 32
#endif

#ifdef HHHH_512_HHHH
#define ALIGNMENT 64
#endif


/*
 cf.
  Programing Language C++ 4th eidtion
  Bjarne Stroustrup
  sec. 34.4

see also
  allocator.h in stdc++ lib

 */
template<typename _Tp>
struct aligned_allocator{
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  aligned_allocator() throw() {}
  ~aligned_allocator() throw() {}

  template<typename _Tp1>
  struct rebind
  { typedef aligned_allocator<_Tp1> other; };

  _Tp* allocate(size_t n){
    void* p=0;
    void **memptr=&p;
    int err=posix_memalign(memptr, ALIGNMENT, sizeof(_Tp) * n);
    if(err){
      throw std::bad_alloc();
    }
    return reinterpret_cast<_Tp*>(p);
  }

  void deallocate(_Tp* p, size_t){
    free(p);
  }
};

  template<typename _T1, typename _T2>
    inline bool
    operator==(const aligned_allocator<_T1>&, const aligned_allocator<_T2>&)
     { return true; }

   template<typename _Tp>
     inline bool
     operator==(const aligned_allocator<_Tp>&, const aligned_allocator<_Tp>&)
     { return true; }

   template<typename _T1, typename _T2>
     inline bool
     operator!=(const aligned_allocator<_T1>&, const aligned_allocator<_T2>&)
     { return false; }

   template<typename _Tp>
     inline bool
     operator!=(const aligned_allocator<_Tp>&, const aligned_allocator<_Tp>&)
     { return false; }

#undef ALIGNMENT
#endif

