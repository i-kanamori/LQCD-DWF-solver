/*!

        @file    vector_shift_double-inc.h

        @brief   shift in simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
#ifndef VECTOR_SHIFT_DOUBLE_INC_H
#define VECTOR_SHIFT_DOUBLE_INC_H

#ifdef IMPLTYPE_ALT_IMPL1

#ifdef HHHH_256_HHHH
#include "vector_shift_double256_impl1-inc.h"
#elif defined HHHH_512_HHHH
#include "vector_shift_double512_impl1-inc.h"
#else
#error("unknown simd type\n")
#endif

#elif defined IMPLTYPE_ALT_GRIDSIMD

#ifdef HHHH_256_HHHH
#include "vector_shift_double256_gridsimd-inc.h"
#elif defined HHHH_512_HHHH
#include "vector_shift_double512_gridsimd-inc.h"
#else
#error("unknown simd type\n")
#endif

#endif

namespace{
  using complex_t = complexD_t;
  using vcomplex_t = vcomplexD_t;
  using real_t = double;
  using vreal_t = vrealD_t;

  void print_simd(const Bridge::VerboseLevel &m_vl){
    if(has_simd_x){
      vout.general(m_vl, "has_simd_x: yes\n");
    } else {
      vout.general(m_vl, "has_simd_x: no\n");
    }
    if(has_simd_y){
      vout.general(m_vl, "has_simd_y: yes\n");
    } else {
      vout.general(m_vl, "has_simd_y: no\n");
    }
    if(has_simd_z){
      vout.general(m_vl, "has_simd_z: yes\n");
    } else {
      vout.general(m_vl, "has_simd_z: no\n");
    }
    if(has_simd_t){
      vout.general(m_vl, "has_simd_t: yes\n");
    } else {
      vout.general(m_vl, "has_simd_t: no\n");
    }
  }
}

#endif  // include guard
