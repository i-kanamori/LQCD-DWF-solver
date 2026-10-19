/*!

        @file    vector_arith-inc.h

        @brief   arithmetics of simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
#ifndef VECTOR_ARITH_INC_SIMD_H
#define VECTOR_ARITH_INC_SIMD_H

#include "valtype.h"


#ifdef HHHH_256_HHHH
#include "vector_arith_256-inc.h"
#elif defined HHHH_512_HHHH
#include "vector_arith_512-inc.h"
#else
#error("unknown simd type\n")
#endif


namespace { // common to implementations

  /*
  static inline vtype_realD operator+(const vtype_ReadD &va, const vtype_RealD &vb){
    //printf("hoge: operator+ for vrealD_t\n");
    return add(va,vb);
  }
  static inline vtype_complexD operator+(const vtype_complexD &va, const vtype_complexD &vb){
    return add(va,vb);
  }

    static inline vtype_realF operator+(const vtype_ReadF &va, const vtype_RealF &vb){
    //printf("hoge: operator+ for vrealD_t\n");
    return add(va,vb);
  }
  static inline vtype_complexF operator+(const vtype_complexF &va, const vtype_complexF &vb){
    return add(va,vb);
  }
  */


  static inline vrealD_t operator+(const vrealD_t &va, const vrealD_t &vb){
    //printf("hoge: operator+ for vrealD_t\n");
    return add(va,vb);
  }
  static inline vcomplexD_t operator+(const vcomplexD_t &va, const vcomplexD_t &vb){
    return add(va,vb);
  }
  inline void operator+=(vcomplexD_t &va,const vcomplexD_t &vb){
    accum_add(va,vb);
  }
  static inline vrealD_t operator-(const vrealD_t &va, const vrealD_t &vb){
    return sub(va,vb);
  }
  static inline vcomplexD_t operator-(const vcomplexD_t &va, const vcomplexD_t &vb){
    return sub(va,vb);
  }
  inline void operator-=(vcomplexD_t &va,const vcomplexD_t &vb){
    accum_sub(va,vb);
  }


  static inline vrealD_t operator*(const vrealD_t &va, const vrealD_t &vb){
    //printf("hoge: operator+ for vrealD_t\n");
    return mult(va,vb);
  }
  static inline vcomplexD_t operator*(const vcomplexD_t &va, const vcomplexD_t &vb){
    return cmult(va,vb);
  }


  static inline vrealF_t operator+(const vrealF_t &va, const vrealF_t &vb){
    return add(va,vb);
  }
  static inline vcomplexF_t operator+(const vcomplexF_t &va, const vcomplexF_t &vb){
    return add(va,vb);
  }
  inline void operator+=(vcomplexF_t &va,const vcomplexF_t &vb){
    accum_add(va,vb);
  }
  static inline vrealF_t operator-(const vrealF_t &va, const vrealF_t &vb){
    return sub(va,vb);
  }
  static inline vcomplexF_t operator-(const vcomplexF_t &va, const vcomplexF_t &vb){
    return sub(va,vb);
  }
  inline void operator-=(vcomplexF_t &va,const vcomplexF_t &vb){
    accum_sub(va,vb);
  }

  static inline vrealF_t operator*(const vrealF_t &va, const vrealF_t &vb){
    return mult(va,vb);
  }
  static inline vcomplexF_t operator*(const vcomplexF_t &va, const vcomplexF_t &vb){
    return cmult(va,vb);
  }

}

#endif
