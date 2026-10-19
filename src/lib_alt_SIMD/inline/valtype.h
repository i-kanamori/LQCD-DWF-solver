/*!
        @file    $Id: valtype.h  #$

        @brief   alias to varables simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */

#ifndef VALTYPE_H
#define VALTYPE_H

  //helpers
#ifdef HHHH_512_HHHH
#undef HHHH_512_HHHH
#endif
#ifdef HHHH_256_HHHH
#undef HHHH_256_HHHH
#endif

#define HHHH_128_HHHH
#if defined(AVX512) || defined(KNC) || defined(IMCI)
#define HHHH_512_HHHH
#undef HHHH_128_HHHH
#endif

#if defined(AVX1) || defined(AVX2) || defined(AVXFMA4)
#define HHHH_256_HHHH
#undef HHHH_128_HHHH
#endif

#ifdef HHHH_128_HHHH
#error("simd length 128 is not supprted\n")
#endif


//////////////////////////////////////////////
//
// headers
//
//////////////////////////////////////////////


#ifdef HHHH_512_HHHH
#include "valtype_512.h"
#elif defined HHHH_256_HHHH
#include "valtype_256.h"
#else
#error("unknown simd type\n")
#endif


// complex number
#include "bridge_complex.h"

// alinged allocator
#include "aligned_allocator.h"

// std library
#include  <vector>


//////////////////////////////////////////////
//
// wrappers
//
//////////////////////////////////////////////

// wrapper struct for vetor varaibles
template<typename ValType> struct ValType_real{
};
template<typename ValType> struct ValType_complex{
};


// aligned vector
template<class T>
using aligned_vector = std::vector<T,aligned_allocator<T> >;


//////////////////////////////////////////////
//
// body of the templated wrappers
//
//////////////////////////////////////////////


template<> struct ValType_complex<double>{
  using vector_type=Double::Vsimd_t;
  using scalar_type=std::complex<double>;
  vector_type vec;
  inline ValType_complex<double>(const vector_type &vv){vec=vv;}

  inline ValType_complex<double>(){}
  static constexpr int Nsimd() {return Double::Nsimd;}
};
typedef ValType_complex<double> vcomplexD_t;
typedef ValType_complex<double>::scalar_type complexD_t;

template<> struct ValType_real<double>
{
  using vector_type=Double::Vsimd_t_r;
  using scalar_type=double;
  vector_type vec;
  inline ValType_real<double>(const vector_type &vv){vec=vv;}
  inline ValType_real<double>(){}
  inline ValType_real<double>(const ValType_complex<double> &vv){vec=vv.vec;}
  static constexpr int Nsimd() {return Double::Nsimd*2;}
};
typedef ValType_real<double> vrealD_t;
typedef ValType_real<double>::scalar_type realD_t;


template<> struct ValType_complex<float>{
  using vector_type=Float::Vsimd_t;
  using scalar_type=std::complex<float>;
  vector_type vec;
  inline ValType_complex<float>(const vector_type &vv){vec=vv;}
  inline ValType_complex<float>(){}
  static constexpr int Nsimd() {return Float::Nsimd;}
};
typedef ValType_complex<float> vcomplexF_t;
typedef ValType_complex<float>::scalar_type complexF_t;



template<> struct ValType_real<float>
{
  using vector_type=Float::Vsimd_t_r;
  using scalar_type=float;
  vector_type vec;
  inline ValType_real<float>(const vector_type &vv){vec=vv;}
  inline ValType_real<float>(const ValType_complex<float> &vv){vec=vv.vec;}
  inline ValType_real<float>(){}
  static constexpr int Nsimd() {return Float::Nsimd*2;}
};
typedef ValType_real<float> vrealF_t;
typedef ValType_real<float>::scalar_type realF_t;



#endif
