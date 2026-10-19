/*!

        @file    vector_arith_256-inc.h

        @brief   arithmetics of simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================

// some arithmetics
#ifndef VECTOT_ARITH_256_SIMD_H
#define VECTOT_ARITH_256_SIMD_H

/* Immediate value, as Quaternary (4 digits) */
#define IMM_QUART4(A,B,C,D) ((A<<6)|(B<<4)|(C<<2)|(D))
/* Immediate value, as Quaternary (2 digits) */
#define IMM_QUART2(A,B) ((A<<2)|(B))
/* Immediate value, as Binary (8 digits) */
#define IMM_BIN8(A,B,C,D,E,F,G,H) ((A<<7)|(B<<6)|(C<<5)|(D<<4)|(E<<3)|(F<<2)|(G<<1)|H)
/* Immediate value, as Binary (4 digits) */
#define IMM_BIN4(A,B,C,D) ((A<<3)|(B<<2)|(C<<1)|D)


namespace {

  ///////////////////////////////////////
  // double
  //////////////////////////////////////

  // broadcast
  static inline void bcast(vrealD_t &va, const realD_t &a){
    va.vec=_mm256_set_pd(a,a,a,a);
  }
  static inline void bcast(vcomplexD_t &va, const complexD_t &a){
    double re=real(a);
    double im=imag(a);
    va.vec=_mm256_set_pd(im,re,im,re);
  }

  // streaming store
  static inline void sstore(vcomplexD_t *vbuf, const vcomplexD_t &va){
    _mm256_stream_pd( (double*) &(vbuf->vec), va.vec);
  }

  static inline void sstore(vrealD_t *vbuf, const vrealD_t &va){
    _mm256_stream_pd( (double*) &(vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealD_t &va){
    va.vec=_mm256_setzero_pd();
  }
  static inline void set_zero(vcomplexD_t &va){
    va.vec=_mm256_setzero_pd();
  }

  // complex conjugate
  static inline vcomplexD_t conj(const vcomplexD_t &in){
    __m256d tmp = _mm256_shuffle_pd(in.vec, in.vec, IMM_BIN4(0,1,0,1)); // re,im --> im, re
    tmp=_mm256_addsub_pd(_mm256_setzero_pd(), tmp); // -im, re
    return vcomplexD_t(_mm256_shuffle_pd(tmp, tmp, IMM_BIN4(0,1,0,1))); // re, -im
  }

  // reduction
  static inline realD_t reduce(const vrealD_t &va){
    __m256d tmp;
    double *sum=(double*)&tmp;

    tmp=(__m256d)_mm256_permute2f128_ps((__m256)va.vec, (__m256)va.vec, IMM_QUART4(0,0,0,1));  // 12 34 --> 34 12
    tmp=_mm256_add_pd(tmp, va.vec);  // (12+34 34+12) = ( 1+3 2+4 3+1 4+2)
    tmp=_mm256_hadd_pd(tmp,tmp);       // (1+3+2+4 1+3+2+4 3+1+4+2 3+1+4+2)
    return sum[0];
  };

  static inline complexD_t reduce(const vcomplexD_t &va){
    __m256d tmp;
    complexD_t *sum=(complexD_t*)&tmp;

    tmp=(__m256d)_mm256_permute2f128_ps((__m256)va.vec, (__m256)va.vec, IMM_QUART4(0,0,0,1));  // 12 34 --> 34 12
    tmp=_mm256_add_pd(tmp, va.vec);  // (12+34 34+12) = ( 1+3 2+4 3+1 4+2)
    return sum[0];

  }

  // add
  static inline vrealD_t add(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(_mm256_add_pd(va.vec,vb.vec));
  }
  static inline void accum_add(vrealD_t &va,const vrealD_t &vb){
    va.vec=_mm256_add_pd(va.vec,vb.vec);
  }
  static inline vcomplexD_t add(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(_mm256_add_pd(va.vec,vb.vec));
  }
  static inline void accum_add(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec=_mm256_add_pd(va.vec,vb.vec);
  }

  // subtraction
  static inline vrealD_t sub(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(_mm256_sub_pd(va.vec,vb.vec));
  }
  static inline void accum_sub(vrealD_t &va,const vrealD_t &vb){
    va.vec=_mm256_sub_pd(va.vec, vb.vec);
  }
  static inline vcomplexD_t sub(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(_mm256_sub_pd(va.vec,vb.vec));
  }
  static inline void accum_sub(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec=_mm256_sub_pd(va.vec, vb.vec);
  }

  // reciprocal number
  static inline vrealD_t rcp(const vrealD_t &va){
    __m256d tmp=_mm256_set_pd(1.0,1.0,1.0,1.0);
    return vrealD_t(_mm256_div_pd(tmp, va.vec));
  }

  // mult
  static inline vrealD_t mult(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t( _mm256_mul_pd(va.vec, vb.vec));
  }

  // mult
  static inline vcomplexD_t mult(const vrealD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t( _mm256_mul_pd(va.vec, vb.vec));
  }


#if defined (AVX2)
  static inline vcomplexD_t cmult(const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    return vcomplexD_t( _mm256_fmaddsub_pd( a_real, vb.vec, a_imag ));
  }

  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    out.vec=_mm256_fmaddsub_pd( a_real, vb.vec, a_imag );
  }

  static inline void accum_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    a_real=_mm256_fmaddsub_pd( a_real, vb.vec, a_imag );
    out.vec=_mm256_add_pd(out.vec, a_real);
  }

  static inline void accum_sub_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    a_real=_mm256_fmaddsub_pd( a_real, vb.vec, a_imag );
    out.vec=_mm256_sub_pd(out.vec, a_real);
  }

#else
  static inline vcomplexD_t cmult(const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    return vcomplexD_t( _mm256_addsub_pd(tmp,a_imag));
  }

  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    out.vec= _mm256_addsub_pd(tmp,a_imag);
  }

  static inline void accum_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    tmp = _mm256_addsub_pd(tmp,a_imag);
    out.vec = _mm256_add_pd(out.vec,tmp);
  }

  static inline void accum_sub_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    tmp = _mm256_addsub_pd(tmp,a_imag);
    out.vec = _mm256_sub_pd(out.vec,tmp);
  }

#endif


  // mult:  conj(a)*b
#if defined (AVX2)
  static inline vcomplexD_t cmult_conj(const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    return vcomplexD_t(_mm256_fmsubadd_pd( a_real, vb.vec, a_imag ));
  }

  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    // slight modification of a*b in Grid
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    out.vec =_mm256_fmsubadd_pd( a_real, vb.vec, a_imag );
  }

  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_fmsubadd_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) , out.vec);
    out.vec =_mm256_fmsubadd_pd( a_real, vb.vec, a_imag );
  }

  static inline void accum_sub_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_real = _mm256_sub_pd(_mm256_setzero_pd(), a_real);
    a_imag = _mm256_fmaddsub_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) , out.vec);
    out.vec =_mm256_fmaddsub_pd( a_real, vb.vec, a_imag );
  }


#else
  static inline vcomplexD_t cmult_conj(const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    a_imag = _mm256_sub_pd( _mm256_setzero_pd(), a_imag);
    return vcomplexD_t(_mm256_addsub_pd(tmp,a_imag));
  }

  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    a_imag = _mm256_mul_pd( a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    a_imag = _mm256_sub_pd( _mm256_setzero_pd(), a_imag);
    out.vec= _mm256_addsub_pd(tmp,a_imag);
  }

  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    out.vec = _mm256_add_pd(out.vec, tmp);
    tmp = _mm256_mul_pd(a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    tmp = _mm256_sub_pd( _mm256_setzero_pd(), tmp);
    out.vec = _mm256_addsub_pd(out.vec, tmp);
  }


  static inline void accum_sub_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m256d a_real = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(0,0,0,0) );
    __m256d a_imag = _mm256_shuffle_pd( va.vec, va.vec, IMM_BIN4(1,1,1,1) );
    __m256d tmp = _mm256_mul_pd(a_real, vb.vec);
    out.vec = _mm256_sub_pd(out.vec, tmp);
    tmp = _mm256_mul_pd(a_imag, _mm256_permute_pd( vb.vec, IMM_BIN4(0,1,0,1) ) );
    out.vec = _mm256_addsub_pd(out.vec, tmp);
  }
#endif


  // sqrt(-1)
  static inline vcomplexD_t TimesI(const vcomplexD_t &in){
    __m256d tmp = _mm256_shuffle_pd(in.vec, in.vec, IMM_BIN4(0,1,0,1));  // re,im --> im, re
    return vcomplexD_t(_mm256_addsub_pd(_mm256_setzero_pd(), tmp));      // 0-im, 0+re
  }

  static void inline accum_TimesI(vcomplexD_t &out, const vcomplexD_t &in){
    __m256d tmp = _mm256_shuffle_pd(in.vec, in.vec, IMM_BIN4(0,1,0,1)); // re,im --> im, re
    out.vec=_mm256_addsub_pd(out.vec, tmp);                             // out.re-im, out.im+re
  }

  static inline vcomplexD_t TimesMinusI(const vcomplexD_t &in){
    __m256d tmp = _mm256_shuffle_pd(in.vec, in.vec, IMM_BIN4(0,1,0,1)); // re,im --> im, re
    tmp = _mm256_sub_pd(_mm256_setzero_pd(), tmp);                      // -im, -re
    return vcomplexD_t(_mm256_addsub_pd(_mm256_setzero_pd(), tmp));     // 0-(-im), 0+(-re)
  }

  static inline void accum_TimesMinusI(vcomplexD_t &out, const vcomplexD_t &in){
    __m256d tmp = _mm256_shuffle_pd(in.vec, in.vec, IMM_BIN4(0,1,0,1)); // re,im --> im, re
    tmp = _mm256_sub_pd(_mm256_setzero_pd(), tmp);                      // -im, -re
    out.vec=_mm256_addsub_pd(out.vec, tmp);                             // out.re-(-im), out.im+(-re)
  }

  static inline vcomplexD_t Times_pm(const vcomplexD_t &in){
    // (c1,c2) --> (c1, -c2)
    __m256d tmp = _mm256_permute4x64_pd(in.vec, IMM_QUART4(0,2,1,3));
    tmp = _mm256_addsub_pd(_mm256_setzero_pd(), tmp);
    return vcomplexD_t(_mm256_permute4x64_pd(tmp,IMM_QUART4(0,2,1,3)));
  }

  static inline vcomplexD_t Times_mp(const vcomplexD_t &in){
    // (c1,c2) --> (-c1, c2)
    __m256d tmp = _mm256_permute4x64_pd(in.vec, IMM_QUART4(3,1,2,0));
    tmp = _mm256_addsub_pd(_mm256_setzero_pd(), tmp);
    return vcomplexD_t(_mm256_permute4x64_pd(tmp,IMM_QUART4(3,1,2,0)));
  }

  // axpy: y := ax + y
  static inline void axpy(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    accum_cmult(vy,va,vx);
  }

#ifdef AVX2
  static inline void axpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec=_mm256_fmadd_pd(va.vec, vx.vec, vy.vec);
  }
  static inline void axpy(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec=_mm256_fmadd_pd(va.vec, vx.vec, vy.vec);
  }
#else
  static inline void axpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    __m256d tmp = _mm256_mul_pd(va.vec, vx.vec);
    vy.vec=_mm256_add_pd(vy.vec, tmp);
  }
  static inline void axpy(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    __m256d tmp = _mm256_mul_pd(va.vec, vx.vec);
    vy.vec=_mm256_add_pd(vy.vec, tmp);
  }
#endif

  // alias
  static inline void vaxpy(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    axpy(vy,va,vx);
  }
  static inline void vaxpy(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    axpy(vy, va, vx);
  }
  static inline void vaxpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    axpy(vy, va, vx);
  }
  static inline void accum_abs2(vrealD_t &vy, const vcomplexD_t &vx){
    axpy(vy, vrealD_t(vx.vec), vrealD_t(vx.vec));
  }


  // y := y-ax
#ifdef AVX2
  static inline void accum_sub_mult(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec=_mm256_fnmadd_pd(va.vec, vx.vec, vy.vec);
  }
  static inline void accum_sub_mult(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec=_mm256_fnmadd_pd(va.vec, vx.vec, vy.vec);
  }
#else
  static inline void accum_sub_mult(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    __m256d tmp = _mm256_mul_pd(va.vec, vx.vec);
    vy.vec=_mm256_sub_pd(vy.vec, tmp);
  }
  static inline void accum_sub_mult(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    __m256d tmp = _mm256_mul_pd(va.vec, vx.vec);
    vy.vec=_mm256_sub_pd(vy.vec, tmp);
  }
#endif



  // aypx: y := a*y+x
#ifdef AVX2
  static inline void aypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec = _mm256_fmadd_pd(va.vec, vy.vec, vx.vec);
  }
  static inline void aypx(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec = _mm256_fmadd_pd(va.vec, vy.vec, vx.vec);
  }
#else
  static inline void aypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec = _mm256_mul_pd(va.vec, vy.vec);
    vy.vec = _mm256_add_pd(vy.vec, vx.vec);
  }
  static inline void aypx(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec = _mm256_mul_pd(va.vec, vy.vec);
    vy.vec = _mm256_add_pd(vy.vec, vx.vec);
  }
#endif
  static inline void aypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    set_cmult(vy, vy, va);
    vy.vec=_mm256_add_pd(vy.vec, vx.vec);
  }

  // alias
  static inline void vaypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    aypx(vy, va, vx);
  }
  static inline void vaypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    aypx(vy, va, vx);
  }


  // scal
  static inline void scal(vcomplexD_t &vx, const vcomplexD_t &va){
    set_cmult(vx,va,vx);
  }
  static inline void scal(vrealD_t &vx, const vrealD_t &va){
    vx.vec = _mm256_mul_pd(vx.vec, va.vec);
  }
  static inline void scal(vcomplexD_t &vx, const vrealD_t &va){
    vx.vec = _mm256_mul_pd(vx.vec, va.vec);
  }

  // (re im) (re im)... --> (re re) (re re)...
  static inline void complex2rere(vcomplexD_t &vx){
    vx.vec=_mm256_movedup_pd(vx.vec); // AVX1
    //    vx.vec.v=_mm256_shuffle_pd(vx.vec.v,vx.vec.v,0x0); // AVX1
    // other
    //    double *tmp=(double*)&vx;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }
  }


  // (re im) (re im)... --> (im im) (im im)
  static inline void complex2imim(vcomplexD_t &vx){
    vx.vec=_mm256_shuffle_pd(vx.vec,vx.vec,0xF); // AVX1/AVX2
    // other
    //    double *tmp=(double*)&vx;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp[2*i]=tmp[2*i+1];
    //    }
  }


  // (re im) (re im)... --> (re 0) (re 0)
  static inline void complex2re(vcomplexD_t &vx){
    __m256d tmp=_mm256_set_pd(0.0, 1.0, 0.0, 1.0);
    vx.vec=_mm256_mul_pd(tmp,vx.vec);
  }

  // (re im) (re im)... --> (0 im) (0 im)
  static inline void complex2im(vcomplexD_t &vx){
    __m256d tmp=_mm256_set_pd(1.0, 0.0, 1.0, 0.0);
    vx.vec=_mm256_mul_pd(tmp,vx.vec);
  }

  // TODO
  // (re1 im1),... (re1' im1')... --> (re1 re1') (re2 re2')...
  static inline void rere2complex(vcomplexD_t &vx1, const vcomplexD_t &vx2){
    // other
    //    double *tmp1=(double*)&vx1;
    //    double *tmp2=(double*)&vx2;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp1[2*i+1]=tmp2[2*i];
    //    }

    //  a0 a1 a2 a3, b0 b1 b2 b3 --> a0 b0 a1 b1
    vx1.vec=_mm256_shuffle_pd(vx1.vec, vx2.vec, IMM_BIN4(0,0,0,0));
  }

  // copy as vector varibles
  static inline void copy(vcomplexD_t &va, const vrealD_t &vb){
    va.vec=vb.vec;
  }

  static inline void copy(vrealD_t &va, const vcomplexD_t &vb){
    va.vec=vb.vec;
  }


  ///////////////////////////////////////
  // float
  //////////////////////////////////////

  // broadcast
  static inline void bcast(vrealF_t &va, const realF_t &a){
    va.vec=_mm256_set_ps(a,a,a,a,a,a,a,a);
  }
  static inline void bcast(vcomplexF_t &va, const complexF_t &a){
    const float re=real(a);
    const float im=imag(a);
    va.vec=_mm256_set_ps(im,re,im,re,im,re,im,re);
  }

  // streaming store
  static inline void sstore(vcomplexF_t *vbuf, const vcomplexF_t &va){
    _mm256_stream_ps( (float*) &(vbuf->vec), va.vec);
  }

  static inline void sstore(vrealF_t *vbuf, const vrealF_t &va){
    _mm256_stream_ps( (float*) &(vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealF_t &va){
    va.vec=_mm256_setzero_ps();
  }
  static inline void set_zero(vcomplexF_t &va){
   va.vec=_mm256_setzero_ps();
  }

  // complex conjugate
  static inline vcomplexF_t conj(const vcomplexF_t &in){
    __m256 tmp = _mm256_shuffle_ps(in.vec, in.vec, IMM_QUART4(2,3,0,1));  // re,im --> im, re
    tmp = _mm256_addsub_ps(_mm256_setzero_ps(), tmp);                         // -im, +re
    return vcomplexF_t(_mm256_shuffle_ps(tmp, tmp, IMM_QUART4(2,3,0,1)));           // re, -im
  }

  // reduction
  static inline realF_t reduce(const vrealF_t &va){
    __m256 tmp;
    float *sum=(float*)&tmp;
    tmp=_mm256_permute2f128_ps(va.vec, va.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
    tmp=_mm256_add_ps(tmp, va.vec);  // (1234+5678 5678+1234) = ( 1+5 2+6 3+7 4+8, 1+5 2+6 3+7 4+8)
    __m256 tmp2= _mm256_permute_ps(tmp, IMM_QUART4(1,0,3,2));             // 1234 5678 --> 3456 7812
    tmp=_mm256_add_ps(tmp,tmp2);       // (1+5)+(3+7), (2+6)+(4+8)...
    tmp=_mm256_hadd_ps(tmp,tmp);       // (1+5+3+7)+(2+6+4+8),....
    return sum[0];
  }
  static inline complexF_t reduce(const vcomplexF_t &va){
    __m256 tmp;
    complexF_t *sum=(complexF_t*)&tmp;
    tmp=_mm256_permute2f128_ps(va.vec, va.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
    tmp=_mm256_add_ps(tmp, va.vec);  // (1234+5678 5678+1234) = ( 1+5 2+6 3+7 4+8, 1+5 2+6 3+7 4+8)
    __m256 tmp2= _mm256_permute_ps(tmp, IMM_QUART4(1,0,3,2));             // 1234 5678 --> 3456 7812
    tmp=_mm256_add_ps(tmp,tmp2);       // (1+5)+(3+7), (2+6)+(4+8)...
    return sum[0];
  }

  // add
  static inline vrealF_t add(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(_mm256_add_ps(va.vec,vb.vec));
  }
  static inline void accum_add(vrealF_t &va,const vrealF_t &vb){
    va.vec=_mm256_add_ps(va.vec,vb.vec);
  }
  static inline vcomplexF_t add(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(_mm256_add_ps(va.vec,vb.vec));
  }
  static inline void accum_add(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec=_mm256_add_ps(va.vec,vb.vec);
  }

  // subtraction
  static inline vrealF_t sub(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(_mm256_sub_ps(va.vec,vb.vec));
  }
  inline void accum_sub(vrealF_t &va,const vrealF_t &vb){
    va.vec=_mm256_sub_ps(va.vec,vb.vec);
  }
  static inline vcomplexF_t sub(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(_mm256_sub_ps(va.vec,vb.vec));
  }
  inline void accum_sub(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec=_mm256_sub_ps(va.vec,vb.vec);
  }

  // reciprocal number
  static inline vrealF_t rcp(const vrealF_t &va){
    __m256 tmp=_mm256_set_ps(1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0);
    return vrealF_t(_mm256_div_ps(tmp, va.vec));
  }

  // mult
  static inline vrealF_t mult(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t( _mm256_mul_ps(va.vec, vb.vec));
  }

  // mult
  static inline vcomplexF_t mult(const vrealF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t( _mm256_mul_ps(va.vec, vb.vec));
  }

#ifdef AVX2
  static inline vcomplexF_t cmult(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    return vcomplexF_t(_mm256_fmaddsub_ps( a_real, vb.vec, a_imag ));
  }

  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    out.vec=_mm256_fmaddsub_ps( a_real, vb.vec, a_imag );
  }

  static inline void accum_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_real = _mm256_fmaddsub_ps( a_real, vb.vec, a_imag );
    out.vec = _mm256_add_ps(out.vec, a_real);
  }

  // out -= va*vb
  static inline void accum_sub_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_real = _mm256_fmaddsub_ps( a_real, vb.vec, a_imag );
    out.vec = _mm256_sub_ps(out.vec, a_real);
  }

#else
  static inline vcomplexF_t cmult(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1)) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br

    // (Ar Br, Ar Bi)-+(Ai Bi, Ai Br) =  (ArBr-AiBi , ArBi+AiBr)
    __m256 tmp = _mm256_mul_ps(a_real, vb.vec);
    return vcomplexF_t( _mm256_addsub_ps(tmp,a_imag) );
  }

  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br

    // (Ar Br, Ar Bi)-+(Ai Bi, Ai Br) =  (ArBr-AiBi , ArBi+AiBr)
    __m256 tmp = _mm256_mul_ps(a_real, vb.vec);
    out.vec= _mm256_addsub_ps(tmp,a_imag);
  }

  static inline void accum_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br

    // (Ar Br, Ar Bi)-+(Ai Bi, Ai Br) =  (ArBr-AiBi , ArBi+AiBr)
    a_real = _mm256_mul_ps(a_real, vb.vec);
    a_real = _mm256_addsub_ps(a_real,a_imag);
    out.vec= _mm256_add_ps(out.vec,a_real);
  }

  static inline void accum_sub_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br

    // (Ar Br, Ar Bi)-+(Ai Bi, Ai Br) =  (ArBr-AiBi , ArBi+AiBr)
    a_real = _mm256_mul_ps(a_real, vb.vec);
    a_real = _mm256_addsub_ps(a_real,a_imag);
    out.vec= _mm256_sub_ps(out.vec,a_real);
  }

#endif

  // mult:  conj(a)*b
#if defined (AVX2)
  static inline vcomplexF_t cmult_conj(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    return vcomplexF_t(_mm256_fmsubadd_ps( a_real, vb.vec, a_imag ));
  }

  static inline void set_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    out.vec=_mm256_fmsubadd_ps( a_real, vb.vec, a_imag );
  }

  static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_real=_mm256_fmsubadd_ps( a_real, vb.vec, a_imag );
    out.vec=_mm256_add_ps(out.vec, a_real);
  }

  static inline void accum_sub_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_real=_mm256_fmsubadd_ps( a_real, vb.vec, a_imag );
    out.vec=_mm256_sub_ps(out.vec, a_real);
  }

#else

  static inline vcomplexF_t cmult_conj(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_imag = _mm256_sub_ps( _mm256_setzero_ps(), a_imag);   // - AiBi, -AiBr
    a_real = _mm256_mul_ps( a_real, vb.vec);              // ArBr, ArBi
    return vcomplexF_t(_mm256_addsub_ps( a_real, a_imag )); // ArBr - (-AiBi), ArBi -AiBr
  }

static inline void set_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_imag = _mm256_sub_ps( _mm256_setzero_ps(), a_imag);   // - AiBi, -AiBr
    a_real = _mm256_mul_ps( a_real, vb.vec);              // ArBr, ArBi
    out.vec=_mm256_addsub_ps( a_real, a_imag ); // ArBr - (-AiBi), ArBi -AiBr
  }

static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_imag = _mm256_sub_ps( _mm256_setzero_ps(), a_imag);   // - AiBi, -AiBr
    a_real = _mm256_mul_ps( a_real, vb.vec);              // ArBr, ArBi
    a_real = _mm256_addsub_ps( a_real, a_imag ); // ArBr - (-AiBi), ArBi +(-AiBr)
    out.vec = _mm256_add_ps(out.vec, a_real);
  }

static inline void accum_sub_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m256 a_real = _mm256_moveldup_ps( va.vec ); // Ar Ar
    __m256 a_imag = _mm256_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm256_mul_ps( a_imag, _mm256_permute_ps( vb.vec, IMM_QUART4(2,3,0,1) ));  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    a_imag = _mm256_sub_ps( _mm256_setzero_ps(), a_imag);   // - AiBi, -AiBr
    a_real = _mm256_mul_ps( a_real, vb.vec);              // ArBr, ArBi
    a_real = _mm256_addsub_ps( a_real, a_imag ); // ArBr - (-AiBi), ArBi +(-AiBr)
    out.vec = _mm256_sub_ps(out.vec, a_real);
  }

#endif

  // sqrt(-1)
  static inline vcomplexF_t TimesI(const vcomplexF_t &in){
    __m256 tmp = _mm256_shuffle_ps(in.vec, in.vec, IMM_QUART4(2,3,0,1));  // re,im --> im, re
    return vcomplexF_t(_mm256_addsub_ps(_mm256_setzero_ps(), tmp));           // 0-im, 0+re
  }

  static inline void accum_TimesI(vcomplexF_t &out, const vcomplexF_t &in){
    __m256 tmp = _mm256_shuffle_ps(in.vec, in.vec, IMM_QUART4(2,3,0,1));  // re,im --> im, re
    out.vec = _mm256_addsub_ps(out.vec, tmp);                             // out_re-im, out_im+re
  }

  static inline vcomplexF_t TimesMinusI(const vcomplexF_t &in){
    __m256 tmp = _mm256_shuffle_ps(in.vec, in.vec, IMM_QUART4(2,3,0,1));  // re,im --> im, re
    tmp = _mm256_addsub_ps(_mm256_setzero_ps(), tmp);                         // -im, +re
    return vcomplexF_t(_mm256_sub_ps(_mm256_setzero_ps(), tmp));              // 0+im, 0-re
  }

  static inline void accum_TimesMinusI(vcomplexF_t &out, const vcomplexF_t &in){
    __m256 tmp = _mm256_shuffle_ps(in.vec, in.vec, IMM_QUART4(2,3,0,1));  // re,im --> im, re
    tmp = _mm256_addsub_ps(_mm256_setzero_ps(), tmp);                         // -im, +re
    out.vec = _mm256_sub_ps(out.vec, tmp);                                // out_re+im, out_im-re
  }

  static inline vcomplexF_t Times_pm(const vcomplexF_t &in){
    // (c1,c2,c3,c4) --> (c1, -c2, c3, -c4)
    __m256 tmp = _mm256_permute_ps(in.vec, IMM_QUART4(0,2,1,3));
    tmp = _mm256_addsub_ps(_mm256_setzero_ps(), tmp);
    return vcomplexF_t(_mm256_permute_ps(tmp, IMM_QUART4(0,2,1,3)));
  }

  static inline vcomplexF_t Times_mp(const vcomplexF_t &in){
    // (c1,c2,c3,c4) --> (-c1, c2, -c3, c4)
    __m256 tmp = _mm256_permute_ps(in.vec, IMM_QUART4(3,1,2,0));
    tmp = _mm256_addsub_ps(_mm256_setzero_ps(), tmp);
    return vcomplexF_t(_mm256_permute_ps(tmp, IMM_QUART4(3,1,2,0)));
  }


  // axpy: y := a*x+y
  static inline void axpy(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    accum_cmult(vy,va,vx);
  }

#ifdef AVX2
  static inline void axpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec=_mm256_fmadd_ps(va.vec, vx.vec, vy.vec);
  }
  static inline void axpy(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec=_mm256_fmadd_ps(va.vec, vx.vec, vy.vec);
  }
#else
  static inline void axpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    __m256 tmp = _mm256_mul_ps(va.vec, vx.vec);
    vy.vec=_mm256_add_ps(vy.vec, tmp);
  }
  static inline void axpy(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    __m256 tmp = _mm256_mul_ps(va.vec, vx.vec);
    vy.vec=_mm256_add_ps(vy.vec, tmp);
  }
#endif

  // alias
  static inline void vaxpy(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    axpy(vy,va,vx);
  }
  static inline void vaxpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    axpy(vy,va,vx);
  }
  static inline void vaxpy(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    axpy(vy,va,vx);
  }
  static inline void accum_abs2(vrealF_t &vy, const vcomplexF_t &vx){
    axpy(vy, vrealF_t(vx.vec), vrealF_t(vx.vec));
  }

  // y:= y-ax
#ifdef AVX2
  static inline void accum_sub_mult(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec=_mm256_fnmadd_ps(va.vec, vx.vec, vy.vec);
  }
  static inline void accum_sub_mult(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec=_mm256_fnmadd_ps(va.vec, vx.vec, vy.vec);
  }
#else
  static inline void accum_sub_mult(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    __m256 tmp = _mm256_mul_ps(va.vec, vx.vec);
    vy.vec=_mm256_sub_ps(vy.vec, tmp);
  }
  static inline void accum_sub_mult(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    __m256 tmp = _mm256_mul_ps(va.vec, vx.vec);
    vy.vec=_mm256_sub_ps(vy.vec, tmp);
  }
#endif


  // aypx: y := a*y+x
#ifdef AVX2
  static inline void aypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec = _mm256_fmadd_ps(va.vec, vy.vec, vx.vec);
  }
  static inline void aypx(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec = _mm256_fmadd_ps(va.vec, vy.vec, vx.vec);
  }
#else
  static inline void aypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec=_mm256_mul_ps(vy.vec, va.vec);
    vy.vec=_mm256_add_ps(vy.vec, vx.vec);
  }
  static inline void aypx(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec=_mm256_mul_ps(vy.vec, va.vec);
    vy.vec=_mm256_add_ps(vy.vec, vx.vec);
  }
#endif
  static inline void aypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    set_cmult(vy, vy, va);
    vy.vec=_mm256_add_ps(vy.vec, vx.vec);
  }
  // alias
  static inline void vaypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    aypx(vy, va, vx);
  }
  static inline void vaypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    aypx(vy, va, vx);
  }

  // scal
  static inline void scal(vcomplexF_t &vx, const vcomplexF_t &va){
    set_cmult(vx,va,vx);
  }
  static inline void scal(vrealF_t &vx, const vrealF_t &va){
    vx.vec = _mm256_mul_ps(vx.vec, va.vec);
  }
  static inline void scal(vcomplexF_t &vx, const vrealF_t &va){
    vx.vec = _mm256_mul_ps(vx.vec, va.vec);
  }

  // (re im) (re im)... --> (re re) (re re)
  static inline void complex2rere(vcomplexF_t &vx){
    vx.vec=_mm256_moveldup_ps(vx.vec); // AVX1,AVX2
    //    vx.vec=_mm256_shuffle_ps(vx.vec, vx.vec, IMM_QUART4(2,2,0,0)); // AVX1
    // other
    //    float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }
  }

  // (re im) (re im)... --> (im im) (im im)
  static inline void complex2imim(vcomplexF_t &vx){
    vx.vec=_mm256_movehdup_ps(vx.vec);  // AVX2,AVX1
    //    vx.vec.v=_mm256_shuffle_ps(vx.vec.v, vx.vec.v, IMM_QUART4(3,3,1,1)); // AVX1
    // other
    //    float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i]=tmp[2*i+1];
    //    }
  }

  // (re im) (re im)... --> (re 0) (re 0)
  static inline void complex2re(vcomplexF_t &vx){
    __m256 tmp=_mm256_set_ps(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    vx.vec=_mm256_mul_ps(tmp,vx.vec);
  }

  // (re im) (re im)... --> (0 im) (0 im)
  static inline void complex2im(vcomplexF_t &vx){
    __m256 tmp=_mm256_set_ps(1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0);
    vx.vec=_mm256_mul_ps(tmp,vx.vec);
  }


  // TODO
  // (re1 im1),... (re1' im1')... --> (re1 re1') (re2 re2')...
  static inline void rere2complex(vcomplexF_t &vx1, const vcomplexF_t &vx2){
    // other
    //    float *tmp1=(float*)&vx1;
    //    float *tmp2=(float*)&vx2;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp1[2*i+1]=tmp2[2*i];
    //    }
    __m256i tmp=_mm256_shuffle_epi32((__m256i)vx2.vec, IMM_QUART4(2,3,0,1));
    vx1.vec=(__m256)_mm256_blend_epi32((__m256i)vx1.vec, tmp, IMM_BIN8(1,0,1,0,1,0,1,0));
  }

  // copy as vector varibles
  static inline void copy(vcomplexF_t &va, const vrealF_t &vb){
    va.vec=vb.vec;
  }

  static inline void copy(vrealF_t &va, const vcomplexF_t &vb){
    va.vec=vb.vec;
  }

} // anonymous

#endif
