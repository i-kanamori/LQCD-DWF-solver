/*!

        @file    vector_arith_512-inc.h

        @brief   arithmetics of simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
// some arithmetics : using Grid
#ifndef VECTOT_ARITH_GRID_H
#define VECTOT_ARITH_GRID_H


/* Immediate value, as Quaternary (4 digits) */
#define IMM_QUART4(A,B,C,D) ((A<<6)|(B<<4)|(C<<2)|(D))
/* Immediate value, as Binary (8 digits) */
#define IMM_BIN8(A,B,C,D,E,F,G,H) ((A<<7)|(B<<6)|(C<<5)|(D<<4)|(E<<3)|(F<<2)|(G<<1)|H)

  /*
    implemetations in this file are largely affected by
    Grid: Data parallel C++ mathematical object library
    https://github.com/paboyle/Grid
  */


namespace {

  ///////////////////////////////////////
  // double
  //////////////////////////////////////

  // broadcast
  static inline void bcast(vrealD_t &va, const realD_t a){
    va.vec=_mm512_set1_pd(a);
  }
#ifdef KNC
  static inline void bcast(vcomplexD_t &va, const complexD_t a){
    double re=real(a);
    double im=imag(a);
    va.vec=_mm512_set_pd(im,re,im,re,im,re,im,re);
  }
#else
  static inline void bcast(vcomplexD_t &va, const complexD_t a){
    double re=real(a);
    double im=imag(a);
    va.vec=_mm512_set4_pd(im,re,im,re);
    //   __m128d tmp=_mm_set_pd(im,re);
    //   va.vec=_mm512_broadcast_f64x2(tmp); // requires  AVX512DQ
  }
#endif

  // streaming store
  static inline void sstore(vcomplexD_t *vbuf, const vcomplexD_t &va){
    _mm512_store_pd((double*) &(vbuf->vec), va.vec);
  }

  static inline void sstore(vrealD_t *vbuf, const vrealD_t &va){
    _mm512_store_pd((double*) &(vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealD_t &va){
    va.vec=_mm512_setzero_pd();
  }
  static inline void set_zero(vcomplexD_t &va){
    va.vec=_mm512_setzero_pd();
  }

  //complex conjuate
  static inline vcomplexD_t conj(const vcomplexD_t &in){
    return vcomplexD_t(_mm512_mask_sub_pd(in.vec, 0xAA, _mm512_setzero_pd(), in.vec));
  }

  // reduction
  static inline realD_t reduce(const vrealD_t &va){
    return _mm512_reduce_add_pd(va.vec);
  }
  static inline complexD_t reduce(const vcomplexD_t &va){
    return complexD_t( _mm512_mask_reduce_add_pd(0x55, va.vec),
                       _mm512_mask_reduce_add_pd(0xAA, va.vec));
  }

  // add
  static inline vrealD_t add(const vrealD_t &va, const vrealD_t &vb){
    return vcomplexD_t(_mm512_add_pd(va.vec, vb.vec));
  }
  static inline void accum_add(vrealD_t &va,const vrealD_t &vb){
    va.vec=_mm512_add_pd(va.vec,vb.vec);
  }
  static inline vcomplexD_t add(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(_mm512_add_pd(va.vec, vb.vec));
  }
  static inline void accum_add(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec=_mm512_add_pd(va.vec,vb.vec);
  }

  // subtraction
  static inline vrealD_t sub(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(_mm512_sub_pd(va.vec,vb.vec));
  }
  static inline void accum_sub(vrealD_t &va,const vrealD_t &vb){
    va.vec=_mm512_sub_pd(va.vec, vb.vec);
  }
  static inline vcomplexD_t sub(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(_mm512_sub_pd(va.vec,vb.vec));
  }
  static inline void accum_sub(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec=_mm512_sub_pd(va.vec, vb.vec);
  }

  // reciprocal number
  static inline vrealD_t rcp(const vrealD_t &va){
    __m512d tmp=_mm512_set1_pd(1.0);
    return vrealD_t(_mm512_div_pd(tmp, va.vec));
  }

  // mult
  static inline vrealD_t mult(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t( _mm512_mul_pd(va.vec, vb.vec));
  }

  // mult
  static inline vcomplexD_t mult(const vrealD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t( _mm512_mul_pd(va.vec, vb.vec));
  }

#ifdef KNC
  static inline vcomplexD_t cmult(const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) -->tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
      __m512d tmp0  = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB ); 
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_real=(r0 r0) (r1 r1) (r2 r2) (r3 r3) :  0xAA =   (F T) (F T) (F T) (F T)
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) --> a_imag = (-i0 i0) (-i1 i1) (-i2 i2) (-i3 i3) : 0x55 = (T F) (T F) (T F) (T F)
      __m512d a_imag = _mm512_mask_sub_pd(va.vec, 0x55, vzero, tmp0);
      // tmp1=(Ar Br, Ar Bi),...
      __m512d tmp1   = _mm512_mul_pd(a_real, vb.vec);
      // vb=(r0 i0)...  --> tmp0=(i0 r0)...
      tmp0 = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      // (Ar Br, Ar Bi) + (Bi,Br)*(-Ai, Ai) = (Ar Br - Bi Ai, Ar Bi + Br Ai)
      return  vcomplexD_t(_mm512_fmadd_pd(tmp0, a_imag, tmp1));
  }
  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      __m512d tmp0  = _mm512_swizzle_pd(va.vec, _MM_SWIZ_REG_CDAB);
      __m512d a_real = (__m512d)_mm512_mask_or_epi64((__m512i)va.vec, 0xAA, (__m512i)vzero, (__m512i)tmp0);
      __m512d a_imag = _mm512_mask_sub_pd(va.vec, 0x55, vzero, tmp0);
      out.vec = _mm512_mul_pd( a_real, vb.vec );
      tmp0      = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }

  static inline void accum_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) -->tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
      __m512d tmp0  = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB );
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_real=(r0 r0) (r1 r1) (r2 r2) (r3 r3) :  0xAA =   (F T) (F T) (F T) (F T)
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) --> a_imag = (-i0 i0) (-i1 i1) (-i2 i2) (-i3 i3) : 0x55 = (T F) (T F) (T F) (T F)
      __m512d a_imag = _mm512_mask_sub_pd(va.vec, 0x55, vzero, tmp0);
      // vb=(r0 i0)...  --> tmp0=(i0 r0)...
      tmp0 = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      // out=(re+Ar Br, im+Ar Bi),...
      out.vec = _mm512_fmadd_pd(a_real, vb.vec, out.vec);
      // (re+Ar Br, im+Ar Bi) + (Bi,Br)*(-Ai, Ai) = (re+Ar Br - Bi Ai, im+Ar Bi + Br Ai)
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }

  static inline void accum_sub_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) -->tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
      __m512d tmp0  = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB ); 
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_real=(r0 r0) (r1 r1) (r2 r2) (r3 r3) :  0xAA =   (F T) (F T) (F T) (F T)
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      //   va=(r0 i0) (r1 i1) (r2 i2) (r3 i3)
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
      //   --> a_imag = (-i0 i0) (-i1 i1) (-i2 i2) (-i3 i3) : 0xAA =   (F T) (F T) (F T) (F T)
      __m512d a_imag = _mm512_mask_sub_pd(tmp0, 0xAA, vzero, va.vec);
      // vb=(r0 i0)...  --> tmp0=(i0 r0)...
      tmp0 = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      // out=(re -Ar Br, im-Ar Bi),...
      out.vec = _mm512_fnmadd_pd(a_real, vb.vec, out.vec);
      // (re-Ar Br, im-Ar Bi) + (Bi,Br)*(Ai, -Ai) = (re-Ar Br + Bi Ai, im-Ar Bi - Br Ai)
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }

#else
  static inline vcomplexD_t cmult(const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ) );
    return vcomplexD_t(_mm512_fmaddsub_pd( a_real, vb.vec, a_imag ));
  }
  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    out.vec = _mm512_mul_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ) );
    out.vec =_mm512_fmaddsub_pd( a_real, vb.vec, out.vec );
  }
  static inline void accum_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    // (Ai Ai)*(Bi Br) -+ (re im) = (-re + Ai Bi, -im+Ai Br)
    out.vec = _mm512_fmaddsub_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ), out.vec );
    // (Ar Ar)*(Br Bi) -+ (-re + Ai Bi, -im+Ai Br)
    // = (re + Ar Br -Ai Bi, im + Ar Bi + Ai Br)
    out.vec = _mm512_fmaddsub_pd( a_real, vb.vec, out.vec );
  }

  static inline void accum_sub_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    a_real = _mm512_sub_pd(_mm512_setzero_pd(), a_real);  // (-Ar -Ar)
    // (Ai Ai)*(Bi Br)+-(re im) = (re + Ai Bi, -im + Ai Br)
    out.vec = _mm512_fmsubadd_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ), out.vec );
    // (-Ar -Ar) * (Br Bi) +- (re + Ai Bi, -im + Ai Br)
    // = ( re - Ar Br + Ai Bi, im -Ar Bi - Ai Br)
    out.vec = _mm512_fmsubadd_pd( a_real, vb.vec, out.vec );
  }
#endif

  // mult:  conj(a)*b
#ifdef KNC
  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) -->tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
      __m512d tmp0  = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB );
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_real=(r0 r0) (r1 r1) (r2 r2) (r3 r3) :  0xAA =  (F T) (F T) (F T) (F T)
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_imag = (i0 -i0) (i1 -i1) (i2 -i2) (i3 -i3) : 0xAA =  (F T) (F T) (F T) (F T)
      __m512d a_imag = _mm512_mask_sub_pd(tmp0, 0xAA, vzero, va.vec);
      // vb=(r0 i0)...  --> tmp0=(i0 r0)...
      tmp0 = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      // out=(Ar Br, Ar Bi),...
      out.vec  = _mm512_mul_pd(a_real, vb.vec);
      // (Ar Br, Ar Bi) + (Bi,Br)*(Ai, -Ai) = (Ar Br + Bi Ai, Ar Bi - Br Ai)
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }
  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      __m512d tmp0  = _mm512_swizzle_pd(va.vec, _MM_SWIZ_REG_CDAB);
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      __m512d a_imag = _mm512_mask_sub_pd(tmp0, 0xAA, vzero, va.vec);
      tmp0      = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      out.vec = _mm512_fmadd_pd( a_real, vb.vec, out.vec );
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }

  static inline void accum_sub_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
      __m512d vzero = _mm512_setzero_pd();
      __m512d tmp0  = _mm512_swizzle_pd(va.vec, _MM_SWIZ_REG_CDAB);
      __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
      __m512d a_imag = _mm512_mask_sub_pd(va.vec, 0x55, vzero, tmp0);
      tmp0      = _mm512_swizzle_pd( vb.vec, _MM_SWIZ_REG_CDAB );
      out.vec = _mm512_fnmadd_pd( a_real, vb.vec, out.vec );
      out.vec = _mm512_fmadd_pd( tmp0, a_imag, out.vec );
  }
#else
  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    // slight modification of a*b in Grid
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ) );
    out.vec =_mm512_fmsubadd_pd( a_real, vb.vec, a_imag );
  }
  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    a_imag = _mm512_fmsubadd_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ) , out.vec);
    out.vec =_mm512_fmsubadd_pd( a_real, vb.vec, a_imag );
  }
  static inline void accum_sub_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    // (Ai Ai) * (Bi Br) -+ (re im) = ( -re + Ai Bi, im + Ai Br)
    a_imag = _mm512_fmaddsub_pd( a_imag, _mm512_permute_pd( vb.vec, 0x55 ) , out.vec);
    a_real = _mm512_sub_pd(_mm512_setzero_pd(), a_real); // (-Ar -Ar)
    // (-Ar -Ar) * (Br Bi) -+ (- re + Ai Bi, im + Ai Br)
    // = (re - Ar Br - Ai Bi, im - Ar Bi  Ai Br)
    out.vec =_mm512_fmaddsub_pd( a_real, vb.vec, a_imag );
  }
#endif
  static inline vcomplexD_t cmult_conj(const vcomplexD_t &va, const vcomplexD_t &vb){
    vcomplexD_t out;
    set_cmult_conj(out,va,vb);
    return out;
  }

  // sqrt(-1)
#ifdef KNC
  static inline vcomplexD_t TimesI(const vcomplexD_t &in){
    __m512d tmp = _mm512_swizzle_pd(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    return vcomplexD_t( _mm512_mask_sub_pd(tmp, 0x55, _mm512_setzero_pd(), tmp) ); // (-im, re)
  }
  static inline vcomplexD_t TimesMinusI(const vcomplexD_t &in){
    __m512d tmp = _mm512_swizzle_pd(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    return vcomplexD_t( _mm512_mask_sub_pd(tmp, 0xAA, _mm512_setzero_pd(), tmp) ); // (im, -re)
  }
  static inline void accum_TimesI(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp = _mm512_swizzle_pd(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    out.vec = _mm512_mask_add_pd(out.vec, 0xAA, out.vec,tmp); // (out_re,    out_im+re)
    out.vec = _mm512_mask_sub_pd(out.vec, 0x55, out.vec,tmp); // (out_re-im, out_im+re)
  }

  static inline void accum_TimesI_compare(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp = _mm512_swizzle_pd(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    tmp=_mm512_mask_sub_pd(tmp, 0x55, _mm512_setzero_pd(), tmp);
    out.vec = _mm512_add_pd(out.vec, tmp);
  }

  static inline void accum_TimesMinusI(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp = _mm512_swizzle_pd(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    out.vec = _mm512_mask_sub_pd(out.vec, 0xAA, out.vec,tmp); // (out_re,    out_im-re)
    out.vec = _mm512_mask_add_pd(out.vec, 0x55, out.vec,tmp); // (out_re+im, out_im-re)
  }
#else // AVX512
  static inline vcomplexD_t TimesI(const vcomplexD_t &in){
    __m512d tmp = _mm512_shuffle_pd(in.vec, in.vec, 0x55); // (im,re)
    return vcomplexD_t( _mm512_mask_sub_pd(tmp, 0x55, _mm512_setzero_pd(), tmp) ); // (-im, re)
  }
  static inline vcomplexD_t TimesMinusI(const vcomplexD_t &in){
    __m512d tmp = _mm512_shuffle_pd(in.vec, in.vec, 0x55); // (im,re)
    return vcomplexD_t( _mm512_mask_sub_pd(tmp, 0xAA, _mm512_setzero_pd(), tmp) ); // (im, -re)
  }
  static void inline accum_TimesI(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp=_mm512_shuffle_pd(in.vec, in.vec, 0x55);
    tmp=_mm512_mask_sub_pd(tmp, 0x55, _mm512_setzero_pd(), tmp);
    out.vec=_mm512_add_pd(out.vec, tmp);
  }
  static void inline accum_TimesI_compare(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp=_mm512_shuffle_pd(in.vec,in.vec,0x55);    // (im, re)
    out.vec = _mm512_mask_add_pd(out.vec, 0xAA, out.vec, tmp); // (out_re,    out_im+re)
    out.vec = _mm512_mask_sub_pd(out.vec, 0x55, out.vec, tmp); // (out_re-im, out_im+re)
  }

  static void inline accum_TimesMinusI(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp=_mm512_shuffle_pd(in.vec,in.vec,0x55);
    tmp = _mm512_mask_sub_pd(tmp, 0xAA, _mm512_setzero_pd(), tmp);
    out.vec=_mm512_add_pd(out.vec, tmp);
  }
  static void inline accum_TimesMinusI_compare(vcomplexD_t &out, const vcomplexD_t &in){
    __m512d tmp=_mm512_shuffle_pd(in.vec,in.vec,0x55);
    out.vec = _mm512_mask_sub_pd(out.vec, 0xAA, out.vec,tmp); // (out_re,    out_im-re)
    out.vec = _mm512_mask_add_pd(out.vec, 0x55, out.vec,tmp); // (out_re+im, out_im-re)
  }

#endif

  static inline vcomplexD_t Times_pm(const vcomplexD_t &in){
    // (-- ++ -- ++)
    //  11 00 11 00   = 0xCC
    __m512d tmp=_mm512_mask_sub_pd(in.vec, 0xCC, _mm512_setzero_pd(), in.vec);
    return vcomplexD_t(tmp);
  }

  static inline vcomplexD_t Times_mp(const vcomplexD_t &in){
    // (++ -- ++ --)
    //  00 11 00 11  = 0x33
    __m512d tmp=_mm512_mask_sub_pd(in.vec, 0x33, _mm512_setzero_pd(), in.vec);
    return vcomplexD_t(tmp);
  }

  // axpy: y := ax + y
  static inline void axpy(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    accum_cmult(vy,va,vx);
  }

  static inline void axpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec = _mm512_fmadd_pd( va.vec, vx.vec, vy.vec );
  }

  static inline void axpy(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec = _mm512_fmadd_pd( va.vec, vx.vec, vy.vec );
  }

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
  static inline void accum_sub_mult(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec = _mm512_fnmadd_pd( va.vec, vx.vec, vy.vec );
  }
  static inline void accum_sub_mult(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec = _mm512_fnmadd_pd( va.vec, vx.vec, vy.vec );
  }

  // aypx: y := a*y+x
#ifdef KNC
  static inline void aypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    //    vy.vec=va.vec*vy.vec + vx.vec;
    __m512d vzero = _mm512_setzero_pd();
    // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) -->tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)
    __m512d tmp0  = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB );
    // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3) --> a_real=(r0 r0) (r1 r1) (r2 r2) (r3 r3) :  0xAA =   (F T) (F T) (F T) (F T)
    __m512d a_real = (__m512d)_mm512_mask_or_epi64( (__m512i)va.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
    // va=(r0 i0) (r1 i1) (r2 i2) (r3 i3) --> a_imag = (-i0 i0) (-i1 i1) (-i2 i2) (-i3 i3) : 0x55 = (T F) (T F) (T F) (T F)
    __m512d a_imag = _mm512_mask_sub_pd(va.vec, 0x55, vzero, tmp0);
    // vb=(r0 i0)...  --> tmp0=(i0 r0)...
    __m512d tmp1 = _mm512_swizzle_pd( vy.vec, _MM_SWIZ_REG_CDAB );
    // out=(Xr+Ar Yr, Xi+Ar Yi),...
    __m512d tmp2 = _mm512_fmadd_pd(a_real, vy.vec, vx.vec);
    // (Xr+Ar Br, Xi+Ar Bi) + (Yi,Yr)*(-Ai, Ai) = (Xr+Ar Yr - Yi Ai, Xi+Ar Yi + Yr Ai)
    vy.vec = _mm512_fmadd_pd( tmp1, a_imag, tmp2 );
  }
#else

  static inline void aypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    __m512d y_real = _mm512_shuffle_pd( vy.vec, vy.vec, 0x00 );
    __m512d y_imag = _mm512_shuffle_pd( vy.vec, vy.vec, 0xFF );

    vy.vec = _mm512_fmaddsub_pd( y_imag, _mm512_permute_pd( va.vec, 0x55 ), vx.vec );
    vy.vec = _mm512_fmaddsub_pd( y_real, va.vec, vy.vec );
  }
#endif

  static inline void aypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec = _mm512_fmadd_pd( va.vec, vy.vec, vx.vec );
  }
  static inline void aypx(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    vy.vec = _mm512_fmadd_pd( va.vec, vy.vec, vx.vec );
  }

  // alias
  static inline void vaypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    aypx(vy, va, vx);
  }
  static inline void vaypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    aypx(vy, va, vx);
  }

  // scal
#ifdef KNC
  static inline void scal(vcomplexD_t &vx, const vcomplexD_t &va){
    __m512d vzero = _mm512_setzero_pd();
    // vx=(Xr, Xi) --> tmp0=(Xi, Xr)
    __m512d tmp0  = _mm512_swizzle_pd( vx.vec, _MM_SWIZ_REG_CDAB );
    //  x_real=(Xr Xr) :  0xAA =  (F T) (F T) (F T) (F T)
    __m512d x_real = (__m512d)_mm512_mask_or_epi64( (__m512i)vx.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );
    // x_imag = (-Xi, Xi)
    __m512d x_imag = _mm512_mask_sub_pd(vx.vec, 0x55, vzero, tmp0);
    //tmp0=(Ai Ar)
    tmp0 = _mm512_swizzle_pd( va.vec, _MM_SWIZ_REG_CDAB );
    // vx =(Ar Xr, Ai Xr),...
    vx.vec  = _mm512_mul_pd(x_real, va.vec);
    // (Ar Xr, Ai Xr) + (-Xi,Xi)*(Ai, Ar) = (Ar Xr - Xi Ai, Ar Xi + Xr Ai)
    vx.vec = _mm512_fmadd_pd( tmp0, x_imag, vx.vec );
  }
#else // AVX512
  static inline void scal(vcomplexD_t &vx, const vcomplexD_t &va){
    __m512d a_real = _mm512_shuffle_pd( va.vec, va.vec, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec, va.vec, 0xFF );
    a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( vx.vec, 0x55 ) );
    vx.vec =_mm512_fmaddsub_pd( a_real, vx.vec, a_imag );
  }
#endif

  static inline void scal(vrealD_t &vx, const vrealD_t &va){
    vx.vec = _mm512_mul_pd(vx.vec, va.vec);
  }
  static inline void scal(vcomplexD_t &vx, const vrealD_t &va){
    vx.vec = _mm512_mul_pd(vx.vec, va.vec);
  }

  // (re im) (re im)... --> (re re) (re re)...
#ifdef KNC
  static inline void complex2rere(vcomplexD_t &vx){
    __m512d vzero = _mm512_setzero_pd();
    __m512d tmp0  = _mm512_swizzle_pd( vx.vec, _MM_SWIZ_REG_CDAB );
    vx.vec = (__m512d)_mm512_mask_or_epi64( (__m512i)vx.vec, 0xAA,(__m512i)vzero, (__m512i)tmp0 );

    //    double *tmp=(double*)&vx;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }

  }
#else // AVX512
  static inline void complex2rere(vcomplexD_t &vx){
    vx.vec=_mm512_shuffle_pd( vx.vec, vx.vec, 0x00 );
  }
#endif

  // (re im) (re im)... --> (im im) (im im)
#ifdef KNC
  static inline void complex2imim(vcomplexD_t &vx){
    __m512d vzero = _mm512_setzero_pd();
    __m512d tmp0  = _mm512_swizzle_pd( vx.vec, _MM_SWIZ_REG_CDAB );
    vx.vec = (__m512d)_mm512_mask_or_epi64( (__m512i)tmp0, 0xAA,(__m512i)vzero, (__m512i)vx.vec );
  //    double *tmp=(double*)&vx;
  //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
  //      tmp[2*i]=tmp[2*i+1];
  }
#else // AVX512
  static inline void complex2imim(vcomplexD_t &vx){
    vx.vec=_mm512_shuffle_pd( vx.vec, vx.vec, 0xFF );
  }
#endif

  // (re im) (re im)... --> (re 0 ) (re 0)
  static inline void complex2re(vcomplexD_t &vx){
    vx.vec=_mm512_mask_mov_pd(vx.vec, 0xAA, _mm512_setzero_pd());
  }

  // (re im) (re im)... --> (0 im) (0 im)
  static inline void complex2im(vcomplexD_t &vx){
    vx.vec=_mm512_mask_mov_pd(vx.vec, 0x55, _mm512_setzero_pd());
  }

  // TODO
  // (re1 im1),... (re1' im1')... --> (re1 re1') (re2 re2')...
  static inline void rere2complex(vcomplexD_t &vx1, const vcomplexD_t &vx2){
    vcomplexD_t tmp=vx2;
    complex2rere(tmp);
    vx1.vec=_mm512_mask_mov_pd(tmp.vec, 0x55, vx1.vec);

    // other
    //double *tmp1=(double*)&vx1;
    //double *tmp2=(double*)&vx2;
    //for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //  tmp1[2*i+1]=tmp2[2*i];
    //}
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
  static inline void bcast(vrealF_t &va, const realF_t a){
    va.vec=_mm512_set1_ps(a);
  }
#ifdef KNC
  static inline void bcast(vcomplexF_t &va, const complexF_t a){
    float re=real(a);
    float im=imag(a);
    va.vec=_mm512_set_ps(im,re,im,re,im,re,im,re,im,re,im,re,im,re,im,re);
  }
#else
  static inline void bcast(vcomplexF_t &va, const complexF_t a){
    float re=real(a);
    float im=imag(a);
    va.vec=_mm512_set4_ps(im,re,im,re);
    //    __m128 tmp=_mm_set_ps(im,re,im,re);
    //    va.vec=_mm512_broadcast_f32x4(tmp);
  }
#endif

  // streaming store
  static inline void sstore(vcomplexF_t *vbuf, const vcomplexF_t &va){
    _mm512_store_ps((float*) &(vbuf->vec), va.vec);
  }

  static inline void sstore(vrealF_t *vbuf, const vrealF_t &va){
    _mm512_store_ps((float*) &(vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealF_t &va){
    va.vec=_mm512_setzero_ps();
  }
  static inline void set_zero(vcomplexF_t &va){
    va.vec=_mm512_setzero_ps();
  }
  //complex conjuate
  static inline vcomplexF_t conj(const vcomplexF_t &in){
    return vcomplexF_t(_mm512_mask_sub_ps(in.vec, 0xAAAA, _mm512_setzero_ps(), in.vec));
  }

  // reduction
  static inline realF_t reduce(const vrealF_t &va){
    return _mm512_reduce_add_ps(va.vec);

  }
  static inline complexF_t reduce(const vcomplexF_t &va){
    return complexF_t( _mm512_mask_reduce_add_ps(0x5555, va.vec),
                       _mm512_mask_reduce_add_ps(0xAAAA, va.vec));
  }

  // add
  static inline vrealF_t add(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(_mm512_add_ps(va.vec, vb.vec));
  }
  static inline void accum_add(vrealF_t &va,const vrealF_t &vb){
    va.vec=_mm512_add_ps(va.vec, vb.vec);
  }
  static inline vcomplexF_t add(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(_mm512_add_ps(va.vec, vb.vec));
  }
  static inline void accum_add(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec=_mm512_add_ps(va.vec, vb.vec);
  }

  // subtraction
  static inline vrealF_t sub(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(_mm512_sub_ps(va.vec, vb.vec));
  }
  inline void accum_sub(vrealF_t &va,const vrealF_t &vb){
    va.vec=_mm512_sub_ps(va.vec, vb.vec);
  }
  static inline vcomplexF_t sub(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(_mm512_sub_ps(va.vec, vb.vec));
  }
  inline void accum_sub(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec=_mm512_sub_ps(va.vec, vb.vec);
  }

  // reciprocal number
  static inline vrealF_t rcp(const vrealF_t &va){
    __m512 tmp=_mm512_set1_ps(1.0);
    return vrealF_t(_mm512_div_ps(tmp, va.vec));
  }

  // mult
  static inline vrealF_t mult(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(_mm512_mul_ps(va.vec, vb.vec) );
  }

  // mult
  static inline vcomplexF_t mult(const vrealF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(_mm512_mul_ps(va.vec, vb.vec) );
  }

  // complex multplication
#ifdef KNC
  static inline vcomplexF_t cmult(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Ai) --> tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // tmp0=(Ai Ar) --> a_real=(Ar Ar) :  0xAAAA =   (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // va=(Ar Ai) --> a_imag = (-Ai Ai) : 0x5555 = (T F) (T F) (T F) (T F)....
    __m512 a_imag = _mm512_mask_sub_ps(va.vec, 0x5555, vzero, tmp0);
    // tmp1=(Ar Br, Ar Bi),...
    __m512 tmp1   = _mm512_mul_ps(a_real, vb.vec);
    // vb=(r0 i0)...  --> tmp0=(i0 r0)...
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // (Ar Br, Ar Bi) + (Bi,Br)*(-Ai, Ai) = (Ar Br - Bi Ai, Ar Bi + Br Ai)
    return  vcomplexF_t(_mm512_fmadd_ps(tmp0, a_imag, tmp1));

  }

  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Ai) --> tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // tmp0=(Ai Ar) --> a_real=(Ar Ar) :  0xAAAA =   (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // va=(Ar Ai) --> a_imag = (-Ai Ai) : 0x5555 = (T F) (T F) (T F) (T F)....
    __m512 a_imag = _mm512_mask_sub_ps(va.vec, 0x5555, vzero, tmp0);
    // out =(Ar Br, Ar Bi),...
    out.vec  = _mm512_mul_ps(a_real, vb.vec);
    // vb=(r0 i0)...  --> tmp0=(i0 r0)...
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // (Ar Br, Ar Bi) + (Bi,Br)*(-Ai, Ai) = (Ar Br - Bi Ai, Ar Bi + Br Ai)
    out.vec = _mm512_fmadd_ps(tmp0, a_imag, out.vec);
  }
  static inline void accum_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Ai) -->tmp0=(Ar Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // tmp0=(Ai Ar) --> a_real=(Ar Ar) :  0xAAAA =   (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // va=(Ar Ai) --> a_imag = (-Ai Ai) : 0x5555 = (T F) (T F) (T F) (T F)....
    __m512 a_imag = _mm512_mask_sub_ps(va.vec, 0x5555, vzero, tmp0);
    // vb=(Ar Ai) --> tmp0=(Bi Br)...
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // out=(re+Ar Br, im+Ar Bi),...
    out.vec = _mm512_fmadd_ps(a_real, vb.vec, out.vec);
    // (re+Ar Br, im+Ar Bi) + (Bi,Br)*(-Ai, Ai) = (re+Ar Br - Bi Ai, im+Ar Bi + Br Ai)
    out.vec = _mm512_fmadd_ps( tmp0, a_imag, out.vec );
  }

  static inline void accum_sub_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Ai) -->tmp0=(Ar Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // tmp0=(Ai Ar) --> a_real=(Ar Ar) :  0xAAAA =   (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    //   va=(r0 i0) (r1 i1) (r2 i2) (r3 i3)...
    // tmp0=(i0 r0) (i1 r1) (i2 r2) (i3 r3)...
    //   --> a_imag = (-i0 i0) (-i1 i1) (-i2 i2) (-i3 i3)... : 0xAAAA =   (F T) (F T) (F T) (F T)...
    __m512 a_imag = _mm512_mask_sub_ps(tmp0, 0xAAAA, vzero, va.vec);
    // vb=(Ar Ai) --> tmp0=(Bi Br)...
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // out=(re-Ar Br, im-Ar Bi),...
    out.vec = _mm512_fnmadd_ps(a_real, vb.vec, out.vec);
    // (re-Ar Br, im-Ar Bi) + (Bi,Br)*(Ai, -Ai) = (re-Ar Br + Bi Ai, im-Ar Bi - Br Ai)
    out.vec = _mm512_fmadd_ps( tmp0, a_imag, out.vec );
  }

#else // AVX512

  static inline vcomplexF_t cmult(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // (Ar, Ar)
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // (Ai, Ai)
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ) ); // (Ai,Ai) (Bi,Br)
    return vcomplexF_t( _mm512_fmaddsub_ps( a_real, vb.vec, a_imag ) );

  }

  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // (Ar, Ar)
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // (Ai, Ai)
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ) );  // (Ai,Ai) (Bi,Br)
    out.vec =_mm512_fmaddsub_ps( a_real, vb.vec, a_imag );

  }

  static inline void accum_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // (Ai, Ai)
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // (Ar, Ar)
    a_imag = _mm512_fmaddsub_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ), out.vec );  // (Ai,Ai) (Bi,Br) + (re, im) = (-re + Ai Bi, im+AiBr)
    out.vec =_mm512_fmaddsub_ps( a_real, vb.vec, a_imag ); // (Ar Br, Ar Bi)+(re-Ai Bi, im+Ai Br)
  }

  static inline void accum_sub_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // (Ai, Ai)
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // (Ar, Ar)
    a_imag = _mm512_fmsubadd_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ), out.vec );  // (Ai,Ai) (Bi,Br) +- (re, im) = (re + Ai Bi, -im+AiBr)
    a_real = _mm512_sub_ps(_mm512_setzero_ps(), a_real); // (-Ar,-Ar)
    out.vec =_mm512_fmsubadd_ps( a_real, vb.vec, a_imag ); // (-Ar Br, -Ar Bi)+-(re+Ai Bi, -im+Ai Br)
  }
#endif


  // complex conjugate multplication
  // mult:  conj(a)*b
#ifdef KNC

  static inline vcomplexF_t cmult_conj(const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Bi) -->tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    //  a_real=(Ar Ar) :  0xAAAA = (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // a_imag = (Ai, -Ai) : 0xAA =  (F T) (F T) (F T) (F T)....
    __m512 a_imag = _mm512_mask_sub_ps(tmp0, 0xAAAA, vzero, va.vec);
    // tmp0=(Bi Br)
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // a_real=(Ar Br, Ar Bi),...
    a_real  = _mm512_mul_ps(a_real, vb.vec);
    // (Ar Br, Ar Bi) + (Bi,Br)*(Ai, -Ai) = (Ar Br + Bi Ai, Ar Bi - Br Ai)
    return vcomplexF_t( _mm512_fmadd_ps( tmp0, a_imag, a_real ) );
  }

  static inline void set_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Bi) -->tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    //  a_real=(Ar Ar) :  0xAAAA = (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // a_imag = (Ai, -Ai) : 0xAA =  (F T) (F T) (F T) (F T)....
    __m512 a_imag = _mm512_mask_sub_ps(tmp0, 0xAAAA, vzero, va.vec);
    // tmp0=(Bi Br)
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // out=( Ar Br,  Ar Bi),...
    out.vec  = _mm512_mul_ps(a_real, vb.vec);
    // (Ar Br, Ar Bi) + (Bi,Br)*(Ai, -Ai) = (Ar Br + Bi Ai, Ar Bi - Br Ai)
    out.vec = _mm512_fmadd_ps( tmp0, a_imag, out.vec );
  }

  static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Bi) -->tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    //  a_real=(Ar Ar) :  0xAAAA = (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // a_imag = (Ai, -Ai) : 0xAA =  (F T) (F T) (F T) (F T)....
    __m512 a_imag = _mm512_mask_sub_ps(tmp0, 0xAAAA, vzero, va.vec);
    // tmp0=(Bi Br)
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // out=(re + Ar Br, im + Ar Bi),...
    out.vec  = _mm512_fmadd_ps(a_real, vb.vec, out.vec);
    // (re + Ar Br, im + Ar Bi) + (Bi,Br)*(Ai, -Ai) = (re + Ar Br + Bi Ai, im + Ar Bi - Br Ai)
    out.vec = _mm512_fmadd_ps( tmp0, a_imag, out.vec );
  }

  static inline void accum_sub_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Bi) -->tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    //  a_real=(Ar Ar) :  0xAAAA = (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // a_imag = (-Ai, Ai) : 0x5555 =  (T F) (T F) (T F)....
    __m512 a_imag = _mm512_mask_sub_ps(va.vec, 0x5555, vzero, tmp0);
    // tmp0=(Bi Br)
    tmp0 = _mm512_swizzle_ps( vb.vec, _MM_SWIZ_REG_CDAB );
    // out=(re - Ar Br, im - Ar Bi),...
    out.vec  = _mm512_fnmadd_ps(a_real, vb.vec, out.vec);
    // (re - Ar Br, im - Ar Bi) + (Bi,Br)*(-Ai, Ai) = (re - Ar Br - Bi Ai, im - Ar Bi + Br Ai)
    out.vec = _mm512_fmadd_ps( tmp0, a_imag, out.vec );
  }

#else // AVX512
  static inline vcomplexF_t cmult_conj(const vcomplexF_t &va, const vcomplexF_t &vb){
  // slight modification of a*b in Grid
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    return vcomplexF_t( _mm512_fmsubadd_ps( a_real, vb.vec, a_imag ) ); // (Ar Br, Ar Bi)+-(Ai Bi, Ai Br) = (ArBr+AiBi, ArBi-AiBr)
  }

  static void inline set_cmult_conj(vcomplexF_t &out, vcomplexF_t va, vcomplexF_t vb){
  // slight modification of a*b in Grid
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    out.vec =_mm512_fmsubadd_ps( a_real, vb.vec, a_imag ); // (Ar Br, Ar Bi)+-(Ai Bi, Ai Br) = (ArBr+AiBi, ArBi-AiBr)
  }

  static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // Ai Ai
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // Ar Ar
    a_imag = _mm512_fmsubadd_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ), out.vec );  // (Ai, Ai) * (Bi, Br) +-(Cr,Ci) = Ai Bi+Cr, Ai Br-Ci
    out.vec =_mm512_fmsubadd_ps( a_real, vb.vec, a_imag ); // (Ar Br, Ar Bi)+-(Ai Bi +Cr, Ai Br -Ci) = (ArBr+AiBi +Cr, ArBi-AiBr +Ci)
  }

  static inline void accum_sub_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // Ai Ai
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // Ar Ar
    a_imag = _mm512_fmaddsub_ps( a_imag, _mm512_permute_ps( vb.vec, 0xB1 ), out.vec );  // (Ai, Ai) * (Bi, Br) -+(Cr,Ci) = Ai Bi-Cr, Ai Br+Ci
    a_real = _mm512_sub_ps( _mm512_setzero_ps(), a_real); // ( -Ar -Ar)
    out.vec =_mm512_fmaddsub_ps( a_real, vb.vec, a_imag ); // (-Ar Br, -Ar Bi)-+(Ai Bi -Cr, Ai Br +Ci) = (-ArBr-AiBi +Cr, -ArBi+AiBr +Ci)
  }

#endif

#ifdef KNC
  static inline vcomplexF_t TimesI(const vcomplexF_t &in){
    __m512 tmp = _mm512_swizzle_ps(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    return vcomplexF_t( _mm512_mask_sub_ps(tmp, 0x5555, _mm512_setzero_ps(), tmp) ); // (-im, re)
  }
  static inline vcomplexF_t TimesMinusI(const vcomplexF_t &in){
    __m512 tmp = _mm512_swizzle_ps(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    return vcomplexF_t( _mm512_mask_sub_ps(tmp, 0xAAAA, _mm512_setzero_ps(), tmp) ); // (im, -re)
  }
  static inline void accum_TimesI(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp = _mm512_swizzle_ps(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    out.vec = _mm512_mask_add_ps(out.vec, 0xAAAA, out.vec,tmp); // (out_re,    out_im+re)
    out.vec = _mm512_mask_sub_ps(out.vec, 0x5555, out.vec,tmp); // (out_re-im, out_im+re)
  }

  static inline void accum_TimesI_compare(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp = _mm512_swizzle_ps(in.vec, _MM_SWIZ_REG_CDAB);  // (im, re)
    tmp=_mm512_mask_sub_ps(tmp, 0x5555, _mm512_setzero_ps(), tmp);
    out.vec = _mm512_add_ps(out.vec, tmp);
  }

  static inline void accum_TimesMinusI(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp = _mm512_swizzle_ps(in.vec, _MM_SWIZ_REG_CDAB);      // (im, re)
    out.vec = _mm512_mask_sub_ps(out.vec, 0xAAAA, out.vec,tmp); // (out_re,    out_im-re)
    out.vec = _mm512_mask_add_ps(out.vec, 0x5555, out.vec,tmp); // (out_re+im, out_im-re)
  }
#else
  // sqrt(-1)
  static inline vcomplexF_t TimesI(const vcomplexF_t &in){
    __m512 tmp = _mm512_shuffle_ps(in.vec, in.vec, 0xB1); // (im,re)
    return vcomplexF_t( _mm512_mask_sub_ps(tmp, 0x5555, _mm512_setzero_ps(), tmp) ); // (-im, re)
  }
  static inline vcomplexF_t TimesMinusI(const vcomplexF_t &in){
    __m512 tmp = _mm512_shuffle_ps(in.vec, in.vec, 0xB1); // (im,re)
    return vcomplexF_t( _mm512_mask_sub_ps(tmp, 0xAAAA, _mm512_setzero_ps(), tmp) ); // (im, -re)
  }

  static void inline accum_TimesI(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp=_mm512_shuffle_ps(in.vec,in.vec,IMM_QUART4(2,3,0,1));
    tmp=_mm512_mask_sub_ps(tmp,0x5555,_mm512_setzero_ps(),tmp);
    out.vec=_mm512_add_ps(out.vec, tmp);
  }
  static void inline accum_TimesI_compare(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp=_mm512_shuffle_ps(in.vec,in.vec,0xB1);    // (im, re)
    out.vec = _mm512_mask_add_ps(out.vec, 0xAAAA, out.vec, tmp); // (out_re,    out_im+re)
    out.vec = _mm512_mask_sub_ps(out.vec, 0x5555, out.vec, tmp); // (out_re-im, out_im+re)
  }

  static void inline accum_TimesMinusI(vcomplexF_t &out, const vcomplexF_t &in){
    __m512 tmp=_mm512_shuffle_ps(in.vec,in.vec,0xB1);
    tmp=_mm512_mask_sub_ps(tmp,0xAAAA,_mm512_setzero_ps(),tmp);
    out.vec=_mm512_add_ps(out.vec, tmp);
  }

#endif

  static inline vcomplexF_t Times_pm(const vcomplexF_t &in){
    // (-- ++ -- ++ -- ++ -- ++)
    //  11 00 11 00 11 00 11 00   = 0xCCCC
    __m512 tmp=_mm512_mask_sub_ps(in.vec, 0xCCCC, _mm512_setzero_ps(), in.vec);
    return vcomplexF_t(tmp);
  }

  static inline vcomplexF_t Times_mp(const vcomplexF_t &in){
    // (++ -- ++ -- ++ -- ++ --)
    //  00 11 00 11 00 11 00 11  = 0x3333
    __m512 tmp=_mm512_mask_sub_ps(in.vec, 0x3333, _mm512_setzero_ps(), in.vec);
    return vcomplexF_t(tmp);
  }

  // axpy: y := a*x+y
  static inline void axpy(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    accum_cmult(vy,va,vx);
  }
  static inline void axpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec = _mm512_fmadd_ps(va.vec, vx.vec, vy.vec);
  }

  static inline void axpy(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec = _mm512_fmadd_ps(va.vec, vx.vec, vy.vec);
  }

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

  // y := y-ax
  static inline void accum_sub_mult(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec = _mm512_fnmadd_ps( va.vec, vx.vec, vy.vec );
  }
  static inline void accum_sub_mult(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec = _mm512_fnmadd_ps( va.vec, vx.vec, vy.vec );
  }

  // aypx: y := a*y+x
#ifdef KNC
  static inline void aypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    //    vy.vec=va.vec*vy.vec + vx.vec;
    __m512 vzero = _mm512_setzero_ps();
    // va=(Ar Ai)  -->tmp0=(Ai Ar)
    __m512 tmp0  = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // a_real=(Ar Ar) :  0xAAAA =   (F T) (F T) (F T) (F T)....
    __m512 a_real = (__m512)_mm512_mask_or_epi32( (__m512i)va.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // a_imag = (-Ai Ai) : 0x5555 = (T F) (T F) (T F) (T F)....
    __m512 a_imag = _mm512_mask_sub_ps(va.vec, 0x5555, vzero, tmp0);
    // tmp1=(Yi Yr)
    __m512 tmp1 = _mm512_swizzle_ps( vy.vec, _MM_SWIZ_REG_CDAB );
    // tmp2=(Xr+Ar Yr, Xi+Ar Yi),...
    __m512 tmp2 = _mm512_fmadd_ps(a_real, vy.vec, vx.vec);
    // (Xr+Ar Br, Xi+Ar Bi) + (Yi,Yr)*(-Ai, Ai) = (Xr+Ar Yr - Yi Ai, Xi+Ar Yi + Yr Ai)
    vy.vec = _mm512_fmadd_ps( tmp1, a_imag, tmp2 );
  }
#else
  // aypx: y := a*y+x
  static inline void aypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    __m512 y_real = _mm512_moveldup_ps( vy.vec ); // Yr Yr
    __m512 y_imag = _mm512_movehdup_ps( vy.vec ); // Yi Yi
    vy.vec = _mm512_fmaddsub_ps( y_imag, _mm512_permute_ps( va.vec, 0xB1 ), vx.vec );
    vy.vec = _mm512_fmaddsub_ps( y_real, va.vec, vy.vec );
  }
#endif

  static inline void aypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec = _mm512_fmadd_ps( va.vec, vy.vec, vx.vec );
  }
  static inline void aypx(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    vy.vec = _mm512_fmadd_ps( va.vec, vy.vec, vx.vec );
  }

  // alias
  static inline void vaypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    aypx(vy, va, vx);
  }
  static inline void vaypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    aypx(vy, va, vx);
  }

  // scal
#ifdef KNC
  static inline void scal(vcomplexF_t &vx, const vcomplexF_t &va){
    __m512 vzero = _mm512_setzero_ps();
    // vx=(Xr, Xi) --> tmp0=(Xi, Xr)
    __m512 tmp0  = _mm512_swizzle_ps( vx.vec, _MM_SWIZ_REG_CDAB );
    //  x_real=(Xr Xr) :  0xAA =  (F T) (F T) (F T) (F T)
    __m512 x_real = (__m512)_mm512_mask_or_epi32( (__m512i)vx.vec, 0xAAAA,(__m512i)vzero, (__m512i)tmp0 );
    // x_imag = (-Xi, Xi)
    __m512 x_imag = _mm512_mask_sub_ps(vx.vec, 0x5555, vzero, tmp0);
    //tmp0=(Ai Ar)
    tmp0 = _mm512_swizzle_ps( va.vec, _MM_SWIZ_REG_CDAB );
    // vx =(Ar Xr, Ai Xr),...
    vx.vec  = _mm512_mul_ps(x_real, va.vec);
    // (Ar Xr, Ai Xr) + (-Xi,Xi)*(Ai, Ar) = (Ar Xr - Xi Ai, Ar Xi + Xr Ai)
    vx.vec = _mm512_fmadd_ps( tmp0, x_imag, vx.vec );
  }
#else // AVX512
  static inline void scal(vcomplexF_t &vx, const vcomplexF_t &va){
    __m512 a_real = _mm512_moveldup_ps( va.vec ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( va.vec ); // Ai Ai
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vx.vec, 0xB1 ) );
    vx.vec =_mm512_fmaddsub_ps( a_real, vx.vec, a_imag );
  }
#endif

  static inline void scal(vrealF_t &vx, const vrealF_t &va){
    vx.vec = _mm512_mul_ps(vx.vec, va.vec);
  }
  static inline void scal(vcomplexF_t &vx, const vrealF_t &va){
    vx.vec = _mm512_mul_ps(vx.vec, va.vec);
  }

  // (re im) (re im)... --> (re re) (re re)
#ifdef KNC
  static inline void complex2rere(vcomplexF_t &vx){
    __m512 tmp0  = _mm512_swizzle_ps( vx.vec, _MM_SWIZ_REG_CDAB );
    vx.vec = (__m512)_mm512_mask_or_epi32( (__m512i)vx.vec, 0xAAAA,(__m512i)_mm512_setzero_ps(), (__m512i)tmp0 );
    //    float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }
  }
#else // AVX512
  static inline void complex2rere(vcomplexF_t &vx){
    vx.vec=_mm512_moveldup_ps( vx.vec );
  }
#endif

  // (re im) (re im)... --> (im im) (im im)
#ifdef KNC
  static inline void complex2imim(vcomplexF_t &vx){
    __m512 tmp0  = _mm512_swizzle_ps( vx.vec, _MM_SWIZ_REG_CDAB );
    vx.vec = (__m512)_mm512_mask_or_epi32( (__m512i)vx.vec, 0x5555,(__m512i)_mm512_setzero_ps(), (__m512i)tmp0 );
    //float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i]=tmp[2*i+1];
    //    }
  }
#else // AVX512
  static inline void complex2imim(vcomplexF_t &vx){
    vx.vec=_mm512_movehdup_ps( vx.vec );
  }
#endif

  // (re im) (re im)... --> (re 0) (re 0)
  static inline void complex2re(vcomplexF_t &vx){
    vx.vec=_mm512_mask_mov_ps(vx.vec, 0xAAAA, _mm512_setzero_ps());
  }


  // (re im) (re im)... --> (0 im) (0 im)
  static inline void complex2im(vcomplexF_t &vx){
    vx.vec=_mm512_mask_mov_ps(vx.vec, 0x5555, _mm512_setzero_ps());
  }

  // TODO
  // (re1 im1),... (re1' im1')... --> (re1 re1') (re2 re2')...
  static inline void rere2complex(vcomplexF_t &vx1, const vcomplexF_t &vx2){
    vcomplexF_t tmp=vx2;
    complex2rere(tmp);
    vx1.vec=_mm512_mask_mov_ps(tmp.vec, 0x5555, vx1.vec);
    // other
    //float *tmp1=(float*)&vx1;
    //    float *tmp2=(float*)&vx2;
    //for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //  tmp1[2*i+1]=tmp2[2*i];
    //}
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
