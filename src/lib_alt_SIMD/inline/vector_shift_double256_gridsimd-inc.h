/*!

        @file    vector_shift_double256_gridsimd-inc.h

        @brief   shift in simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
/*

 Helper functions for AVX/AVX2

 assumes
  simd_layout={1,1,1,2}

*/

#ifndef VECTOR_SHIFT_DOUBLE256_GRIDSIMD_INC_H
#define VECTOR_SHIFT_DOUBLE256_GRIDSIMD_INC_H


#define SET_MASK(mask, flag)						\
  __m256d mask;								\
  if(flag){								\
    mask=_mm256_cmp_pd(_mm256_setzero_pd(),_mm256_setzero_pd(),_CMP_TRUE_US); \
  }else {								\
    mask=_mm256_setzero_pd();						\
  }


#ifdef AFOPR_HAS_SIMD_X
#undef AFOPR_HAS_SIMD_X
#endif
#ifdef AFOPR_HAS_SIMD_Y
#undef AFOPR_HAS_SIMD_Y
#endif
#ifdef AFOPR_HAS_SIMD_Z
#undef AFOPR_HAS_SIMD_Z
#endif
#ifdef AFOPR_HAS_SIMD_T
#undef AFOPR_HAS_SIMD_T
#endif

#define AFOPR_HAS_SIMD_T

namespace{

  
  //  using mask_t=bool;  // TODO: rewrite with available masks
  using mask_t = __m256d;
  
  static const bool has_simd_x=false;
  static const bool has_simd_y=false;
  static const bool has_simd_z=false;
  static const bool has_simd_t=true;


  ////////////////////////////////////////////////////////////////
  //
  // "shift" in simd
  //
  ///////////////////////////////////////////////////////////////

  // +x direction, -x direction: no simd
  static inline void shift_xp(vcomplexD_t &v){}
  static inline void shift_xp(vcomplexD_t &v, const mask_t &need_shift){}
  static inline void shift_xm(vcomplexD_t &v){}
  static inline void shift_xm(vcomplexD_t &v, const mask_t &need_shift){}

  // +y direction, -y direction: no simd
  static inline void shift_yp(vcomplexD_t &v){}
  static inline void shift_yp(vcomplexD_t &v, const mask_t &need_shift){}
  static inline void shift_ym(vcomplexD_t &v){}
  static inline void shift_ym(vcomplexD_t &v, const mask_t &need_shift){}

  // +z drection: no simd in double  
  static inline void shift_zp(vcomplexD_t &v){}
  static inline void shift_zp(vcomplexD_t &v, const mask_t &need_shift){}

  // -z direction: the same as +z
  static inline void shift_zm(vcomplexD_t &v){}
  static inline void shift_zm(vcomplexD_t &v, const mask_t &need_shift){}

  // +t direction
  static inline void shift_tp(vcomplexD_t &v){
    v.vec = (__m256d)_mm256_permute2f128_ps((__m256)v.vec, (__m256)v.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
  }
  static inline void shift_tp(vcomplexD_t &v, const mask_t &need_shift){
    __m256d tmp = (__m256d)_mm256_permute2f128_ps((__m256)v.vec, (__m256)v.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
    v.vec = _mm256_blendv_pd(v.vec, tmp, need_shift);
  }

  // -t direction: the same as +t
  static inline void shift_tm(vcomplexD_t &v){
    shift_tp(v);
  }
  static inline void shift_tm(vcomplexD_t &v, const mask_t &need_shift){
    shift_tp(v,need_shift);
  }

  //////////////////////////////////////////////////////////////////////
  //
  // packing simd to communication buffer
  //  
  static inline void pack_nosimd(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sstore(&buf[offset  ],sp1);
    sstore(&buf[offset+1],sp2);
  }


  static inline void pack_xp(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_xm(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_yp(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_ym(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_zp(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_zm(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }


  static inline void pack_tp(vcomplexD_t *buf, vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    sp1.vec=_mm256_permute2f128_pd(sp1.vec, sp2.vec, IMM_QUART4(0,2,0,0));
    sstore(buf+offset, sp1);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[1]=bb[0];
  }


  static inline void pack_tm(vcomplexD_t *buf, vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    sp1.vec=_mm256_permute2f128_pd(sp1.vec, sp2.vec, IMM_QUART4(0,1,0,3));
    sstore(buf+offset, sp1);
    
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[0]=bb[1];
    //    sstore(&buf[offset], sp1);
  }

  //////////////////////////////////////////////////////////////////////
  //
  // unpacking simd from communication buffer
  //  
  static inline void unpack_nosimd(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sp1=rbuf[offset];
    sp2=rbuf[offset+1];
  }

  static inline void unpack_xp(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_xm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_yp(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_ym(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_zp(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_zm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }


  static inline void unpack_tp(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){

    const int offset=offset_base+a;
    __m256d b=_mm256_load_pd((double*)(rbuf+offset));
    sp1.vec=_mm256_permute2f128_pd(sp1.vec, b, IMM_QUART4(0,1,0,2));
    sp2.vec=_mm256_permute2f128_pd(sp2.vec, b, IMM_QUART4(0,1,0,3));

    // vComplexD tmp=rbuf[offset];
    // ComplexD *p =(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[0]=p[0];
    // bb[0]=p[1];
  }

  static inline void unpack_tm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m256d b=_mm256_load_pd((double*)(rbuf+offset));
    sp1.vec=_mm256_permute2f128_pd(sp1.vec, b, IMM_QUART4(0,3,0,0));
    sp2.vec=_mm256_permute2f128_pd(sp2.vec, b, IMM_QUART4(0,2,0,0));

    // vComplexD tmp=rbuf[offset];
    // ComplexD *p =(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[1]=p[1];
    // bb[1]=p[0];
  }

  // boundary condition: multiply a real number
  static inline void apply_boundary(vcomplexD_t &val, const vrealD_t &boundary_value){
    val.vec = _mm256_mul_pd(val.vec, boundary_value.vec);
  }

  static inline void set_boundary_xp(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_xm(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_yp(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_ym(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_zp(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_zm(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_tp(vrealD_t &boundary_value, const double rval){
    //    bcast(boundary_value,rval);
    vout.crucial("set_bourndary_tp: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_tm(vrealD_t &boundary_value, const double rval){
    //    bcast(boundary_value,rval);
    vout.crucial("set_bourndary_tm: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  // dummy
  static inline void shift_and_merge_xp(vcomplexD_t &sp1, vcomplexD_t &sp2){
  }
  static inline void shift_and_merge_xm(vcomplexD_t &sp1, vcomplexD_t &sp2){
  }

  // boundary condition: mulpliy a complex number
  // not implemented


  
}

#endif
