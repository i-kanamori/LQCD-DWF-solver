/*!

        @file    vector_shift_float256_gridsimd-inc.h

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
  simd_layout={1,1,2,2}

*/

#ifndef VECTOR_SHIFT_FLOAT256_GRIDSIMD_INC_H
#define VECTOR_SHIFT_FLOAT256_GRIDSIMD_INC_H


#define SET_MASK(mask, flag)						\
  __m256 mask;								\
  if(flag){								\
    mask=_mm256_cmp_ps(_mm256_setzero_ps(),_mm256_setzero_ps(),_CMP_TRUE_US); \
  }else {								\
    mask=_mm256_setzero_ps();						\
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

#define AFOPR_HAS_SIMD_Z
#define AFOPR_HAS_SIMD_T


namespace{

  //  using mask_t=bool;  // TODO: rewrite with available masks
  using mask_t=__m256;
  
  static const bool has_simd_x=false;
  static const bool has_simd_y=false;
  static const bool has_simd_z=true;
  static const bool has_simd_t=true;


  ////////////////////////////////////////////////////////////////
  //
  // "shift" in simd
  //
  ///////////////////////////////////////////////////////////////

  // +x direction, -x direction: no simd
  static inline void shift_xp(vcomplexF_t &v){}
  static inline void shift_xp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_xm(vcomplexF_t &v){}
  static inline void shift_xm(vcomplexF_t &v, const mask_t &need_shift){}
  
  // +y direction, -y direction: no simd
  static inline void shift_yp(vcomplexF_t &v){}
  static inline void shift_yp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_ym(vcomplexF_t &v){}
  static inline void shift_ym(vcomplexF_t &v, const mask_t &need_shift){}
			      
  // +z drection: 
  static inline void shift_zp(vcomplexF_t &v){
    v.vec=_mm256_shuffle_ps(v.vec,v.vec,IMM_QUART4(1,0,3,2)); 
  }
  static inline void shift_zp(vcomplexF_t &v, const mask_t &need_shift){
    __m256 tmp = _mm256_shuffle_ps(v.vec,v.vec,IMM_QUART4(1,0,3,2)); 
    v.vec = _mm256_blendv_ps(v.vec, tmp, need_shift);
  }

  // -z direction: the same as +z
  static inline void shift_zm(vcomplexF_t &v){
    shift_zp(v);
  }
  static inline void shift_zm(vcomplexF_t &v, const mask_t &need_shift){
    shift_zp(v, need_shift);
  }

  // +t direction
  static inline void shift_tp(vcomplexF_t &v){
    v.vec = _mm256_permute2f128_ps(v.vec, v.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
  }
  static inline void shift_tp(vcomplexF_t &v, const mask_t &need_shift){
    __m256 tmp = _mm256_permute2f128_ps(v.vec, v.vec, IMM_QUART4(0,0,0,1));  // 1234 5678 --> 5678 1234
    v.vec = _mm256_blendv_ps(v.vec, tmp, need_shift);
  }

  // -t direction: the same as +t
  static inline void shift_tm(vcomplexF_t &v){
    shift_tp(v);
  }
  static inline void shift_tm(vcomplexF_t &v, const mask_t &need_shift){
    shift_tp(v,need_shift);
  }

  //////////////////////////////////////////////////////////////////////
  //
  // packing simd to communication buffer
  //  
  static inline void pack_nosimd(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sstore(&buf[offset  ],sp1);
    sstore(&buf[offset+1],sp2);
  }

  static inline void pack_xp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_xm(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_yp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_ym(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_zp(vcomplexF_t *buf, vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=_mm256_shuffle_ps(sp1.vec, sp2.vec, IMM_QUART4(1,0,1,0));
    sstore(buf+offset, sp1);
    //    ComplexF_t *aa=(ComplexF_t*)&sp1;
    //    ComplexF_t *bb=(ComplexF_t*)&sp2;
    //    aa[1]=bb[0];
    //    aa[3]=bb[2];
  }

  static inline void pack_zm(vcomplexF_t *buf, vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=_mm256_shuffle_ps(sp2.vec, sp1.vec, IMM_QUART4(3,2,3,2));
    sstore(buf+offset, sp1);
    //    ComplexF_t *aa=(ComplexF_t*)&sp1;
    //    ComplexF_t *bb=(ComplexF_t*)&sp2;
    //    aa[0]=bb[1];
    //    aa[2]=bb[3];
  }


  static inline void pack_tp(vcomplexF_t *buf, vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    sp1.vec=_mm256_permute2f128_ps(sp1.vec, sp2.vec, IMM_QUART4(0,2,0,0));
    sstore(buf+offset, sp1);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[2]=bb[0];
    //    aa[3]=bb[1];

  }


  static inline void pack_tm(vcomplexF_t *buf, vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    sp1.vec=_mm256_permute2f128_ps(sp1.vec, sp2.vec, IMM_QUART4(0,1,0,3));
    sstore(buf+offset, sp1);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[0]=bb[2];
    //    aa[1]=bb[3];
  }

  //////////////////////////////////////////////////////////////////////
  //
  // unpacking simd from communication buffer
  //  
  static inline void unpack_nosimd(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sp1=rbuf[offset];
    sp2=rbuf[offset+1];
  }

  static inline void unpack_xp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_xm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_yp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_ym(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_zp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m256 b=_mm256_load_ps((float*)(rbuf+offset));
    sp1.vec=_mm256_shuffle_ps(b, sp1.vec, IMM_QUART4(3,2,1,0));
    sp2.vec=_mm256_shuffle_ps(b, sp2.vec, IMM_QUART4(3,2,3,2));

    // vComplexF tmp=rbuf[offset];
    // ComplexF *p =(ComplexF*)&tmp;
    // ComplexF *aa=(ComplexF*)&sp1;
    // ComplexF *bb=(ComplexF*)&sp2;
    // aa[0]=p[0];  // a:63--0, p:63--0
    // aa[2]=p[2];  // a:(128+)63--0, p:(128+)63--0
    // bb[0]=p[1];  // b:63--0, p:127--64
    // bb[2]=p[3];  // a:(128+)63--0, p:(128+)127--64

  }
  static inline void unpack_zm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m256 b=_mm256_load_ps((float*)(rbuf+offset));
    sp1.vec=_mm256_shuffle_ps(sp1.vec, b, IMM_QUART4(3,2,1,0));
    sp2.vec=_mm256_shuffle_ps(sp2.vec, b, IMM_QUART4(1,0,1,0));

    // vComplexF tmp=rbuf[offset];
    // ComplexF *p =(ComplexF*)&tmp;
    // ComplexF *aa=(ComplexF*)&sp1;
    // ComplexF *bb=(ComplexF*)&sp2;
    // aa[1]=p[1];
    // aa[3]=p[3];
    // bb[1]=p[0];
    // bb[3]=p[2];
  }


  static inline void unpack_tp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){

    const int offset=offset_base+a;
    __m256 b=_mm256_load_ps((float*)(rbuf+offset));
    sp1.vec=_mm256_permute2f128_ps(sp1.vec, b, IMM_QUART4(0,1,0,2));
    sp2.vec=_mm256_permute2f128_ps(sp2.vec, b, IMM_QUART4(0,1,0,3));

    // vComplexF tmp=rbuf[offset];
    // ComplexF *p =(ComplexF*)&tmp;
    // ComplexF *aa=(ComplexF*)&sp1;
    // ComplexF *bb=(ComplexF*)&sp2;
    // aa[0]=p[0];
    // aa[1]=p[1];
    // bb[0]=p[2];
    // bb[1]=p[3];
  }

  static inline void unpack_tm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){

    const int offset=offset_base+a;
    __m256 b=_mm256_load_ps((float*)(rbuf+offset));
    sp1.vec=_mm256_permute2f128_ps(sp1.vec, b, IMM_QUART4(0,3,0,0));
    sp2.vec=_mm256_permute2f128_ps(sp2.vec, b, IMM_QUART4(0,2,0,0));

    // vComplexF tmp=rbuf[offset];
    // ComplexF *p =(ComplexF*)&tmp;
    // ComplexF *aa=(ComplexF*)&sp1;
    // ComplexF *bb=(ComplexF*)&sp2;
    // aa[2]=p[2];
    // aa[3]=p[3];
    // bb[2]=p[0];
    // bb[3]=p[1];
  }

  // boundary condition: mulpliy a real number
  static inline void apply_boundary(vcomplexF_t &val, const vrealF_t &boundary_value){
    val.vec = _mm256_mul_ps(val.vec, boundary_value.vec);
  }

  static inline void set_boundary_xp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_xm(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_yp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_ym(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_zp(vrealF_t &boundary_value, const float rval){
    vout.crucial("set_boundary_zp: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_zm(vrealF_t &boundary_value, const float rval){
    vout.crucial("set_boundary_zm: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_tp(vrealF_t &boundary_value, const float rval){
    vout.crucial("set_boundary_tp: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_tm(vrealF_t &boundary_value, const float rval){ 
    vout.crucial("set_boundary_tm: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  // dummy
  static inline void shift_and_merge_xp(vcomplexF_t &sp1, vcomplexF_t &sp2){
  }
  static inline void shift_and_merge_xm(vcomplexF_t &sp1, vcomplexF_t &sp2){
  }

  // boundary condition: mulpliy a complex number
  // not implemented

} // end of anonymous namespace

#endif
