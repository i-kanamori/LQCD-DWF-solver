/*!

        @file    vector_shift_float256_impl1-inc.h

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
  simd_layout={4,1,1,1}
  Ncol=3

*/

#ifndef VECTOR_SHIFT_FLOAT256_IMPL1_INC_H
#define VECTOR_SHIFT_FLOAT256_IMPL1_INC_H


#define SET_MASK(mask, flag) const bool &mask = flag;

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

#define AFOPR_HAS_SIMD_X


namespace{

  using mask_t=bool;  // irreleavant for impl1

  static const bool has_simd_x=true;
  static const bool has_simd_y=false;
  static const bool has_simd_z=false;
  static const bool has_simd_t=false;


  ////////////////////////////////////////////////////////////////
  //
  // "shift" in simd
  //
  ///////////////////////////////////////////////////////////////

  // +x direction, -x direction
  static inline void shift_and_merge_xp(vcomplexF_t &sp1, vcomplexF_t &sp2){
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    complexF_t tmp=bb[0];  // needed for aa==bb case
    //    aa[0]=aa[1];
    //    aa[1]=aa[2];
    //    aa[2]=aa[3];
    //    aa[3]=tmp;

    //  a0 a1 a2 a3, b0 b1 b2 b3 --> a1 a2 a3 b0
    __m256 tmp=_mm256_blend_ps(sp1.vec, sp2.vec,IMM_BIN8(0,0,0,0,0,0,1,1));
    sp1.vec=(__m256)_mm256_permute4x64_epi64((__m256i)tmp, IMM_QUART4(0,3,2,1));    // 1234 5678 --> 3456 7812
  }

  static inline void shift_and_merge_xm(vcomplexF_t &sp1, vcomplexF_t &sp2){
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    complexF_t tmp=bb[3];  // needed for aa==bb case
    //    aa[3]=aa[2];
    //    aa[2]=aa[1];
    //    aa[1]=aa[0];
    //    aa[0]=tmp;

    //  b0 b1 b2 b3, a0 a1 a2 a3  --> b3 a0 a1 a2
    __m256 tmp=_mm256_blend_ps(sp1.vec, sp2.vec,IMM_BIN8(1,1,0,0,0,0,0,0));
    sp1.vec=(__m256)_mm256_permute4x64_epi64((__m256i)tmp, IMM_QUART4(2,1,0,3));
  }


  // +x direction, -x direction: dummy
  static inline void shift_xp(vcomplexF_t &v){}
  static inline void shift_xp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_xm(vcomplexF_t &v){}
  static inline void shift_xm(vcomplexF_t &v, const mask_t &need_shift){}

  // +y direction, -y direction: no simd
  static inline void shift_yp(vcomplexF_t &v){}
  static inline void shift_yp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_ym(vcomplexF_t &v){}
  static inline void shift_ym(vcomplexF_t &v, const mask_t &need_shift){}

  // +z direction, -z direction: no simd
  static inline void shift_zp(vcomplexF_t &v){}
  static inline void shift_zp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_zm(vcomplexF_t &v){}
  static inline void shift_zm(vcomplexF_t &v, const mask_t &need_shift){}

  // +t direction, -t direction: no simd
  static inline void shift_tp(vcomplexF_t &v){}
  static inline void shift_tp(vcomplexF_t &v, const mask_t &need_shift){}
  static inline void shift_tm(vcomplexF_t &v){}
  static inline void shift_tm(vcomplexF_t &v, const mask_t &need_shift){}


  //////////////////////////////////////////////////////////////////////
  //
  // packing simd to communication buffer
  //
  static inline void pack_nosimd(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sstore(&buf[offset  ],sp1);
    sstore(&buf[offset+1],sp2);
  }

  //  static inline void pack_xp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a){
  //    pack_nosimd(buf,sp1,sp2,offset_base,a);
  //  }

  //  static inline void pack_xm(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a){
  //    pack_nosimd(buf,sp1,sp2,offset_base,a);
  //  }

  static inline void pack_yp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_ym(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_zp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_zm(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_tp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_tm(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
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

  //  static inline void unpack_xp(vcomplexF_t &sp1, vcomplexF_t &sp2, vcomplexF_t *rbuf, const int offset_base, const int a){
  //    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  //  }
  //  static inline void unpack_xm(vcomplexF_t &sp1, vcomplexF_t &sp2, vcomplexF_t *rbuf, const int offset_base, const int a){
  //    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  //  }

  static inline void unpack_yp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_ym(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_zp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_zm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_tp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }
  static inline void unpack_tm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }


  // boundary condition: mulpliy a real number
  static inline void apply_boundary(vcomplexF_t &val, const vrealF_t &boundary_value){
    val.vec = _mm256_mul_ps(val.vec, boundary_value.vec);
  }

  static inline void set_boundary_xp(vrealF_t &boundary_value, const float rval){
    boundary_value.vec = _mm256_set_ps(rval,rval,1,1,1,1,1,1);
  }

  static inline void set_boundary_xm(vrealF_t &boundary_value, const float rval){
    boundary_value.vec = _mm256_set_ps(1,1,1,1,1,1,rval,rval);
  }

  static inline void set_boundary_yp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_ym(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_zp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_zm(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_tp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_tm(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  // boundary condition: mulpliy a complex number
  // not implemented

} // end of anonymous namespace

#endif
