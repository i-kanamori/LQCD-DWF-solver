/*!

        @file    vector_shift_double512_gridsimd-inc.h

        @brief   shift in simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
/*

 Helper functions for AVX512

 assumes
  simd_layout={1,1,2,2}

*/


#ifndef VECTOR_SHIFT_DOUBLE512_GRIDSIMD_INC_H
#define VECTOR_SHIFT_DOUBLE512_GRIDSIMD_INC_H

#define SET_MASK(mask, flag) const __mmask16 mask = (flag?0xffff:0x0000);

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

  using mask_t=__mmask16;
  
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
  static inline void shift_xp(vcomplexD_t &v){}
  static inline void shift_xp(vcomplexD_t &v, const mask_t &need_shift){}
  static inline void shift_xm(vcomplexD_t &v){}
  static inline void shift_xm(vcomplexD_t &v, const mask_t &need_shift){}

  // +y direction, -y direction: no simd
  static inline void shift_yp(vcomplexD_t &v){}
  static inline void shift_yp(vcomplexD_t &v, const mask_t &need_shift){}
  static inline void shift_ym(vcomplexD_t &v){}
  static inline void shift_ym(vcomplexD_t &v, const mask_t &need_shift){}

  // +z drection
  static inline void shift_zp(vcomplexD_t &v){
#ifdef KNC
    v.vec = _mm512_swizzle_pd(v.vec,_MM_SWIZ_REG_BADC);
#else
    v.vec = _mm512_shuffle_f64x2(v.vec, v.vec,IMM_QUART4(2,3,0,1));
#endif
  }
  static inline void shift_zp(vcomplexD_t &v, const mask_t &need_shift){
#ifdef KNC
    v.vec = _mm512_mask_swizzle_pd(v.vec, need_shift, v.vec ,_MM_SWIZ_REG_BADC);
#else
    v.vec = _mm512_mask_shuffle_f64x2(v.vec, need_shift, v.vec, v.vec, IMM_QUART4(2,3,0,1));
#endif

  }

  // -z direction: the same as +z
  static inline void shift_zm(vcomplexD_t &v){
    shift_zp(v);
  }
  static inline void shift_zm(vcomplexD_t &v, const mask_t &need_shift){
    shift_zp(v,need_shift);
  }

  // +t direction
  static inline void shift_tp(vcomplexD_t &v){
    v.vec = (__m512d)_mm512_alignr_epi32((__m512i)v.vec, (__m512i)v.vec, 8);
  }
  static inline void shift_tp(vcomplexD_t &v, const mask_t &need_shift){
    v.vec = (__m512d)_mm512_mask_alignr_epi32((__m512i)v.vec, need_shift, (__m512i)v.vec, (__m512i)v.vec, 8);
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
    const int offset=2*(offset_base+a);
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

  static inline void pack_zp(vcomplexD_t *buf, vcomplexD_t &sp1, vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=(__m512d)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0xF0F0, (__m512i)sp2.vec,(__m512i)sp2.vec,12);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[1]=bb[0];
    //    aa[3]=bb[2];
    sstore(buf+offset, sp1);
 }


  static inline void pack_zm(vcomplexD_t *buf, vcomplexD_t &sp1, vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=(__m512d)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0x0F0F, (__m512i)sp2.vec,(__m512i)sp2.vec,4);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[0]=bb[1];
    //    aa[2]=bb[3];
    sstore(buf+offset, sp1);
  }

  static inline void pack_tp(vcomplexD_t *buf, vcomplexD_t &sp1, vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=(__m512d)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0xFF00, (__m512i)sp2.vec,(__m512i)sp2.vec,8);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[2]=bb[0];
    //    aa[3]=bb[1];
    sstore(buf+offset, sp1);
 }


  static inline void pack_tm(vcomplexD_t *buf, vcomplexD_t &sp1, vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec=(__m512d)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0x00FF, (__m512i)sp2.vec,(__m512i)sp2.vec,8);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    aa[0]=bb[2];
    //    aa[1]=bb[3];
    sstore(buf+offset, sp1);
  }

  //////////////////////////////////////////////////////////////////////
  //
  // unpacking simd from communication buffer
  //  
  static inline void unpack_nosimd(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=2*(offset_base+a);
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
    const int offset=offset_base+a;

    __m512d b=_mm512_load_pd((double*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_pd(sp1.vec, 0x33, b);
    sp2.vec =(__m512d)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0x0F0F,(__m512i)b,(__m512i)b,4);

    // vComplexD tmp=rbuf[offset];
    // ComplexD *p=(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[0]=p[0];
    // aa[2]=p[2];
    // bb[0]=p[1];
    // bb[2]=p[3];
  }

  static inline void unpack_zm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    __m512d b=_mm512_load_pd((double*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_pd(sp1.vec, 0xCC, b);
    sp2.vec=(__m512d)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0xF0F0, (__m512i)b, (__m512i)b,12);          
    // vComplexD tmp=rbuf[offset];
    // ComplexD *p=(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[1]=p[1];
    // aa[3]=p[3];
    // bb[1]=p[0];
    // bb[3]=p[2];

  }


  static inline void unpack_tp(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    __m512d b=_mm512_load_pd((double*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_pd(sp1.vec, 0x0F, b);
    sp2.vec = (__m512d)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0x00FF, (__m512i)b,(__m512i)b,8);          
    // vComplexD tmp=rbuf[offset];
    // ComplexD *p=(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[0]=p[0];
    // aa[1]=p[1];
    // bb[0]=p[2];
    // bb[1]=p[3];
  }

  static inline void unpack_tm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    __m512d b =_mm512_load_pd((double*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_pd(sp1.vec, 0xF0, b);
    sp2.vec = (__m512d)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0xFF00, (__m512i)b,(__m512i)b,8);       

    // vComplexD tmp=rbuf[offset];
    // ComplexD *p=(ComplexD*)&tmp;
    // ComplexD *aa=(ComplexD*)&sp1;
    // ComplexD *bb=(ComplexD*)&sp2;
    // aa[2]=p[2];
    // aa[3]=p[3];
    // bb[2]=p[0];
    // bb[3]=p[1];
  }



  // boundary condition: multiply a real number
  static inline void apply_boundary(vcomplexD_t &val, const vrealD_t &boundary_value){
    val.vec = _mm512_mul_pd(val.vec, boundary_value.vec);
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
    vout.crucial("set_bourndary_zp: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_zm(vrealD_t &boundary_value, const double rval){
    vout.crucial("set_bourndary_zm: not yet ready\n");
    exit(EXIT_FAILURE);
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
