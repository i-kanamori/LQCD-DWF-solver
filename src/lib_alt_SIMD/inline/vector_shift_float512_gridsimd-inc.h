/*!

        @file    vector_shift_float512_gridsimd-inc.h

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
  simd_layout={1,2,2,2}

*/

#ifndef VECTOR_SHIFT_FLOAT512_GRIDSIMD_INC_H
#define VECTOR_SHIFT_FLOAT512_GRIDSIMD_INC_H

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

#define AFOPR_HAS_SIMD_Y
#define AFOPR_HAS_SIMD_Z
#define AFOPR_HAS_SIMD_T

namespace{

  using mask_t=__mmask16;
  
  static const bool has_simd_x=false;
  static const bool has_simd_y=true;
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

  // +y direction
  static inline void shift_yp(vcomplexF_t &v){
#ifdef KNC
    v.vec = _mm512_swizzle_ps(v.vec, _MM_SWIZ_REG_BADC);
#else
    v.vec = _mm512_shuffle_ps(v.vec, v.vec, IMM_QUART4(1,0,3,2));
#endif
  }
  static inline void shift_yp(vcomplexF_t &v, const mask_t &need_shift){
#ifdef KNC
    v.vec = _mm512_mask_swizzle_ps(v.vec, need_shift, v.vec, _MM_SWIZ_REG_BADC);
#else
    v.vec = _mm512_mask_shuffle_ps(v.vec, need_shift, v.vec, v.vec, IMM_QUART4(1,0,3,2));
#endif
  }

  // -y direction: the same as +y
  static void inline shift_ym(vcomplexF_t &v){
    shift_yp(v);
  }
  static void inline shift_ym(vcomplexF_t &v, const mask_t &need_shift){
    shift_yp(v, need_shift);
  }

  // +z drection
  static inline void shift_zp(vcomplexF_t &v){
#ifdef KNC
    v.vec = _mm512_permute4f128_ps(v.vec, (_MM_PERM_ENUM)IMM_QUART4(2,3,0,1));
#else
    v.vec = _mm512_shuffle_f32x4(v.vec, v.vec,IMM_QUART4(2,3,0,1));
#endif
  }
  static inline void shift_zp(vcomplexF_t &v, const mask_t &need_shift){
#ifdef KNC
    v.vec = _mm512_mask_permute4f128_ps(v.vec, need_shift, v.vec, (_MM_PERM_ENUM)IMM_QUART4(2,3,0,1));
#else
    v.vec = _mm512_mask_shuffle_f32x4(v.vec, need_shift, v.vec, v.vec,IMM_QUART4(2,3,0,1));
#endif
  }

  // -z direction: the same as +z
  static inline void shift_zm(vcomplexF_t &v){
    shift_zp(v);
  }
  static inline void shift_zm(vcomplexF_t &v, const mask_t &need_shift){
    shift_zp(v,need_shift);
  }

  // +t direction
  static inline void shift_tp(vcomplexF_t &v){
    v.vec=(__m512)_mm512_alignr_epi32((__m512i)v.vec,(__m512i)v.vec, 8);
  }
  static inline void shift_tp(vcomplexF_t &v, const mask_t &need_shift){
    v.vec=(__m512)_mm512_mask_alignr_epi32((__m512i)v.vec, need_shift, (__m512i)v.vec,(__m512i)v.vec, 8);
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
    const int offset=2*(offset_base+a);
    sstore(&buf[offset  ],sp1);
    sstore(&buf[offset+1],sp2);
  }

  static inline void pack_xp(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_xm(vcomplexF_t *buf, const vcomplexF_t &sp1, const vcomplexF_t &sp2,  const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_yp(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2, const int offset_base, const int a=0){

    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0xCCCC, (__m512i)sp2.vec,(__m512i)sp2.vec,14);          
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[1]=bb[0];
    //    aa[3]=bb[2];
    //    aa[5]=bb[4];
    //    aa[7]=bb[6];
    sstore(buf+offset, sp1);
  }

  static inline void pack_ym(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2,  const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0x3333, (__m512i)sp2.vec,(__m512i)sp2.vec,2);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[0]=bb[1];
    //    aa[2]=bb[3];
    //    aa[4]=bb[5];
    //    aa[6]=bb[7];
    sstore(buf+offset, sp1);
  }

  static inline void pack_zp(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0xF0F0, (__m512i)sp2.vec,(__m512i)sp2.vec,12);          
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[2]=bb[0];
    //    aa[3]=bb[1];
    //    aa[6]=bb[4];
    //    aa[7]=bb[5];
    sstore(buf+offset, sp1);
 }


  static inline void pack_zm(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0x0F0F, (__m512i)sp2.vec,(__m512i)sp2.vec,4);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[0]=bb[2];
    //    aa[1]=bb[3];
    //    aa[4]=bb[6];
    //    aa[5]=bb[7];
    sstore(buf+offset, sp1);
  }

  static inline void pack_tp(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0xFF00, (__m512i)sp2.vec,(__m512i)sp2.vec,8);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[4]=bb[0];
    //    aa[5]=bb[1];
    //    aa[6]=bb[2];
    //    aa[7]=bb[3];
    sstore(buf+offset, sp1);
  }


  static inline void pack_tm(vcomplexF_t *buf, vcomplexF_t &sp1, vcomplexF_t &sp2, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    sp1.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp1.vec, 0x00FF, (__m512i)sp2.vec,(__m512i)sp2.vec,8);
    //    complexF_t *aa=(complexF_t*)&sp1;
    //    complexF_t *bb=(complexF_t*)&sp2;
    //    aa[0]=bb[4];
    //    aa[1]=bb[5];
    //    aa[2]=bb[6];
    //    aa[3]=bb[7];
    sstore(buf+offset, sp1);
  }

  //////////////////////////////////////////////////////////////////////
  //
  // unpacking simd from communication buffer
  //
  static inline void unpack_nosimd(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=2*(offset_base+a);
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
    const int offset=offset_base+a;
    __m512 b=_mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0x3333, b);
    sp2.vec  = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0x3333, (__m512i)b,(__m512i)b,2);

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p= (complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    // aa[0]=p[0];
    // aa[2]=p[2];
    // aa[4]=p[4];
    // aa[6]=p[6];
    // bb[0]=p[1];
    // bb[2]=p[3];
    // bb[4]=p[5];
    // bb[6]=p[7];

  }
  static inline void unpack_ym(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m512 b=_mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0xCCCC, b);
    sp2.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0xCCCC, (__m512i)b,(__m512i)b,14);

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p= (complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    // aa[1]=p[1];
    // aa[3]=p[3];
    // aa[5]=p[5];
    // aa[7]=p[7];
    // bb[1]=p[0];
    // bb[3]=p[2];
    // bb[5]=p[4];
    // bb[7]=p[6];
  }

  static inline void unpack_zp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){

    const int offset=offset_base+a;
    __m512 b=_mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0x0F0F, b);
    sp2.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0x0F0F, (__m512i)b,(__m512i)b,4); 

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p= (complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    // aa[0]=p[0];
    // aa[1]=p[1];
    // aa[4]=p[4];
    // aa[5]=p[5];
    // bb[0]=p[2];
    // bb[1]=p[3];
    // bb[4]=p[6];
    // bb[5]=p[7];
  }

  static inline void unpack_zm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m512 b=_mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0xF0F0, b);
    sp2.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0xF0F0, (__m512i)b,(__m512i)b,12);

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p=(complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    // aa[2]=p[2];
    // aa[3]=p[3];
    // aa[6]=p[6];
    // aa[7]=p[7];
    // bb[2]=p[0];
    // bb[3]=p[1];
    // bb[6]=p[4];
    // bb[7]=p[5];
  }


  static inline void unpack_tp(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;

    __m512 b = _mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0x00FF, b);
    sp2.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0x00FF, (__m512i)b,(__m512i)b,8);

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p=(complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    //
    // aa[0]=p[0];
    // aa[1]=p[1];
    // aa[2]=p[2];
    // aa[3]=p[3];
    // bb[0]=p[4];
    // bb[1]=p[5];
    // bb[2]=p[6];
    // bb[3]=p[7];
  }

  static inline void unpack_tm(vcomplexF_t &sp1, vcomplexF_t &sp2, const vcomplexF_t *rbuf, const int offset_base, const int a=0){
    const int offset=offset_base+a;
    __m512 b = _mm512_load_ps((float*)(rbuf+offset));
    sp1.vec = _mm512_mask_mov_ps(sp1.vec, 0xFF00, b);
    sp2.vec = (__m512)_mm512_mask_alignr_epi32((__m512i)sp2.vec, 0xFF00, (__m512i)b,(__m512i)b,8);

    // vcomplexF_t tmp=rbuf[offset];
    // complexF_t *p=(complexF_t*)&tmp;
    // complexF_t *aa=(complexF_t*)&sp1;
    // complexF_t *bb=(complexF_t*)&sp2;
    // aa[4]=p[4];
    // aa[5]=p[5];
    // aa[6]=p[6];
    // aa[7]=p[7];
    // bb[4]=p[0];
    // bb[5]=p[1];
    // bb[6]=p[2];
    // bb[7]=p[3];
  }

  // boundary condition: mulpliy a real number
  static inline void apply_boundary(vcomplexF_t &val, const vrealF_t &boundary_value){
    val.vec = _mm512_mul_ps(val.vec, boundary_value.vec);
  }

  static inline void set_boundary_xp(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_xm(vrealF_t &boundary_value, const float rval){
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_yp(vrealF_t &boundary_value, const float rval){
    vout.crucial("set_boundary_yp: not yet ready\n");
    exit(EXIT_FAILURE);
  }

  static inline void set_boundary_ym(vrealF_t &boundary_value, const float rval){
    vout.crucial("set_boundary_ym: not yet ready\n");
    exit(EXIT_FAILURE);
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
