/*!

        @file    vector_shift_double512_impl1-inc.h

        @brief   shift in simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
/*
 Helper functions for AVX512 / KNC

 assumes
  simd_layout={4,1,1,1}


*/

#ifndef VECTOR_SHIFT_DOUBLE512_IMPL1_INC_H
#define VECTOR_SHIFT_DOUBLE512_IMPL1_INC_H

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

#define AFOPR_HAS_SIMD_X

namespace{

  using mask_t=__mmask16;

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
  static inline void shift_and_merge_xp(vcomplexD_t &sp1, vcomplexD_t &sp2){
    //  a0 a1 a2 a3, b0 b1 b2 b3 --> a1 a2 a3 b0
    sp1.vec=(__m512d)_mm512_alignr_epi32((__m512i)sp2.vec,(__m512i)sp1.vec, 4);  
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    complexD_t  tmp=bb[0];  // needed for aa==bb case
    //    aa[0]=aa[1];
    //    aa[1]=aa[2];
    //    aa[2]=aa[3];
    //    aa[3]=tmp;

  }
  static inline void shift_and_merge_xm(vcomplexD_t &sp1, vcomplexD_t &sp2){
    //  b0 b1 b2 b3, a0 a1 a2 a3  --> b3 a0 a1 a2
    sp1.vec=(__m512d)_mm512_alignr_epi32((__m512i)sp1.vec,(__m512i)sp2.vec, 12);
    //    complexD_t *aa=(complexD_t*)&sp1;
    //    complexD_t *bb=(complexD_t*)&sp2;
    //    complexD_t  tmp=bb[3];  // needed for aa==bb case
    //    aa[3]=aa[2];
    //    aa[2]=aa[1];
    //    aa[1]=aa[0];
    //    aa[0]=tmp;
  }

  // +x direction, -x direction: dummy
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

  // +t drection: no simd in double
  static inline void shift_tp(vcomplexD_t &v){}
  static inline void shift_tp(vcomplexD_t &v, const mask_t &need_shift){}

  // -t direction: the same as +t
  static inline void shift_tm(vcomplexD_t &v){}
  static inline void shift_tm(vcomplexD_t &v, const mask_t &need_shift){}


  //////////////////////////////////////////////////////////////////////
  //
  // packing simd to communication buffer
  //
  static inline void pack_nosimd(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    const int offset=2*offset_base+2*a;
    sstore(&buf[offset  ],sp1);
    sstore(&buf[offset+1],sp2);
  }


  //  static inline void pack_xp(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
  //    pack_nosimd(buf,sp1,sp2,offset_base,a);
  //  }

  //  static inline void pack_xm(vcomplexD_t *buf, const vcomplexD_t &sp1, const vcomplexD_t &sp2,  const int offset_base, const int a=0){
  //    pack_nosimd(buf,sp1,sp2,offset_base,a);
  //  }

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
    pack_nosimd(buf,sp1,sp2,offset_base,a);
  }

  static inline void pack_tm(vcomplexD_t *buf, vcomplexD_t &sp1, const vcomplexD_t &sp2, const int offset_base, const int a=0){
    pack_nosimd(buf,sp1,sp2,offset_base,a);
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

  //  static inline void unpack_xp(vcomplexD_t &sp1, vcomplexD_t &sp2, vcomplexD_t *rbuf, const int offset_base, const int a=0){
  //    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  //  }
  //  static inline void unpack_xm(vcomplexD_t &sp1, vcomplexD_t &sp2, vcomplexD_t *rbuf, const int offset_base, const int a=0){
  //    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  //  }

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
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }

  static inline void unpack_tm(vcomplexD_t &sp1, vcomplexD_t &sp2, const vcomplexD_t *rbuf, const int offset_base, const int a=0){
    unpack_nosimd(sp1,sp2,rbuf,offset_base,a);
  }


  // boundary condition: mulpliy a real number
  static inline void apply_boundary(vcomplexD_t &val, const vrealD_t &boundary_value){
    val.vec = _mm512_mul_pd(val.vec, boundary_value.vec);
  }

  static inline void set_boundary_xp(vrealD_t &boundary_value, const double rval){
    boundary_value.vec = _mm512_set_pd(rval,rval,1,1,1,1,1,1);
  }

  static inline void set_boundary_xm(vrealD_t &boundary_value, const double rval){
    boundary_value.vec = _mm512_set_pd(1,1,1,1,1,1,rval,rval);
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
    bcast(boundary_value,rval);
  }

  static inline void set_boundary_tm(vrealD_t &boundary_value, const double rval){
    bcast(boundary_value,rval);
  }


} // anonymous namespace

#endif
