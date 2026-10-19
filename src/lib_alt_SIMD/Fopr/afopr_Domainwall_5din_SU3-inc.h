#ifndef AFOPR_DOMAINWALL_5DIN_SIMD_SU3_INC_H
#define AFOPR_DOMAINWALL_5DIN_SIMD_SU3_INC_H
/*!

        @file    $Id: afopr_Domainwall_5din_SU3-inc.h #$

        @brief   SU(3) specials of the Domainwall_5din Dirac operator

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */

//====================================================================

#define SPINPROJ_USE_CHIRAL_REP
#include "Wilson_inc/spinproj_macro-inc.h"
#include "Wilson_inc/matvec_su3_macro-inc.h"
#undef SPINPROJ_USE_CHIRAL_REP

//====================================================================

// macros
#define GET_GAUGE(u, p, s0, mu )   \
  sss=s0*36+9*mu;                  \
  vcomplex_t u ## 00=(p)[sss];     \
  vcomplex_t u ## 01=(p)[sss+1];   \
  vcomplex_t u ## 02=(p)[sss+2];   \
  vcomplex_t u ## 10=(p)[sss+3];   \
  vcomplex_t u ## 11=(p)[sss+4];   \
  vcomplex_t u ## 12=(p)[sss+5];   \
  vcomplex_t u ## 20=(p)[sss+6];   \
  vcomplex_t u ## 21=(p)[sss+7];   \
  vcomplex_t u ## 22=(p)[sss+8];

// dir must be xp, xm, yp, ym,...
#define SHIFT_AND_MERGE_GAUGE(dir, u, un)      \
  shift_and_merge_ ## dir(u ## 00, un ## 00) ; \
  shift_and_merge_ ## dir(u ## 01, un ## 01) ; \
  shift_and_merge_ ## dir(u ## 02, un ## 02) ; \
  shift_and_merge_ ## dir(u ## 10, un ## 10) ; \
  shift_and_merge_ ## dir(u ## 11, un ## 11) ; \
  shift_and_merge_ ## dir(u ## 12, un ## 12) ; \
  shift_and_merge_ ## dir(u ## 20, un ## 20) ; \
  shift_and_merge_ ## dir(u ## 21, un ## 21) ; \
  shift_and_merge_ ## dir(u ## 22, un ## 22) ;

// dir must be xp, xm, yp, ym,...
#define SHIFT_AND_ZERO_GAUGE(dir, u) {          \
    vcomplex_t vzero;                           \
    set_zero(vzero);                            \
    shift_and_merge_ ## dir(u ## 00, vzero) ;   \
    shift_and_merge_ ## dir(u ## 01, vzero) ;   \
    shift_and_merge_ ## dir(u ## 02, vzero) ;   \
    shift_and_merge_ ## dir(u ## 10, vzero) ;   \
    shift_and_merge_ ## dir(u ## 11, vzero) ;   \
    shift_and_merge_ ## dir(u ## 12, vzero) ;   \
    shift_and_merge_ ## dir(u ## 20, vzero) ;   \
    shift_and_merge_ ## dir(u ## 21, vzero) ;   \
    shift_and_merge_ ## dir(u ## 22, vzero) ;   \
  }

#define GET_SPINOR( p, s0 )        \
  sss=s0*12;                       \
  vcomplex_t wc0_sp1=(p)[sss];     \
  vcomplex_t wc1_sp1=(p)[sss+1];   \
  vcomplex_t wc2_sp1=(p)[sss+2];   \
  vcomplex_t wc0_sp2=(p)[sss+3];   \
  vcomplex_t wc1_sp2=(p)[sss+4];   \
  vcomplex_t wc2_sp2=(p)[sss+5];   \
  vcomplex_t wc0_sp3=(p)[sss+6];   \
  vcomplex_t wc1_sp3=(p)[sss+7];   \
  vcomplex_t wc2_sp3=(p)[sss+8];   \
  vcomplex_t wc0_sp4=(p)[sss+9];   \
  vcomplex_t wc1_sp4=(p)[sss+10];  \
  vcomplex_t wc2_sp4=(p)[sss+11];

#define DEFINE_SPINOR \
  vcomplex_t wc0_sp1; \
  vcomplex_t wc1_sp1; \
  vcomplex_t wc2_sp1; \
  vcomplex_t wc0_sp2; \
  vcomplex_t wc1_sp2; \
  vcomplex_t wc2_sp2; \
  vcomplex_t wc0_sp3; \
  vcomplex_t wc1_sp3; \
  vcomplex_t wc2_sp3; \
  vcomplex_t wc0_sp4; \
  vcomplex_t wc1_sp4; \
  vcomplex_t wc2_sp4;



#define SET_SPINOR( p, s0 )\
  {                        \
    int sss=s0*12;         \
    (p)[sss]   =wc0_sp1;   \
    (p)[sss+1] =wc1_sp1;   \
    (p)[sss+2] =wc2_sp1;   \
    (p)[sss+3] =wc0_sp2;   \
    (p)[sss+4] =wc1_sp2;   \
    (p)[sss+5] =wc2_sp2;   \
    (p)[sss+6] =wc0_sp3;   \
    (p)[sss+7] =wc1_sp3;   \
    (p)[sss+8] =wc2_sp3;   \
    (p)[sss+9] =wc0_sp4;   \
    (p)[sss+10]=wc1_sp4;   \
    (p)[sss+11]=wc2_sp4;   \
  }

namespace {
  using real_t=typename complex_t::value_type;
  static inline void site_gm5_Dirac(vcomplex_t *out, const vcomplex_t *in){

    const vcomplex_t sp3_c0=in[6];
    const vcomplex_t sp3_c1=in[7];
    const vcomplex_t sp3_c2=in[8];
    const vcomplex_t sp4_c0=in[9];
    const vcomplex_t sp4_c1=in[10];
    const vcomplex_t sp4_c2=in[11];
    const vcomplex_t sp1_c0=in[0];
    const vcomplex_t sp1_c1=in[1];
    const vcomplex_t sp1_c2=in[2];
    const vcomplex_t sp2_c0=in[3];
    const vcomplex_t sp2_c1=in[4];
    const vcomplex_t sp2_c2=in[5];

    out[0] = sp3_c0;
    out[1] = sp3_c1;
    out[2] = sp3_c2;
    out[3] = sp4_c0;
    out[4] = sp4_c1;
    out[5] = sp4_c2;
    out[6] = sp1_c0;
    out[7] = sp1_c1;
    out[8] = sp1_c2;
    out[9] = sp2_c0;
    out[10] = sp2_c1;
    out[11] = sp2_c2;
    // streaming store version
    //    sstore(&out[0], sp3_c0);
    //    sstore(&out[1], sp3_c1);
    //    sstore(&out[2], sp3_c2);
    //    sstore(&out[3], sp4_c0);
    //    sstore(&out[4], sp4_c1);
    //    sstore(&out[5], sp4_c2);
    //    sstore(&out[6], sp1_c0);
    //    sstore(&out[7], sp1_c1);
    //    sstore(&out[8], sp1_c2);
    //    sstore(&out[9], sp2_c0);
    //    sstore(&out[10], sp2_c1);
    //    sstore(&out[11], sp2_c2);
  }

  static inline void site_gm5_Chiral(vcomplex_t *out, const vcomplex_t *in){
    vcomplex_t vzero;
    set_zero(vzero);

    const vcomplex_t sp1_c0= in[0];
    const vcomplex_t sp1_c1= in[1];
    const vcomplex_t sp1_c2= in[2];
    const vcomplex_t sp2_c0= in[3];
    const vcomplex_t sp2_c1= in[4];
    const vcomplex_t sp2_c2= in[5];
    const vcomplex_t sp3_c0= sub(vzero, in[6]);
    const vcomplex_t sp3_c1= sub(vzero, in[7]);
    const vcomplex_t sp3_c2= sub(vzero, in[8]);
    const vcomplex_t sp4_c0= sub(vzero, in[9]);
    const vcomplex_t sp4_c1= sub(vzero, in[10]);
    const vcomplex_t sp4_c2= sub(vzero, in[11]);

    out[0] = sp1_c0;
    out[1] = sp1_c1;
    out[2] = sp1_c2;
    out[3] = sp2_c0;
    out[4] = sp2_c1;
    out[5] = sp2_c2;
    out[6] = sp3_c0;
    out[7] = sp3_c1;
    out[8] = sp3_c2;
    out[9] = sp4_c0;
    out[10] = sp4_c1;
    out[11] = sp4_c2;
    // streaming store version
    //    sstore(&out[0], sp1_c0);
    //    sstore(&out[1], sp1_c1);
    //    sstore(&out[2], sp1_c2);
    //    sstore(&out[3], sp2_c0);
    //    sstore(&out[4], sp2_c1);
    //    sstore(&out[5], sp2_c2);
    //    sstore(&out[6], sp3_c0);
    //    sstore(&out[7], sp3_c1);
    //    sstore(&out[8], sp3_c2);
    //    sstore(&out[9], sp4_c0);
    //    sstore(&out[10], sp4_c1);
    //    sstore(&out[11], sp4_c2);
  }


  // this function is specific to the chiral rep.
  static inline void prewilson_at_site(vcomplex_t *pv, const vcomplex_t *pw,
                                       const int Ns,
                                       const real_t *c_array, // c_i := c_i/b_i
                                       const real_t m){
    constexpr int Ndf=12;
    assert(Ns % 2 == 0);


    { // chirality +
      const vcomplex_t *in=pw;
      vcomplex_t tmp_sp1_c0=in[0];
      vcomplex_t tmp_sp1_c1=in[1];
      vcomplex_t tmp_sp1_c2=in[2];
      vcomplex_t tmp_sp2_c0=in[3];
      vcomplex_t tmp_sp2_c1=in[4];
      vcomplex_t tmp_sp2_c2=in[5];
      vcomplex_t sp1_c0, sp1_c1, sp1_c2;
      vcomplex_t sp2_c0, sp2_c1, sp2_c2;
      {// to s5=0
        vcomplex_t *out=pv;
        out[0] = tmp_sp1_c0;
        out[1] = tmp_sp1_c1;
        out[2] = tmp_sp1_c2;
        out[3] = tmp_sp2_c0;
        out[4] = tmp_sp2_c1;
        out[5] = tmp_sp2_c2;
      }
      for(int s5=1; s5<Ns-1; s5+=2){ // from s5=1,2,...,Ns-1
        in += Ndf;
        vreal_t vc;

        // to s5
        vcomplex_t *out=pv+ Ndf*s5;
        bcast(vc, c_array[s5]);
        sp1_c0=in[0];
        sp1_c1=in[1];
        sp1_c2=in[2];
        sp2_c0=in[3];
        sp2_c1=in[4];
        sp2_c2=in[5];
        aypx(tmp_sp1_c0, vc, sp1_c0);
        aypx(tmp_sp1_c1, vc, sp1_c1);
        aypx(tmp_sp1_c2, vc, sp1_c2);
        aypx(tmp_sp2_c0, vc, sp2_c0);
        aypx(tmp_sp2_c1, vc, sp2_c1);
        aypx(tmp_sp2_c2, vc, sp2_c2);
        out[0] = tmp_sp1_c0;
        out[1] = tmp_sp1_c1;
        out[2] = tmp_sp1_c2;
        out[3] = tmp_sp2_c0;
        out[4] = tmp_sp2_c1;
        out[5] = tmp_sp2_c2;

        // to s5+1
        in += Ndf;
        out += Ndf;

        bcast(vc, c_array[s5+1]);
        tmp_sp1_c0=in[0];
        tmp_sp1_c1=in[1];
        tmp_sp1_c2=in[2];
        tmp_sp2_c0=in[3];
        tmp_sp2_c1=in[4];
        tmp_sp2_c2=in[5];
        aypx(sp1_c0, vc, tmp_sp1_c0);
        aypx(sp1_c1, vc, tmp_sp1_c1);
        aypx(sp1_c2, vc, tmp_sp1_c2);
        aypx(sp2_c0, vc, tmp_sp2_c0);
        aypx(sp2_c1, vc, tmp_sp2_c1);
        aypx(sp2_c2, vc, tmp_sp2_c2);
        out[0] = sp1_c0;
        out[1] = sp1_c1;
        out[2] = sp1_c2;
        out[3] = sp2_c0;
        out[4] = sp2_c1;
        out[5] = sp2_c2;
      } // s5

      {
        in += Ndf;
        vreal_t vc;

        // to Ns-1
        vcomplex_t *out=pv+ Ndf*(Ns-1);
        bcast(vc, c_array[Ns-1]);
        sp1_c0=in[0];
        sp1_c1=in[1];
        sp1_c2=in[2];
        sp2_c0=in[3];
        sp2_c1=in[4];
        sp2_c2=in[5];
        aypx(tmp_sp1_c0, vc, sp1_c0);
        aypx(tmp_sp1_c1, vc, sp1_c1);
        aypx(tmp_sp1_c2, vc, sp1_c2);
        aypx(tmp_sp2_c0, vc, sp2_c0);
        aypx(tmp_sp2_c1, vc, sp2_c1);
        aypx(tmp_sp2_c2, vc, sp2_c2);
        out[0] = tmp_sp1_c0;
        out[1] = tmp_sp1_c1;
        out[2] = tmp_sp1_c2;
        out[3] = tmp_sp2_c0;
        out[4] = tmp_sp2_c1;
        out[5] = tmp_sp2_c2;

        // to s5=0
        out=pv;
        bcast(vc, -m*c_array[0]);
        axpy(out[0], vc, sp1_c0);
        axpy(out[1], vc, sp1_c1);
        axpy(out[2], vc, sp1_c2);
        axpy(out[3], vc, sp2_c0);
        axpy(out[4], vc, sp2_c1);
        axpy(out[5], vc, sp2_c2);
      }
    } // chirality+


    {// chirality-
      const vcomplex_t *in=pw;
      vcomplex_t tmp_sp3_c0=in[6];
      vcomplex_t tmp_sp3_c1=in[7];
      vcomplex_t tmp_sp3_c2=in[8];
      vcomplex_t tmp_sp4_c0=in[9];
      vcomplex_t tmp_sp4_c1=in[10];
      vcomplex_t tmp_sp4_c2=in[11];
      vcomplex_t sp3_c0, sp3_c1, sp3_c2;
      vcomplex_t sp4_c0, sp4_c1, sp4_c2;

      // from s5=0
      { // to s5=Ns-1
        vreal_t vc;
        bcast(vc, -m*c_array[Ns-1]);
        vcomplex_t *out=pv+Ndf*(Ns-1);
        out[6]=mult(vc, tmp_sp3_c0);
        out[7]=mult(vc, tmp_sp3_c1);
        out[8]=mult(vc, tmp_sp3_c2);
        out[9]=mult(vc, tmp_sp4_c0);
        out[10]=mult(vc, tmp_sp4_c1);
        out[11]=mult(vc, tmp_sp4_c2);
      }
      for(int s5=1; s5<Ns-1; s5+=2){ // from s5=1,2,...,Ns-1
        in += Ndf;
        vreal_t vc;

        // to s5-1
        vcomplex_t *out=pv+ Ndf*(s5-1);
        bcast(vc,  c_array[s5-1]);
        sp3_c0=in[6];
        sp3_c1=in[7];
        sp3_c2=in[8];
        sp4_c0=in[9];
        sp4_c1=in[10];
        sp4_c2=in[11];
        axpy(tmp_sp3_c0, vc, sp3_c0);
        axpy(tmp_sp3_c1, vc, sp3_c1);
        axpy(tmp_sp3_c2, vc, sp3_c2);
        axpy(tmp_sp4_c0, vc, sp4_c0);
        axpy(tmp_sp4_c1, vc, sp4_c1);
        axpy(tmp_sp4_c2, vc, sp4_c2);
        out[6]=tmp_sp3_c0;
        out[7]=tmp_sp3_c1;
        out[8]=tmp_sp3_c2;
        out[9]=tmp_sp4_c0;
        out[10]=tmp_sp4_c1;
        out[11]=tmp_sp4_c2;

        // to s5
        in += Ndf;
        out += Ndf;

        bcast(vc,  c_array[s5]);
        tmp_sp3_c0=in[6];
        tmp_sp3_c1=in[7];
        tmp_sp3_c2=in[8];
        tmp_sp4_c0=in[9];
        tmp_sp4_c1=in[10];
        tmp_sp4_c2=in[11];
        axpy(sp3_c0, vc, tmp_sp3_c0);
        axpy(sp3_c1, vc, tmp_sp3_c1);
        axpy(sp3_c2, vc, tmp_sp3_c2);
        axpy(sp4_c0, vc, tmp_sp4_c0);
        axpy(sp4_c1, vc, tmp_sp4_c1);
        axpy(sp4_c2, vc, tmp_sp4_c2);
        out[6]=sp3_c0;
        out[7]=sp3_c1;
        out[8]=sp3_c2;
        out[9]=sp4_c0;
        out[10]=sp4_c1;
        out[11]=sp4_c2;
      } //s5
      {
        in += Ndf;
        vreal_t vc;

        // to s5-2
        vcomplex_t *out=pv+ Ndf*(Ns-2);
        bcast(vc,  c_array[Ns-1]);
        sp3_c0=in[6];
        sp3_c1=in[7];
        sp3_c2=in[8];
        sp4_c0=in[9];
        sp4_c1=in[10];
        sp4_c2=in[11];
        axpy(tmp_sp3_c0, vc, sp3_c0);
        axpy(tmp_sp3_c1, vc, sp3_c1);
        axpy(tmp_sp3_c2, vc, sp3_c2);
        axpy(tmp_sp4_c0, vc, sp4_c0);
        axpy(tmp_sp4_c1, vc, sp4_c1);
        axpy(tmp_sp4_c2, vc, sp4_c2);
        out[6]=tmp_sp3_c0;
        out[7]=tmp_sp3_c1;
        out[8]=tmp_sp3_c2;
        out[9]=tmp_sp4_c0;
        out[10]=tmp_sp4_c1;
        out[11]=tmp_sp4_c2;

        out+= Ndf;

        accum_add(out[6],sp3_c0);
        accum_add(out[7],sp3_c1);
        accum_add(out[8],sp3_c2);
        accum_add(out[9],sp4_c0);
        accum_add(out[10],sp4_c1);
        accum_add(out[11],sp4_c2);
      }

    } // chirality-
  }


  // this function is specific to the chiral rep.
  static inline void prewilson_at_site_hide(vcomplex_t *pv, const vcomplex_t *pw,
                                            const int Ns,
                                            const real_t *c_array, // c_i := c_i/b_i
                                            const real_t m){
    constexpr int Ndf=12;
    const vcomplex_t *in=pw;

    vcomplex_t tmp_sp1_c0=in[0];
    vcomplex_t tmp_sp1_c1=in[1];
    vcomplex_t tmp_sp1_c2=in[2];
    vcomplex_t tmp_sp2_c0=in[3];
    vcomplex_t tmp_sp2_c1=in[4];
    vcomplex_t tmp_sp2_c2=in[5];
    vcomplex_t tmp_sp3_c0=in[6];
    vcomplex_t tmp_sp3_c1=in[7];
    vcomplex_t tmp_sp3_c2=in[8];
    vcomplex_t tmp_sp4_c0=in[9];
    vcomplex_t tmp_sp4_c1=in[10];
    vcomplex_t tmp_sp4_c2=in[11];

    {    // from s5=0
      { // to s5=Ns-1
        vreal_t vc;
        bcast(vc, -m*c_array[Ns-1]);
        vcomplex_t *out=pv+Ndf*(Ns-1);
        out[6]=mult(vc, tmp_sp3_c0);
        out[7]=mult(vc, tmp_sp3_c1);
        out[8]=mult(vc, tmp_sp3_c2);
        out[9]=mult(vc, tmp_sp4_c0);
        out[10]=mult(vc, tmp_sp4_c1);
        out[11]=mult(vc, tmp_sp4_c2);
      }
      {// to s5=0
        vcomplex_t *out=pv;
        out[0] = tmp_sp1_c0;
        out[1] = tmp_sp1_c1;
        out[2] = tmp_sp1_c2;
        out[3] = tmp_sp2_c0;
        out[4] = tmp_sp2_c1;
        out[5] = tmp_sp2_c2;
      }
    }
    for(int s5=1; s5<Ns; ++s5){ // from s5=1,2,...,Ns-1
      in += Ndf;
      {
        vreal_t vc;
        bcast(vc, c_array[s5]);
        const vcomplex_t sp1_c0=in[0];
        const vcomplex_t sp1_c1=in[1];
        const vcomplex_t sp1_c2=in[2];
        const vcomplex_t sp2_c0=in[3];
        const vcomplex_t sp2_c1=in[4];
        const vcomplex_t sp2_c2=in[5];

        // to s5
        vcomplex_t *out=pv+ Ndf*s5;
        aypx(tmp_sp1_c0, vc, sp1_c0);
        aypx(tmp_sp1_c1, vc, sp1_c1);
        aypx(tmp_sp1_c2, vc, sp1_c2);
        aypx(tmp_sp2_c0, vc, sp2_c0);
        aypx(tmp_sp2_c1, vc, sp2_c1);
        aypx(tmp_sp2_c2, vc, sp2_c2);
        out[0] = tmp_sp1_c0;
        out[1] = tmp_sp1_c1;
        out[2] = tmp_sp1_c2;
        out[3] = tmp_sp2_c0;
        out[4] = tmp_sp2_c1;
        out[5] = tmp_sp2_c2;

        // to s5+1
        tmp_sp1_c0 = sp1_c0;
        tmp_sp1_c1 = sp1_c1;
        tmp_sp1_c2 = sp1_c2;
        tmp_sp2_c0 = sp2_c0;
        tmp_sp2_c1 = sp2_c1;
        tmp_sp2_c2 = sp2_c2;
      }
      {
        vreal_t vc;
        bcast(vc,  c_array[s5-1]);
        const vcomplex_t sp3_c0=in[6];
        const vcomplex_t sp3_c1=in[7];
        const vcomplex_t sp3_c2=in[8];
        const vcomplex_t sp4_c0=in[9];
        const vcomplex_t sp4_c1=in[10];
        const vcomplex_t sp4_c2=in[11];

        // to s5-1
        vcomplex_t *out=pv+ Ndf*(s5-1);
        axpy(tmp_sp3_c0, vc, sp3_c0);
        axpy(tmp_sp3_c1, vc, sp3_c1);
        axpy(tmp_sp3_c2, vc, sp3_c2);
        axpy(tmp_sp4_c0, vc, sp4_c0);
        axpy(tmp_sp4_c1, vc, sp4_c1);
        axpy(tmp_sp4_c2, vc, sp4_c2);
        out[6]=tmp_sp3_c0;
        out[7]=tmp_sp3_c1;
        out[8]=tmp_sp3_c2;
        out[9]=tmp_sp4_c0;
        out[10]=tmp_sp4_c1;
        out[11]=tmp_sp4_c2;

        // to s5
        tmp_sp3_c0 = sp3_c0;
        tmp_sp3_c1 = sp3_c1;
        tmp_sp3_c2 = sp3_c2;
        tmp_sp4_c0 = sp4_c0;
        tmp_sp4_c1 = sp4_c1;
        tmp_sp4_c2 = sp4_c2;
      }
    } //s5

    { // to s5=Ns-1
      vcomplex_t *out=pv+ Ndf*(Ns-1);
      accum_add(out[6],tmp_sp3_c0);
      accum_add(out[7],tmp_sp3_c1);
      accum_add(out[8],tmp_sp3_c2);
      accum_add(out[9],tmp_sp4_c0);
      accum_add(out[10],tmp_sp4_c1);
      accum_add(out[11],tmp_sp4_c2);
    }
    {// to s5=0
      vcomplex_t *out=pv;
      vreal_t vc;
      bcast(vc, -m*c_array[0]);
      axpy(out[0], vc, tmp_sp1_c0);
      axpy(out[1], vc, tmp_sp1_c1);
      axpy(out[2], vc, tmp_sp1_c2);
      axpy(out[3], vc, tmp_sp2_c0);
      axpy(out[4], vc, tmp_sp2_c1);
      axpy(out[5], vc, tmp_sp2_c2);
    }
  }

  // this function is specific to the chiral rep.
  static inline void prewilson_dag_at_site(vcomplex_t *pv, const vcomplex_t *pw,
                                           const int Ns,
                                           const real_t *cm_array, // c[s5-1]/b[s5]
                                           const real_t *cp_array, // c[s5+1]/b[s5]
                                           const real_t m){
    constexpr int Ndf=12;

    { // chirality +
      const vcomplex_t *in=pw;

      vcomplex_t tmp_sp1_c0=in[0];
      vcomplex_t tmp_sp1_c1=in[1];
      vcomplex_t tmp_sp1_c2=in[2];
      vcomplex_t tmp_sp2_c0=in[3];
      vcomplex_t tmp_sp2_c1=in[4];
      vcomplex_t tmp_sp2_c2=in[5];
      vcomplex_t sp1_c0, sp1_c1, sp1_c2;
      vcomplex_t sp2_c0, sp2_c1, sp2_c2;
      {    // from s5=0
        vreal_t vc;
        bcast(vc, -m*cp_array[Ns-1]); // c[0]/b[Ns-1]

        // to s5=Ns-1
        vcomplex_t *out=pv+Ndf*(Ns-1);
        out[0]=mult(vc, tmp_sp1_c0);
        out[1]=mult(vc, tmp_sp1_c1);
        out[2]=mult(vc, tmp_sp1_c2);
        out[3]=mult(vc, tmp_sp2_c0);
        out[4]=mult(vc, tmp_sp2_c1);
        out[5]=mult(vc, tmp_sp2_c2);
      }
      vcomplex_t *out=pv;
      for(int s5=1; s5<Ns-1; s5+=2){ // from s5=1,2,...,Ns-1
        in += Ndf;
        vreal_t vc;

        // to s5-1
        vcomplex_t *out=pv+ Ndf*(s5-1);
        bcast(vc,  cp_array[s5-1]);  // c[s5]/b[s5-1]

        sp1_c0=in[0];
        sp1_c1=in[1];
        sp1_c2=in[2];
        sp2_c0=in[3];
        sp2_c1=in[4];
        sp2_c2=in[5];

        axpy(tmp_sp1_c0, vc, sp1_c0);
        axpy(tmp_sp1_c1, vc, sp1_c1);
        axpy(tmp_sp1_c2, vc, sp1_c2);
        axpy(tmp_sp2_c0, vc, sp2_c0);
        axpy(tmp_sp2_c1, vc, sp2_c1);
        axpy(tmp_sp2_c2, vc, sp2_c2);

        out[0]=tmp_sp1_c0;
        out[1]=tmp_sp1_c1;
        out[2]=tmp_sp1_c2;
        out[3]=tmp_sp2_c0;
        out[4]=tmp_sp2_c1;
        out[5]=tmp_sp2_c2;

        // to s5
        in += Ndf;
        out += Ndf;

        bcast(vc,  cp_array[s5]); // c[s5+1]/b[s5]

        tmp_sp1_c0=in[0];  // read from s5+1
        tmp_sp1_c1=in[1];
        tmp_sp1_c2=in[2];
        tmp_sp2_c0=in[3];
        tmp_sp2_c1=in[4];
        tmp_sp2_c2=in[5];

        axpy(sp1_c0, vc, tmp_sp1_c0);
        axpy(sp1_c1, vc, tmp_sp1_c1);
        axpy(sp1_c2, vc, tmp_sp1_c2);
        axpy(sp2_c0, vc, tmp_sp2_c0);
        axpy(sp2_c1, vc, tmp_sp2_c1);
        axpy(sp2_c2, vc, tmp_sp2_c2);

        out[0]=sp1_c0;
        out[1]=sp1_c1;
        out[2]=sp1_c2;
        out[3]=sp2_c0;
        out[4]=sp2_c1;
        out[5]=sp2_c2;

      } //s5

      {
        in += Ndf;

        // to s5=Ns-2
        vcomplex_t *out=pv+ Ndf*(Ns-2);
        vreal_t vc;
        bcast(vc,  cp_array[Ns-2]); // c[Ns-1]/b[Ns-2]
        sp1_c0=in[0];
        sp1_c1=in[1];
        sp1_c2=in[2];
        sp2_c0=in[3];
        sp2_c1=in[4];
        sp2_c2=in[5];

        axpy(tmp_sp1_c0, vc, sp1_c0);
        axpy(tmp_sp1_c1, vc, sp1_c1);
        axpy(tmp_sp1_c2, vc, sp1_c2);
        axpy(tmp_sp2_c0, vc, sp2_c0);
        axpy(tmp_sp2_c1, vc, sp2_c1);
        axpy(tmp_sp2_c2, vc, sp2_c2);

        out[0]=tmp_sp1_c0;
        out[1]=tmp_sp1_c1;
        out[2]=tmp_sp1_c2;
        out[3]=tmp_sp2_c0;
        out[4]=tmp_sp2_c1;
        out[5]=tmp_sp2_c2;

        // to s5=Ns-1
        out += Ndf;
        accum_add(out[0], sp1_c0);
        accum_add(out[1], sp1_c1);
        accum_add(out[2], sp1_c2);
        accum_add(out[3], sp2_c0);
        accum_add(out[4], sp2_c1);
        accum_add(out[5], sp2_c2);
      }
    } // chirality +


    { // chirality -

      const vcomplex_t *in=pw;

      vcomplex_t tmp_sp3_c0=in[6];
      vcomplex_t tmp_sp3_c1=in[7];
      vcomplex_t tmp_sp3_c2=in[8];
      vcomplex_t tmp_sp4_c0=in[9];
      vcomplex_t tmp_sp4_c1=in[10];
      vcomplex_t tmp_sp4_c2=in[11];
      vcomplex_t sp3_c0, sp3_c1, sp3_c2;
      vcomplex_t sp4_c0, sp4_c1, sp4_c2;

      {// to s5=0
        vcomplex_t *out=pv;
        out[6]  = tmp_sp3_c0;
        out[7]  = tmp_sp3_c1;
        out[8]  = tmp_sp3_c2;
        out[9]  = tmp_sp4_c0;
        out[10] = tmp_sp4_c1;
        out[11] = tmp_sp4_c2;
      }
      {
        // to Ns-1
        vcomplex_t vc;
        bcast(vc, cm_array[Ns-1]); // c[0]/b[Ns-1]
        vcomplex_t *out=pv+(Ns-1)*Ndf;
        out[6] = mult(vc, tmp_sp3_c0);
        out[7] = mult(vc, tmp_sp3_c1);
        out[8] = mult(vc, tmp_sp3_c2);
        out[9] = mult(vc, tmp_sp4_c0);
        out[10] = mult(vc, tmp_sp4_c1);
        out[11] = mult(vc, tmp_sp4_c2);

      }
      for(int s5=1; s5<Ns-1; s5+=2){ // from s5=1,2,...,Ns-1
        in += Ndf;
        vcomplex_t *out=pv+ Ndf*s5;
        vreal_t vc;

        // to s5
        bcast(vc,  cm_array[s5]); // c[s5-1]/b[s5]
        sp3_c0=in[6];
        sp3_c1=in[7];
        sp3_c2=in[8];
        sp4_c0=in[9];
        sp4_c1=in[10];
        sp4_c2=in[11];

        aypx(tmp_sp3_c0, vc, sp3_c0);
        aypx(tmp_sp3_c1, vc, sp3_c1);
        aypx(tmp_sp3_c2, vc, sp3_c2);
        aypx(tmp_sp4_c0, vc, sp4_c0);
        aypx(tmp_sp4_c1, vc, sp4_c1);
        aypx(tmp_sp4_c2, vc, sp4_c2);
        out[6] = tmp_sp3_c0;
        out[7] = tmp_sp3_c1;
        out[8] = tmp_sp3_c2;
        out[9] = tmp_sp4_c0;
        out[10] = tmp_sp4_c1;
        out[11] = tmp_sp4_c2;

        // to s5+1
        in+=Ndf;
        out+=Ndf;

        bcast(vc,  cm_array[s5+1]); // c[s5]/b[s5+1]
        tmp_sp3_c0=in[6];
        tmp_sp3_c1=in[7];
        tmp_sp3_c2=in[8];
        tmp_sp4_c0=in[9];
        tmp_sp4_c1=in[10];
        tmp_sp4_c2=in[11];

        aypx(sp3_c0, vc, tmp_sp3_c0);
        aypx(sp3_c1, vc, tmp_sp3_c1);
        aypx(sp3_c2, vc, tmp_sp3_c2);
        aypx(sp4_c0, vc, tmp_sp4_c0);
        aypx(sp4_c1, vc, tmp_sp4_c1);
        aypx(sp4_c2, vc, tmp_sp4_c2);
        out[6] = sp3_c0;
        out[7] = sp3_c1;
        out[8] = sp3_c2;
        out[9] = sp4_c0;
        out[10] = sp4_c1;
        out[11] = sp4_c2;

      } //s5

      {// to s5=Ns-1
        in += Ndf;

        vreal_t vc;
        vcomplex_t *out=pv+Ndf*(Ns-1);

        bcast(vc,  cm_array[Ns-1]); // c[Ns-2]/b[Ns-1]
        sp3_c0=in[6];
        sp3_c1=in[7];
        sp3_c2=in[8];
        sp4_c0=in[9];
        sp4_c1=in[10];
        sp4_c2=in[11];

        aypx(tmp_sp3_c0, vc, sp3_c0);
        aypx(tmp_sp3_c1, vc, sp3_c1);
        aypx(tmp_sp3_c2, vc, sp3_c2);
        aypx(tmp_sp4_c0, vc, sp4_c0);
        aypx(tmp_sp4_c1, vc, sp4_c1);
        aypx(tmp_sp4_c2, vc, sp4_c2);
        out[6] = tmp_sp3_c0;
        out[7] = tmp_sp3_c1;
        out[8] = tmp_sp3_c2;
        out[9] = tmp_sp4_c0;
        out[10] = tmp_sp4_c1;
        out[11] = tmp_sp4_c2;

        out += Ndf;

        // s5=0
        out=pv;
        bcast(vc,  -m*cm_array[0]);  // c[Ns-1]/b[0]
        axpy(out[6],  vc, sp3_c0);
        axpy(out[7],  vc, sp3_c1);
        axpy(out[8],  vc, sp3_c2);
        axpy(out[9],  vc, sp4_c0);
        axpy(out[10], vc, sp4_c1);
        axpy(out[11], vc, sp4_c2);
      }
    } // chirality-
  }
} // anonymous namespace


//====================================================================
template<>
void AFopr_Domainwall_5din<AFIELDspec>::mult_R(AFIELDspec& v, const AFIELDspec& w)
{
  constexpr int Ndf_ferm=12;
  assert(m_Ndf_ferm == Ndf_ferm);
  assert(w.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(v.check_size(2*m_Ndf_ferm5, m_Nst, 1));

  const vcomplex_t *pw= w.ptr_vcomplex();
  vcomplex_t *pv= v.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Rst);

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

  for(int site = is; site < ns; ++site){
    const size_t offset=site*Ndf_ferm5;
    for(int i=0; i<Ns; ++i){
      const vcomplex_t *in=pw + offset + i*Ndf_ferm;
      vcomplex_t *out =pv + offset + (Ns-i-1)*Ndf_ferm;
      for(int a=0; a<Ndf_ferm; ++a){ // expect loop unroling by compiler
        out[a]=in[a];
        //vcomplex_t val=in[a];
        //sstore(&out[a], val);
      }
    }
  }
}


//====================================================================
template<>
void AFopr_Domainwall_5din<AFIELDspec>::mult_gm5R(AFIELDspec& v, const AFIELDspec& w)
{
  constexpr int Ndf_ferm=12;
  assert(m_Ndf_ferm == Ndf_ferm);
  assert(w.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(v.check_size(2*m_Ndf_ferm5, m_Nst, 1));

  const vcomplex_t *pw= w.ptr_vcomplex();
  vcomplex_t *pv= v.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Rst);

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

  if(m_repr == "Dirac"){
    for(size_t site = is; site < ns; ++site){
      const size_t offset=site*Ndf_ferm5;
      for(int i=0; i<Ns; ++i){
        const vcomplex_t *in=pw + offset + i*Ndf_ferm;
        vcomplex_t *out =pv + offset + (Ns-i-1)*Ndf_ferm;
        site_gm5_Dirac(out,in);
      }
    }
  } else if(m_repr == "Chiral"){
    for(size_t site = is; site < ns; ++site){
      const size_t offset=site*Ndf_ferm5;
      for(int i=0; i<Ns; ++i){
        const vcomplex_t *in=pw + offset + i*Ndf_ferm;
        vcomplex_t *out =pv + offset + (Ns-i-1)*Ndf_ferm;
        site_gm5_Chiral(out,in);
      }
    }
  } else {
    vout.crucial("%s: unknown rep in mult_gm5R, cannot happen\n", class_name.c_str());
  }
}

//====================================================================
template<>
void AFopr_Domainwall_5din<AFIELDspec>::mult_gm5(AFIELDspec& v, const AFIELDspec& w)
{
  constexpr int Ndf_ferm=12;
  assert(m_Ndf_ferm == Ndf_ferm);
  assert(w.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(v.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  const vcomplex_t *pw= w.ptr_vcomplex();
  vcomplex_t *pv= v.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Rst);
  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

  if(m_repr == "Dirac"){
    for(size_t site = is; site < ns; ++site){
      const size_t offset=site*Ndf_ferm5;
      for(int i=0; i<Ns; ++i){
        const vcomplex_t *in=pw + offset + i*Ndf_ferm;
        vcomplex_t *out =pv + offset + i*Ndf_ferm;
        site_gm5_Dirac(out,in);
      }
    }
  } else if(m_repr == "Chiral"){
    for(size_t site = is; site < ns; ++site){
      const size_t offset=site*Ndf_ferm5;
      for(int i=0; i<Ns; ++i){
        const vcomplex_t *in=pw + offset + i*Ndf_ferm;
        vcomplex_t *out =pv + offset + i*Ndf_ferm;
        site_gm5_Chiral(out,in);
      }
    }
  } else {
    vout.crucial("%s: unknown rep in mult_gm5, cannot happen\n", class_name.c_str());
  }
}

//====================================================================
template<>
void AFopr_Domainwall_5din<AFIELDspec>::D(AFIELDspec& v, const AFIELDspec& w)
{
  constexpr int Ndf_ferm=12;
  assert(m_Ndf_ferm == Ndf_ferm);
  assert(w.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(v.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(CommonParameters::Nc() == 3);

  if(debug) {
    vout.general("afopr_Domainwall_5din: D, started [SU(3) version]\n");
    double w2=w.norm2();
    vout.general("    input norm2 = %f\n", w2);
  }

  const vcomplex_t *pw_full= w.ptr_vcomplex();
  vcomplex_t *pv_full= v.ptr_vcomplex();

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

  // reduced, after divided by the simd layout
  const int Rst=m_Rst;
  const int Rx=m_Rx;
  const int Ry=m_Ry;
  const int Rz=m_Rz;
  const int Rt=m_Rt;

  const int bulk_size=m_Rst;
  const int boundary_size=m_Rst_boundary;

  int sss; // working index for macro

  // initial barrier
#pragma omp barrier
// thread assignments
  int ith, nth, is_bulk, ns_bulk;

#ifdef USE_MPI
  mpi_bools bools=m_bools;
  if(bools.has_mpi_t || bools.has_mpi_z || bools.has_mpi_y){
    set_threadtask_without_master(ith, nth, is_bulk, ns_bulk, bulk_size);
  }else{
    set_threadtask(ith, nth, is_bulk, ns_bulk, bulk_size);
  }
  int is_boundary, ns_boundary;
  int isx,nsx,isy,nsy,isz,nsz,ist,nst;
  set_threadtask_without_master(ith, nth, isy, nsy, bulk_size/Ry);
  set_threadtask_without_master(ith, nth, isz, nsz, bulk_size/Rz);
  set_threadtask_without_master(ith, nth, ist, nst, bulk_size/Rt);
  set_threadtask(ith, nth, is_boundary, ns_boundary, boundary_size);
#else
  set_threadtask(ith, nth, is_bulk, ns_bulk, bulk_size);
#endif // USE_MPI

  // work vector
  vcomplex_t *pv_site  = &m_v_at_site_array[ith][0];
  vcomplex_t *pw_site  = &m_w_at_site_array[ith][0];
  vcomplex_t *pw0_site = &m_w0_at_site_array[ith][0];

  // coefficients
  const real_t *c_array = &m_c_over_b_array[ith][0];
  const real_t *c_diag_array = &m_c_over_b_diag_array[ith][0];
  const real_t *b_array = &m_b_array[ith][0];
  const real_t mq = m_mq;
  const real_t M0 = m_M0;
  const vcomplex_t *Uptr=m_Uptr;

#ifdef USE_MPI
  // prepare sending buffer
  comm_buffer.start_receiving();

#define PACK_TO_BUFFER

  if(bools.has_mpi_y){
#ifdef RUN_HOP_YP
    for(int i=isy; i<nsy; ++i){
      const int ss = m_site_list_pack_y[i].here;
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_yp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_YM
    for(int i=isy; i<nsy; ++i){
      const int ss = m_site_list_pack_y[i].here + Rx*(Ry-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      GET_GAUGE(u,  Uptr, ss, DirY );
      const vcomplex_t *pw = pw_site;
      vcomplex_t *pv = buf_ym + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_ym-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: y-direction

  if(bools.has_mpi_z){
#ifdef RUN_HOP_ZP
    for(int i=isz; i<nsz; ++i){
      const int ss = m_site_list_pack_z[i].here;
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_zp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_ZM
    for(int i=isz; i<nsz; ++i){
      const int ss=m_site_list_pack_z[i].here+Rx*Ry*(Rz-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      GET_GAUGE(u,  Uptr, ss, DirZ );
      const vcomplex_t *pw = pw_site;
      vcomplex_t *pv = buf_zm + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_zm-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: z-direction

  if(bools.has_mpi_t){
#ifdef RUN_HOP_TP
    for(int i=ist; i<nst; ++i){
      const int ss = m_site_list_pack_t[i].here;
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_tp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_TM
    for(int i=ist; i<nst; ++i){
      const int ss = m_site_list_pack_t[i].here + Rx*Ry*Rz*(Rt-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      GET_GAUGE(u,  Uptr, ss, DirT );
      const vcomplex_t *pw = pw_site;
      vcomplex_t *pv = buf_tm + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_tm-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: t-direction
#undef PACK_TO_BUFFER

#pragma omp barrier
comm_buffer.start_sending();
comm_buffer.test_all();
#endif  // USE_MPI

for(int s0=is_bulk; s0<ns_bulk; ++s0){

    // get coordinates
    const int s=m_site_list_bulk[s0].here;
    const int sn_xp=m_site_list_bulk[s0].neighbor[DirXp];
    const int sn_xm=m_site_list_bulk[s0].neighbor[DirXm];
    const int sn_yp=m_site_list_bulk[s0].neighbor[DirYp];
    const int sn_ym=m_site_list_bulk[s0].neighbor[DirYm];
    const int sn_zp=m_site_list_bulk[s0].neighbor[DirZp];
    const int sn_zm=m_site_list_bulk[s0].neighbor[DirZm];
    const int sn_tp=m_site_list_bulk[s0].neighbor[DirTp];
    const int sn_tm=m_site_list_bulk[s0].neighbor[DirTm];

#ifndef RUN_HOP_XP
    { // initialize work vector to store hoppings
      vcomplex_t vzero;
      set_zero(vzero);
      for(int i=0; i<Ndf_ferm5; ++i){
        pv_site[i]=vzero;
      }
    }
#endif

    const vcomplex_t *in_base= pw_full;
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
#include "Wilson_inc/afopr_Domainwall_5din_SU3_impl1_at_site_bulk-inc.h"
  } // s0 for bulk



 #ifdef USE_MPI
// finalize the communication, all direction
comm_buffer.wait_receiving();
#pragma omp barrier


for(int s0=is_boundary; s0<ns_boundary; ++s0){

  // get coordinates
  const int s=m_site_list_boundary[s0].here;
  const int s_packed_y = m_site_list_boundary[s0].packed[DirY];
  const int s_packed_z = m_site_list_boundary[s0].packed[DirZ];
  const int s_packed_t = m_site_list_boundary[s0].packed[DirT];
  const int sn_xp=m_site_list_boundary[s0].neighbor[DirXp];
  const int sn_xm=m_site_list_boundary[s0].neighbor[DirXm];
  const int sn_yp=m_site_list_boundary[s0].neighbor[DirYp];
  const int sn_ym=m_site_list_boundary[s0].neighbor[DirYm];
  const int sn_zp=m_site_list_boundary[s0].neighbor[DirZp];
  const int sn_zm=m_site_list_boundary[s0].neighbor[DirZm];
  const int sn_tp=m_site_list_boundary[s0].neighbor[DirTp];
  const int sn_tm=m_site_list_boundary[s0].neighbor[DirTm];

#define NEED_TO_READ_BUFFER

  if(sn_yp>-1){
#ifdef RUN_HOP_YP
    vcomplex_t *rbuf = rbuf_yp + s_packed_y * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirY );
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
#include "Wilson_inc/multi_hop_yp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // RUN_HOP_YP
  }
  if(sn_ym>-1){
#ifdef RUN_HOP_YM
    vcomplex_t *rbuf = rbuf_ym + s_packed_y * Ndf_ferm5/2;
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      #include "Wilson_inc/multi_hop_ym-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // N_HOP_YM
  }

  if(sn_zp>-1){
#ifdef RUN_HOP_ZP
    vcomplex_t *rbuf=rbuf_zp + s_packed_z * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirZ );
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
#include "Wilson_inc/multi_hop_zp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // RUN_HOP_ZP
  }
  if(sn_zm>-1){
#ifdef RUN_HOP_ZM
    vcomplex_t *rbuf=rbuf_zm + s_packed_z * Ndf_ferm5/2;
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      #include "Wilson_inc/multi_hop_zm-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // N_HOP_ZM
  }

  if(sn_tp>-1){
#ifdef RUN_HOP_TP
    vcomplex_t *rbuf=rbuf_tp + s_packed_t * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirT );
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
#include "Wilson_inc/multi_hop_tp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // RUN_HOP_TP
  }
  if(sn_tm>-1){
#ifdef RUN_HOP_TM
    vcomplex_t *rbuf=rbuf_tm + s_packed_t * Ndf_ferm5/2;
    vcomplex_t *pv = pv_full + s*Ndf_ferm5;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      #include "Wilson_inc/multi_hop_tm-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
#endif // N_HOP_TM
  }

 } // s0
comm_buffer.wait_sending();

#undef NEED_TO_READ_BUFFER

#endif // USE_MPI 
  // final barrier
#pragma omp barrier

}

template<>
void AFopr_Domainwall_5din<AFIELDspec>::Ddag(AFIELDspec& v, const AFIELDspec& w)
{
  constexpr int Ndf_ferm=12;
  assert(m_Ndf_ferm == Ndf_ferm);
  assert(w.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  assert(v.check_size(2*m_Ndf_ferm5, m_Nst, 1));
  if(debug) {
    vout.general("afopr_Domainwall_5din: Ddag, started [SU(3) version]\n");
  }
#define USE_DAGGER

  const vcomplex_t *pw_full = w.ptr_vcomplex();
  vcomplex_t *pv_full = v.ptr_vcomplex();

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

  // reduced, after divided by the simd layout
  const int Rst=m_Rst;
  const int Rx=m_Rx;
  const int Ry=m_Ry;
  const int Rz=m_Rz;
  const int Rt=m_Rt;

  const int bulk_size=m_Rst;
  const int boundary_size=m_Rst_boundary;
  assert(CommonParameters::Nc() == 3);

  int sss; // working index for macro

  // initial barrier
#pragma omp barrier
// thread assignments
  int ith, nth, is_bulk, ns_bulk;

#ifdef USE_MPI

  mpi_bools bools=m_bools;
  if(bools.has_mpi_t || bools.has_mpi_z || bools.has_mpi_y){
    set_threadtask_without_master(ith, nth, is_bulk, ns_bulk, bulk_size);
  }else{
    set_threadtask(ith, nth, is_bulk, ns_bulk, bulk_size);
  }
  int is_boundary, ns_boundary;
  int isx,nsx,isy,nsy,isz,nsz,ist,nst;
  set_threadtask_without_master(ith, nth, isy, nsy, bulk_size/Ry);
  set_threadtask_without_master(ith, nth, isz, nsz, bulk_size/Rz);
  set_threadtask_without_master(ith, nth, ist, nst, bulk_size/Rt);
  set_threadtask(ith, nth, is_boundary, ns_boundary, boundary_size);

#else
  set_threadtask(ith, nth, is_bulk, ns_bulk, bulk_size);
#endif // USE_MPI

  // work vector
  vcomplex_t *pv_site  = &m_v_at_site_array[ith][0];
  vcomplex_t *pw_site  = &m_w_at_site_array[ith][0];
  vcomplex_t *pw0_site = &m_w0_at_site_array[ith][0];

  // coefficients
  const real_t *cm_array = &m_cm_over_b_array[ith][0];
  const real_t *cp_array = &m_cp_over_b_array[ith][0];
  const real_t *cm_diag_array = &m_cm_over_b_diag_array[ith][0];
  const real_t *cp_diag_array = &m_cp_over_b_diag_array[ith][0];
  const real_t *b_array = &m_b_array[ith][0];
  const real_t mq = m_mq;
  const real_t M0 = m_M0;

  const vcomplex_t *Uptr=m_Uptr;


#ifdef USE_MPI
  // prepare sending buffer
  comm_buffer.start_receiving();

#define PACK_TO_BUFFER
#define SPINPROJ_USE_DAGGER

  if(bools.has_mpi_y){
#ifdef RUN_HOP_YP
    for(int i=isy; i<nsy; ++i){
      const int ss = m_site_list_pack_y[i].here;
      const vcomplex_t *pw = pw_full + ss*Ndf_ferm5;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_yp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_YM
    for(int i=isy; i<nsy; ++i){
      const int ss = m_site_list_pack_y[i].here + Rx*(Ry-1);
      GET_GAUGE(u,  Uptr, ss, DirY );
      const vcomplex_t *pw =pw_full + ss*Ndf_ferm5;
      vcomplex_t *pv = buf_ym + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_ym-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: y-direction

  if(bools.has_mpi_z){
#ifdef RUN_HOP_ZP
    for(int i=isz; i<nsz; ++i){
      const int ss = m_site_list_pack_z[i].here;
      const vcomplex_t *pw = pw_full + ss*Ndf_ferm5;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_zp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_ZM
    for(int i=isz; i<nsz; ++i){
      const int ss = m_site_list_pack_z[i].here + Rx*Ry*(Rz-1);
      GET_GAUGE(u,  Uptr, ss, DirZ );
      const vcomplex_t *pw =pw_full + ss*Ndf_ferm5;
      vcomplex_t *pv = buf_zm + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_zm-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: z-direction

  if(bools.has_mpi_t){
#ifdef RUN_HOP_TP
    for(int i=ist; i<nst; ++i){
      const int ss = m_site_list_pack_t[i].here;
      const vcomplex_t *pw = pw_full + ss*Ndf_ferm5;
      int s_packed=i*Ns;
      for(int s5=0; s5<Ns; ++s5){
        const int s=s5;
#include "Wilson_inc/pack_tp-inc.h"
        s_packed++;
      }
    } //i
#endif
#ifdef RUN_HOP_TM
    for(int i=ist; i<nst; ++i){
      const int ss = m_site_list_pack_t[i].here + Rx*Ry*Rz*(Rt-1);
      GET_GAUGE(u,  Uptr, ss, DirT );
      const vcomplex_t *pw =pw_full + ss*Ndf_ferm5;
      vcomplex_t *pv = buf_tm + i*Ndf_ferm5/2;
      for(int s5=0; s5<Ns; ++s5){
#include "Wilson_inc/multi_hop_tm-inc.h"
        pw+=Ndf_ferm;
        pv+=Ndf_ferm/2;
      }
    } //i
#endif
  } // packing: t-direction

#undef PACK_TO_BUFFER
#undef SPINPROJ_USE_DAGGER

#pragma omp barrier
comm_buffer.start_sending();
comm_buffer.test_all();
#endif  // USE_MPI

  for(int s0=is_bulk; s0<ns_bulk; ++s0){

    // get coordinates
    const int s=m_site_list_bulk[s0].here;
    const int sn_xp=m_site_list_bulk[s0].neighbor[DirXp];
    const int sn_xm=m_site_list_bulk[s0].neighbor[DirXm];
    const int sn_yp=m_site_list_bulk[s0].neighbor[DirYp];
    const int sn_ym=m_site_list_bulk[s0].neighbor[DirYm];
    const int sn_zp=m_site_list_bulk[s0].neighbor[DirZp];
    const int sn_zm=m_site_list_bulk[s0].neighbor[DirZm];
    const int sn_tp=m_site_list_bulk[s0].neighbor[DirTp];
    const int sn_tm=m_site_list_bulk[s0].neighbor[DirTm];

#ifndef RUN_HOP_XP
    { // initialize work vector to store hoppings
      vcomplex_t vzero;
      set_zero(vzero);
      for(int i=0; i<Ndf_ferm5; ++i){
        pv_site[i]=vzero;
      }
    }
#endif
    const vcomplex_t *in_base= pw_full;
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
#include "Wilson_inc/afopr_Domainwall_5din_SU3_impl1_at_site_bulk-inc.h"
#undef USE_DAGGER
  } // s0 for bulk


 #ifdef USE_MPI
// finalize the communication, all direction
comm_buffer.wait_receiving();
#pragma omp barrier

for(int s0=is_boundary; s0<ns_boundary; ++s0){

  // get coordinates
  const int s=m_site_list_boundary[s0].here;
  const int s_packed_y = m_site_list_boundary[s0].packed[DirY];
  const int s_packed_z = m_site_list_boundary[s0].packed[DirZ];
  const int s_packed_t = m_site_list_boundary[s0].packed[DirT];
  const int sn_xp=m_site_list_boundary[s0].neighbor[DirXp];
  const int sn_xm=m_site_list_boundary[s0].neighbor[DirXm];
  const int sn_yp=m_site_list_boundary[s0].neighbor[DirYp];
  const int sn_ym=m_site_list_boundary[s0].neighbor[DirYm];
  const int sn_zp=m_site_list_boundary[s0].neighbor[DirZp];
  const int sn_zm=m_site_list_boundary[s0].neighbor[DirZm];
  const int sn_tp=m_site_list_boundary[s0].neighbor[DirTp];
  const int sn_tm=m_site_list_boundary[s0].neighbor[DirTm];


#define NEED_TO_READ_BUFFER
#define SPINPROJ_USE_DAGGER
#define HOP_DOES_NOT_ACCUM
  if(sn_yp>-1){
#ifdef RUN_HOP_YP
    vcomplex_t *rbuf = rbuf_yp + s_packed_y * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirY );
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      //
#include "Wilson_inc/multi_hop_yp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_YP
  }
  if(sn_ym>-1){
#ifdef RUN_HOP_YM
    vcomplex_t *rbuf = rbuf_ym + s_packed_y * Ndf_ferm5/2;
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }

      #include "Wilson_inc/multi_hop_ym-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_YM
  }

  if(sn_zp>-1){
#ifdef RUN_HOP_ZP
    vcomplex_t *rbuf = rbuf_zp + s_packed_z * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirZ );
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      //
#include "Wilson_inc/multi_hop_zp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_ZP
  }
  if(sn_zm>-1){
#ifdef RUN_HOP_ZM
    vcomplex_t *rbuf = rbuf_zm + s_packed_z * Ndf_ferm5/2;
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }

      #include "Wilson_inc/multi_hop_zm-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_ZM
  }

  if(sn_tp>-1){
#ifdef RUN_HOP_TP
    vcomplex_t *rbuf = rbuf_tp + s_packed_t * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirT );
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }
      //
#include "Wilson_inc/multi_hop_tp-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_TP
  }
  if(sn_tm>-1){
#ifdef RUN_HOP_TM
    vcomplex_t *rbuf = rbuf_tm + s_packed_t * Ndf_ferm5/2;
    vcomplex_t *pv = pv_site;
    const vcomplex_t *pw = pw_full + s*Ndf_ferm5;  // not actually used in the impl1
    for(int s5=0; s5<Ns; ++s5){
      vreal_t vb;
      bcast(vb, -0.5*b_array[s5]);
      for(int i=0; i<Ndf_ferm/2; ++i){
        scal(rbuf[i], vb);
      }

      #include "Wilson_inc/multi_hop_tm-inc.h"
      rbuf += Ndf_ferm/2;
      pv += Ndf_ferm;
      pw += Ndf_ferm;
    }
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
#endif // RUN_HOP_YM
  }


 } // s0

#undef HOP_DOES_NOT_ACCUM
#undef NEED_TO_READ_BUFFER
#undef SPINPROJ_USE_DAGGER
 comm_buffer.wait_sending();

#endif // USE_MPI 

  // final barrier
#pragma omp barrier

}


//====================================================================
template<>
double AFopr_Domainwall_5din<AFIELDspec>::flop_count(std::string mode)
{
  int    Lvol = CommonParameters::Lvol();
  double flop_site_D, flop_site_Ddag;

  double flop_site_prec=0;  // prec is not yet ready

  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  if (m_repr == "Dirac") {
    vout.crucial("flop_count for Dirac rep. is not yet ready\n");
    return 0;
  } else if (m_repr == "Chiral") {
    flop_site_D = static_cast<double>(
                                      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))  // Wilson
                                      + m_Nc * m_Nd * 2 * (2*9+1+2)   // algebra in 5th direction
                                      );
    flop_site_D *= m_Ns;

    flop_site_Ddag = static_cast<double>(
                                         m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))  // Wilson
                                         + m_Nc * m_Nd * 2 * (2*2+1)   // algebra in 5th direction
                                         );
    flop_site_Ddag *= m_Ns;

  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n", class_name.c_str());
    abort();
  }

  double flop=0;
  if ( (mode == "D") || (mode == "H") ){
    flop=flop_site_D;
  } else if ( (mode == "Ddag") || (mode=="Hdag")) {
    flop=flop_site_Ddag;
  } else if ((mode == "DdagD") || (mode == "DDdag")) {
    flop=flop_site_D + flop_site_Ddag;
  } else if (mode == "D_prec") {
    flop=flop_site_D+flop_site_prec;
  } else if (mode == "Ddag_prec") {
    flop=flop_site_D+flop_site_prec;
  } else if (mode == "DdagD_prec") {
    flop=flop_site_D+flop_site_Ddag+flop_site_prec;
  } else if ( (mode == "Prec") || (mode == "Precdag") ) {
    flop=flop_site_prec;
  } else {
    vout.general("warning: %s unknown mode for flop_count(): %s\n", class_name.c_str(), mode.c_str());
  }
  flop*=static_cast<double>(Lvol);
  return flop;
}


#endif  // end of the icnlude gurad

