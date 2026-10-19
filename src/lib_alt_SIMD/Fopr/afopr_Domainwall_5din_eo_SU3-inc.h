/*!
        @file    afopr_Domainwall_5din_eo_SU3-inc.h

        @brief   Domainwall Dirac operator, even-odd SU(3) implemntation (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $

 */
//====================================================================
#ifndef _AFOPR_DOMAINWALL_5DIN_EO_SU3_INC_H_
#define _AFOPR_DOMAINWALL_5DIN_EO_SU3_INC_H_

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
  using vcomplex_t = AFIELDspec::vcomplex_t;
  using vreal_t = AFIELDspec::vreal_t;

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

  //====================================================================
double flop_count(const std::string &repr,
                  const std::string &mode,
                  const Bridge::VerboseLevel &vl,
                  const std::string &class_name)
{
  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]
  // Modified for optimal clover mult implementation [23 Feb 2018 I.Kanamori]

  static const int Nc=3;
  static const int Nd=4;
  if (repr == "Dirac") {
    flop_site = static_cast<double>(
      Nc * Nd * (4 + 6 * (4 * Nc + 2) + 2 * (4 * Nc + 1))
      + 8 * Nc * Nc * Nd * Nd/2  // <- clover term as chiral
      - 4 * Nc * Nd                  // <- diagonal is real
      + 4 * Nc * Nd );               // <- trans. chiral <--> dirac
  } else if (repr == "Chiral") {
    flop_site = static_cast<double>(
      Nc * Nd * (4 + 8 * (4 * Nc + 2))
      + 8 * Nc * Nc * Nd * Nd/ 2   // <- clover term
      - 4 * Nc * Nd);                  // <- diagonal is real
  } else {
    vout.crucial(vl, "%s: input repr is undefined: repr=%s\n",
                 class_name.c_str(), repr.c_str());
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop=0;
}


  void accum_aPm5_Chiral(vcomplex_t &out_sp1, vcomplex_t &out_sp2,
                         vcomplex_t &out_sp3, vcomplex_t &out_sp4,
                         const vcomplex_t &va,
                         const vcomplex_t sp1, const vcomplex_t sp2,
                         const vcomplex_t sp3, const vcomplex_t sp4) {
    vcomplex_t tmp1 = sub(out_sp1, sp1);;
    vcomplex_t tmp2 = sub(out_sp2, sp2);;
    vcomplex_t tmp3 = add(out_sp3, sp3);;
    vcomplex_t tmp4 = add(out_sp4, sp4);;
    axpy(out_sp1, va, tmp1);
    axpy(out_sp2, va, tmp2);
    axpy(out_sp3, va, tmp3);
    axpy(out_sp4, va, tmp4);
  }

  void accum_scal_Pm5_Chiral(vcomplex_t &out_sp1, vcomplex_t &out_sp2,
                             vcomplex_t &out_sp3, vcomplex_t &out_sp4,
                             const vreal_t &va,
                             const vcomplex_t sp1, const vcomplex_t sp2,
                             const vcomplex_t sp3, const vcomplex_t sp4) {
    // out += a (1-gm5)in
    vcomplex_t tmp3 = add(sp3, sp3);
    vcomplex_t tmp4 = add(sp4, sp4);
    axpy(out_sp3, va, tmp3);
    axpy(out_sp4, va, tmp4);
  }

  void accum_scal_Pp5_Chiral(vcomplex_t &out_sp1, vcomplex_t &out_sp2,
                             vcomplex_t &out_sp3, vcomplex_t &out_sp4,
                             const vreal_t &va,
                             const vcomplex_t sp1, const vcomplex_t sp2,
                             const vcomplex_t sp3, const vcomplex_t sp4) {
    vcomplex_t tmp1 = add(sp1, sp1);
    vcomplex_t tmp2 = add(sp2, sp2);
    axpy(out_sp1, va, tmp1);
    axpy(out_sp2, va, tmp2);
  }

} // anonymous namespace


//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_U_inv(AFIELDspec &v,
                                                      const AFIELDspec &w)
{
  using AFIELD=AFIELDspec;

  const vcomplex_t *pw=w.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  constexpr int NS=4;
  constexpr int NCOL=3;
  const int Ns=m_Ns;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, isite, nsite;
  set_threadtask(ith, nth, isite, nsite, Rst2);

  const real_t* dpinv = &m_dpinv_array[ith][0];
  const real_t*    dm = &m_dm_array[ith][0];
  const real_t*     f = &m_f_array[ith][0];

  for(int site = isite; site < nsite; ++site){
    size_t offset=site*Ns*NS*NCOL;
    const vcomplex_t *in0=pw+offset;
    vcomplex_t *out0=pv+offset;

    real_t a=dpinv[Ns-1]; // chiral
    vreal_t va;
    bcast(va, a);
    const vcomplex_t *in = in0 + (Ns-1)*NS*NCOL;
    vcomplex_t *out = out0 + (Ns-1)*NS*NCOL;
    vcomplex_t xc0_sp1 = in[0];
    vcomplex_t xc1_sp1 = in[1];
    vcomplex_t xc2_sp1 = in[2];
    vcomplex_t xc0_sp2 = in[3];
    vcomplex_t xc1_sp2 = in[4];
    vcomplex_t xc2_sp2 = in[5];
    vcomplex_t xc0_sp3 = in[6];
    vcomplex_t xc1_sp3 = in[7];
    vcomplex_t xc2_sp3 = in[8];
    vcomplex_t xc0_sp4 = in[9];
    vcomplex_t xc1_sp4 = in[10];
    vcomplex_t xc2_sp4 = in[11];


    scal(xc0_sp1, va);
    scal(xc1_sp1, va);
    scal(xc2_sp1, va);
    scal(xc0_sp2, va);
    scal(xc1_sp2, va);
    scal(xc2_sp2, va);
    scal(xc0_sp3, va);
    scal(xc1_sp3, va);
    scal(xc2_sp3, va);
    scal(xc0_sp4, va);
    scal(xc1_sp4, va);
    scal(xc2_sp4, va);

    out[0] = xc0_sp1;
    out[1] = xc1_sp1;
    out[2] = xc2_sp1;
    out[3] = xc0_sp2;
    out[4] = xc1_sp2;
    out[5] = xc2_sp2;
    out[6] = xc0_sp3;
    out[7] = xc1_sp3;
    out[8] = xc2_sp3;
    out[9] = xc0_sp4;
    out[10]= xc1_sp4;
    out[11]= xc2_sp4;


    for(int is = Ns-2; is >= 0; --is){
      in -= NS*NCOL;
      out -= NS*NCOL;
      real_t a  = dpinv[is];
      bcast(va, a);
      { // chirality + part
        real_t ap = -f[is];
        vreal_t vap;
        bcast(vap, ap);
        vcomplex_t yc0_sp1=in[0];
        vcomplex_t yc1_sp1=in[1];
        vcomplex_t yc2_sp1=in[2];
        vcomplex_t yc0_sp2=in[3];
        vcomplex_t yc1_sp2=in[4];
        vcomplex_t yc2_sp2=in[5];

        axpy(yc0_sp1, vap, xc0_sp1); // yc +=ap(1+gm5)/2 xn
        axpy(yc1_sp1, vap, xc1_sp1);
        axpy(yc2_sp1, vap, xc2_sp1);
        axpy(yc0_sp2, vap, xc0_sp2);
        axpy(yc1_sp2, vap, xc1_sp2);
        axpy(yc2_sp2, vap, xc2_sp2);
        scal(yc0_sp1, va);
        scal(yc1_sp1, va);
        scal(yc2_sp1, va);
        scal(yc0_sp2, va);
        scal(yc1_sp2, va);
        scal(yc2_sp2, va);

        out[0] = yc0_sp1;
        out[1] = yc1_sp1;
        out[2] = yc2_sp1;
        out[3] = yc0_sp2;
        out[4] = yc1_sp2;
        out[5] = yc2_sp2;
      }

      {  // chirality - part
        real_t am = dm[is];
        vreal_t vam;
        bcast(vam, am);

        vcomplex_t yc0_sp3=in[6];
        vcomplex_t yc1_sp3=in[7];
        vcomplex_t yc2_sp3=in[8];
        vcomplex_t yc0_sp4=in[9];
        vcomplex_t yc1_sp4=in[10];
        vcomplex_t yc2_sp4=in[11];

        aypx(xc0_sp3, vam, yc0_sp3); // yc +=am(1-gm5)/2 vc
        aypx(xc1_sp3, vam, yc1_sp3);
        aypx(xc2_sp3, vam, yc2_sp3);
        aypx(xc0_sp4, vam, yc0_sp4);
        aypx(xc1_sp4, vam, yc1_sp4);
        aypx(xc2_sp4, vam, yc2_sp4);

        scal(xc0_sp3, va);
        scal(xc1_sp3, va);
        scal(xc2_sp3, va);
        scal(xc0_sp4, va);
        scal(xc1_sp4, va);
        scal(xc2_sp4, va);

        out[6]  = xc0_sp3;
        out[7]  = xc1_sp3;
        out[8]  = xc2_sp3;
        out[9]  = xc0_sp4;
        out[10] = xc1_sp4;
        out[11] = xc2_sp4;
      }
    } // is
  } //site

}


//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_Udag_inv(AFIELDspec &v,
                                                      const AFIELDspec &w)
{
  using AFIELD=AFIELDspec;

  const vcomplex_t *pw=w.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  constexpr int NS=4;
  constexpr int NCOL=3;
  const int Ns=m_Ns;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, isite, nsite;
  set_threadtask(ith, nth, isite, nsite, Rst2);

  const real_t* dpinv = &m_dpinv_array[ith][0];
  const real_t*    dm = &m_dm_array[ith][0];
  const real_t*     f = &m_f_array[ith][0];

  for(int site = isite; site < nsite; ++site){
    size_t offset=site*Ns*NS*NCOL;
    const vcomplex_t *in0 = pw+offset;
    vcomplex_t *out0 = pv+offset;

    real_t a=dpinv[0];
    vreal_t va;
    bcast(va, a);
    const vcomplex_t *in = in0;
    vcomplex_t *out = out0;
    vcomplex_t xc0_sp1 = in[0];
    vcomplex_t xc1_sp1 = in[1];
    vcomplex_t xc2_sp1 = in[2];
    vcomplex_t xc0_sp2 = in[3];
    vcomplex_t xc1_sp2 = in[4];
    vcomplex_t xc2_sp2 = in[5];
    vcomplex_t xc0_sp3 = in[6];
    vcomplex_t xc1_sp3 = in[7];
    vcomplex_t xc2_sp3 = in[8];
    vcomplex_t xc0_sp4 = in[9];
    vcomplex_t xc1_sp4 = in[10];
    vcomplex_t xc2_sp4 = in[11];

    scal(xc0_sp1, va);
    scal(xc1_sp1, va);
    scal(xc2_sp1, va);
    scal(xc0_sp2, va);
    scal(xc1_sp2, va);
    scal(xc2_sp2, va);
    scal(xc0_sp3, va);
    scal(xc1_sp3, va);
    scal(xc2_sp3, va);
    scal(xc0_sp4, va);
    scal(xc1_sp4, va);
    scal(xc2_sp4, va);

    out[0]  = xc0_sp1;
    out[1]  = xc1_sp1;
    out[2]  = xc2_sp1;
    out[3]  = xc0_sp2;
    out[4]  = xc1_sp2;
    out[5]  = xc2_sp2;
    out[6]  = xc0_sp3;
    out[7]  = xc1_sp3;
    out[8]  = xc2_sp3;
    out[9]  = xc0_sp4;
    out[10] = xc1_sp4;
    out[11] = xc2_sp4;


    real_t ap = f[0];
    vreal_t vap;
    bcast(vap, ap);
    scal(xc0_sp1, vap);
    scal(xc1_sp1, vap);
    scal(xc2_sp1, vap);
    scal(xc0_sp2, vap);
    scal(xc1_sp2, vap);
    scal(xc2_sp2, vap);

    for(int is = 1; is < Ns-1; ++is){
      in += NS*NCOL;
      out += NS*NCOL;
      real_t a  = dpinv[is];
      bcast(va, a);

      { // chirality + part
        vcomplex_t yc0_sp1 = in[0];
        vcomplex_t yc1_sp1 = in[1];
        vcomplex_t yc2_sp1 = in[2];
        vcomplex_t yc0_sp2 = in[3];
        vcomplex_t yc1_sp2 = in[4];
        vcomplex_t yc2_sp2 = in[5];

        scal(yc0_sp1, va);
        scal(yc1_sp1, va);
        scal(yc2_sp1, va);
        scal(yc0_sp2, va);
        scal(yc1_sp2, va);
        scal(yc2_sp2, va);

        out[0] = yc0_sp1;
        out[1] = yc1_sp1;
        out[2] = yc2_sp1;
        out[3] = yc0_sp2;
        out[4] = yc1_sp2;
        out[5] = yc2_sp2;

        real_t ap = f[is];
        vreal_t vap;
        bcast(vap, ap);
        axpy(xc0_sp1, vap, yc0_sp1); // x +=ap(1+gm5)/2 y
        axpy(xc1_sp1, vap, yc1_sp1);
        axpy(xc2_sp1, vap, yc2_sp1);
        axpy(xc0_sp2, vap, yc0_sp2);
        axpy(xc1_sp2, vap, yc1_sp2);
        axpy(xc2_sp2, vap, yc2_sp2);
      }

      {  // chirality - part
        real_t am = dm[is-1];
        vreal_t vam;
        bcast(vam, am);

        vcomplex_t yc0_sp3 = in[6];
        vcomplex_t yc1_sp3 = in[7];
        vcomplex_t yc2_sp3 = in[8];
        vcomplex_t yc0_sp4 = in[9];
        vcomplex_t yc1_sp4 = in[10];
        vcomplex_t yc2_sp4 = in[11];

        aypx(xc0_sp3, vam, yc0_sp3); // vc = am vc + vy
        aypx(xc1_sp3, vam, yc1_sp3);
        aypx(xc2_sp3, vam, yc2_sp3);
        aypx(xc0_sp4, vam, yc0_sp4);
        aypx(xc1_sp4, vam, yc1_sp4);
        aypx(xc2_sp4, vam, yc2_sp4);
        scal(xc0_sp3, va);
        scal(xc1_sp3, va);
        scal(xc2_sp3, va);
        scal(xc0_sp4, va);
        scal(xc1_sp4, va);
        scal(xc2_sp4, va);

        out[6]  = xc0_sp3;
        out[7]  = xc1_sp3;
        out[8]  = xc2_sp3;
        out[9]  = xc0_sp4;
        out[10] = xc1_sp4;
        out[11] = xc2_sp4;
      }
    } // is


    in += NS*NCOL;
    out += NS*NCOL;
    const int is=Ns-1;
    a = dpinv[is];
    bcast(va, a);
    { // chirality + part
      bcast(vap, ap);
      vcomplex_t yc0_sp1=in[0];
      vcomplex_t yc1_sp1=in[1];
      vcomplex_t yc2_sp1=in[2];
      vcomplex_t yc0_sp2=in[3];
      vcomplex_t yc1_sp2=in[4];
      vcomplex_t yc2_sp2=in[5];

      accum_sub(yc0_sp1, xc0_sp1);
      accum_sub(yc1_sp1, xc1_sp1);
      accum_sub(yc2_sp1, xc2_sp1);
      accum_sub(yc0_sp2, xc0_sp2);
      accum_sub(yc1_sp2, xc1_sp2);
      accum_sub(yc2_sp2, xc2_sp2);

      scal(yc0_sp1, va);
      scal(yc1_sp1, va);
      scal(yc2_sp1, va);
      scal(yc0_sp2, va);
      scal(yc1_sp2, va);
      scal(yc2_sp2, va);

      out[0] = yc0_sp1;
      out[1] = yc1_sp1;
      out[2] = yc2_sp1;
      out[3] = yc0_sp2;
      out[4] = yc1_sp2;
      out[5] = yc2_sp2;
    }

    {  // chirality - part
      real_t am = dm[is-1];
      vreal_t vam;
      bcast(vam, am);

      vcomplex_t yc0_sp3=in[6];
      vcomplex_t yc1_sp3=in[7];
      vcomplex_t yc2_sp3=in[8];
      vcomplex_t yc0_sp4=in[9];
      vcomplex_t yc1_sp4=in[10];
      vcomplex_t yc2_sp4=in[11];

      aypx(xc0_sp3, vam, yc0_sp3); // xc = am xc + yc
      aypx(xc1_sp3, vam, yc1_sp3);
      aypx(xc2_sp3, vam, yc2_sp3);
      aypx(xc0_sp4, vam, yc0_sp4);
      aypx(xc1_sp4, vam, yc1_sp4);
      aypx(xc2_sp4, vam, yc2_sp4);
      scal(xc0_sp3, va);
      scal(xc1_sp3, va);
      scal(xc2_sp3, va);
      scal(xc0_sp4, va);
      scal(xc1_sp4, va);
      scal(xc2_sp4, va);

      out[6]  = xc0_sp3;
      out[7]  = xc1_sp3;
      out[8]  = xc2_sp3;
      out[9]  = xc0_sp4;
      out[10] = xc1_sp4;
      out[11] = xc2_sp4;

    }

  } //site

}


//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_L_inv(AFIELDspec &v,
                                                      const AFIELDspec &w)
{
  using AFIELD=AFIELDspec;

  const vcomplex_t *pw=w.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  constexpr int NS=4;
  constexpr int NCOL=3;
  const int Ns=m_Ns;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, isite, nsite;
  set_threadtask(ith, nth, isite, nsite, Rst2);

  const real_t* dpinv = &m_dpinv_array[ith][0];
  const real_t*    dm = &m_dm_array[ith][0];
  const real_t*     e = &m_e_array[ith][0];

  for(int site = isite; site < nsite; ++site){
    size_t offset=site*Ns*NS*NCOL;
    const vcomplex_t *in0=pw+offset;
    vcomplex_t *out0=pv+offset;

    const vcomplex_t *in = in0;
    vcomplex_t *out = out0;
    vcomplex_t xc0_sp1 = in[0];
    vcomplex_t xc1_sp1 = in[1];
    vcomplex_t xc2_sp1 = in[2];
    vcomplex_t xc0_sp2 = in[3];
    vcomplex_t xc1_sp2 = in[4];
    vcomplex_t xc2_sp2 = in[5];
    vcomplex_t xc0_sp3 = in[6];
    vcomplex_t xc1_sp3 = in[7];
    vcomplex_t xc2_sp3 = in[8];
    vcomplex_t xc0_sp4 = in[9];
    vcomplex_t xc1_sp4 = in[10];
    vcomplex_t xc2_sp4 = in[11];

    out[0] = xc0_sp1;
    out[1] = xc1_sp1;
    out[2] = xc2_sp1;
    out[3] = xc0_sp2;
    out[4] = xc1_sp2;
    out[5] = xc2_sp2;
    out[6] = xc0_sp3;
    out[7] = xc1_sp3;
    out[8] = xc2_sp3;
    out[9] = xc0_sp4;
    out[10]= xc1_sp4;
    out[11]= xc2_sp4;

    // P_- for chiral rep.
    real_t a=e[0];
    vreal_t va;
    bcast(va, a);
    scal(xc0_sp3, va);
    scal(xc1_sp3, va);
    scal(xc2_sp3, va);
    scal(xc0_sp4, va);
    scal(xc1_sp4, va);
    scal(xc2_sp4, va);

    for(int is = 1; is < Ns-1; ++is){
      in += NS*NCOL;
      out += NS*NCOL;
      { // chirality + part
        real_t ap = dm[is]*dpinv[is-1];
        vreal_t vap;
        bcast(vap, ap);
        vcomplex_t yc0_sp1=in[0];
        vcomplex_t yc1_sp1=in[1];
        vcomplex_t yc2_sp1=in[2];
        vcomplex_t yc0_sp2=in[3];
        vcomplex_t yc1_sp2=in[4];
        vcomplex_t yc2_sp2=in[5];

        aypx(xc0_sp1, vap, yc0_sp1); // yc +=ap(1+gm5)/2 xc
        aypx(xc1_sp1, vap, yc1_sp1);
        aypx(xc2_sp1, vap, yc2_sp1);
        aypx(xc0_sp2, vap, yc0_sp2);
        aypx(xc1_sp2, vap, yc1_sp2);
        aypx(xc2_sp2, vap, yc2_sp2);

        out[0] = xc0_sp1;
        out[1] = xc1_sp1;
        out[2] = xc2_sp1;
        out[3] = xc0_sp2;
        out[4] = xc1_sp2;
        out[5] = xc2_sp2;
      }
      {  // chirality - part
        real_t am = e[is];
        vreal_t vam;
        bcast(vam, am);

        vcomplex_t yc0_sp3 = in[6];
        vcomplex_t yc1_sp3 = in[7];
        vcomplex_t yc2_sp3 = in[8];
        vcomplex_t yc0_sp4 = in[9];
        vcomplex_t yc1_sp4 = in[10];
        vcomplex_t yc2_sp4 = in[11];

        out[6]  = yc0_sp3;
        out[7]  = yc1_sp3;
        out[8]  = yc2_sp3;
        out[9]  = yc0_sp4;
        out[10] = yc1_sp4;
        out[11] = yc2_sp4;

        axpy(xc0_sp3, vam, yc0_sp3);
        axpy(xc1_sp3, vam, yc1_sp3);
        axpy(xc2_sp3, vam, yc2_sp3);
        axpy(xc0_sp4, vam, yc0_sp4);
        axpy(xc1_sp4, vam, yc1_sp4);
        axpy(xc2_sp4, vam, yc2_sp4);
      }
    } // is

    in += NS*NCOL;
    out += NS*NCOL;
    const int is=Ns-1;
    { // chirality + part
      real_t ap = dm[is]*dpinv[is-1];
      vreal_t vap;
      bcast(vap, ap);
      vcomplex_t yc0_sp1 = in[0];
      vcomplex_t yc1_sp1 = in[1];
      vcomplex_t yc2_sp1 = in[2];
      vcomplex_t yc0_sp2 = in[3];
      vcomplex_t yc1_sp2 = in[4];
      vcomplex_t yc2_sp2 = in[5];

      aypx(xc0_sp1, vap, yc0_sp1);
      aypx(xc1_sp1, vap, yc1_sp1);
      aypx(xc2_sp1, vap, yc2_sp1);
      aypx(xc0_sp2, vap, yc0_sp2);
      aypx(xc1_sp2, vap, yc1_sp2);
      aypx(xc2_sp2, vap, yc2_sp2);

      out[0] = xc0_sp1;
      out[1] = xc1_sp1;
      out[2] = xc2_sp1;
      out[3] = xc0_sp2;
      out[4] = xc1_sp2;
      out[5] = xc2_sp2;
    }
    {  // chirality - part
        vcomplex_t yc0_sp3 = in[6];
        vcomplex_t yc1_sp3 = in[7];
        vcomplex_t yc2_sp3 = in[8];
        vcomplex_t yc0_sp4 = in[9];
        vcomplex_t yc1_sp4 = in[10];
        vcomplex_t yc2_sp4 = in[11];

        accum_sub(yc0_sp3, xc0_sp3);
        accum_sub(yc1_sp3, xc1_sp3);
        accum_sub(yc2_sp3, xc2_sp3);
        accum_sub(yc0_sp4, xc0_sp4);
        accum_sub(yc1_sp4, xc1_sp4);
        accum_sub(yc2_sp4, xc2_sp4);

        out[6]  = yc0_sp3;
        out[7]  = yc1_sp3;
        out[8]  = yc2_sp3;
        out[9]  = yc0_sp4;
        out[10] = yc1_sp4;
        out[11] = yc2_sp4;
    }

  } //site

}


//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_Ldag_inv(AFIELDspec &v,
                                                      const AFIELDspec &w)
{
  using AFIELD=AFIELDspec;

  const vcomplex_t *pw=w.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  constexpr int NS=4;
  constexpr int NCOL=3;
  const int Ns=m_Ns;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, isite, nsite;
  set_threadtask(ith, nth, isite, nsite, Rst2);

  const real_t* dpinv = &m_dpinv_array[ith][0];
  const real_t*    dm = &m_dm_array[ith][0];
  const real_t*     e = &m_e_array[ith][0];

  for(int site = isite; site < nsite; ++site){
    size_t offset=site*Ns*NS*NCOL;
    const vcomplex_t *in0=pw+offset;
    vcomplex_t *out0=pv+offset;

    const vcomplex_t *in = in0 + (Ns-1)*NS*NCOL;
    vcomplex_t *out = out0 + (Ns-1)*NS*NCOL;
    vcomplex_t xc0_sp1 = in[0];
    vcomplex_t xc1_sp1 = in[1];
    vcomplex_t xc2_sp1 = in[2];
    vcomplex_t xc0_sp2 = in[3];
    vcomplex_t xc1_sp2 = in[4];
    vcomplex_t xc2_sp2 = in[5];
    vcomplex_t xc0_sp3 = in[6];
    vcomplex_t xc1_sp3 = in[7];
    vcomplex_t xc2_sp3 = in[8];
    vcomplex_t xc0_sp4 = in[9];
    vcomplex_t xc1_sp4 = in[10];
    vcomplex_t xc2_sp4 = in[11];

    out[0]  = xc0_sp1;
    out[1]  = xc1_sp1;
    out[2]  = xc2_sp1;
    out[3]  = xc0_sp2;
    out[4]  = xc1_sp2;
    out[5]  = xc2_sp2;
    out[6]  = xc0_sp3;
    out[7]  = xc1_sp3;
    out[8]  = xc2_sp3;
    out[9]  = xc0_sp4;
    out[10] = xc1_sp4;
    out[11] = xc2_sp4;


    for(int is = Ns-2; is > -1 ; --is){
      in -= NS*NCOL;
      out -= NS*NCOL;
      { // chirality + part
        real_t ap = dm[is+1]*dpinv[is];
        vreal_t vap;
        bcast(vap, ap);
        vcomplex_t yc0_sp1=in[0];
        vcomplex_t yc1_sp1=in[1];
        vcomplex_t yc2_sp1=in[2];
        vcomplex_t yc0_sp2=in[3];
        vcomplex_t yc1_sp2=in[4];
        vcomplex_t yc2_sp2=in[5];

        aypx(xc0_sp1, vap, yc0_sp1); // x = ap *x  + y
        aypx(xc1_sp1, vap, yc1_sp1);
        aypx(xc2_sp1, vap, yc2_sp1);
        aypx(xc0_sp2, vap, yc0_sp2);
        aypx(xc1_sp2, vap, yc1_sp2);
        aypx(xc2_sp2, vap, yc2_sp2);

        out[0] = xc0_sp1;
        out[1] = xc1_sp1;
        out[2] = xc2_sp1;
        out[3] = xc0_sp2;
        out[4] = xc1_sp2;
        out[5] = xc2_sp2;

      }

      {  // chirality - part
        real_t am = -e[is];
        vreal_t vam;
        bcast(vam, am);

        vcomplex_t yc0_sp3=in[6];
        vcomplex_t yc1_sp3=in[7];
        vcomplex_t yc2_sp3=in[8];
        vcomplex_t yc0_sp4=in[9];
        vcomplex_t yc1_sp4=in[10];
        vcomplex_t yc2_sp4=in[11];

        axpy(yc0_sp3, vam, xc0_sp3);
        axpy(yc1_sp3, vam, xc1_sp3);
        axpy(yc2_sp3, vam, xc2_sp3);
        axpy(yc0_sp4, vam, xc0_sp4);
        axpy(yc1_sp4, vam, xc1_sp4);
        axpy(yc2_sp4, vam, xc2_sp4);

        out[6]  = yc0_sp3;
        out[7]  = yc1_sp3;
        out[8]  = yc2_sp3;
        out[9]  = yc0_sp4;
        out[10] = yc1_sp4;
        out[11] = yc2_sp4;
      }
    } // is

  } //site

}

//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_D_eo(AFIELDspec &v, const AFIELDspec &w,
                                                     const int ieo)
{
  typedef AFIELDspec AFIELD;

#ifndef AFOPR_DOMAINWALL_5DIN_EO_FORCE_COMMUNICATE
  assert( !m_bools.has_mpi_x);
#endif

  // reduced, after divided by the simd layout
  const int Rst=m_Rst;
  const int Rx=m_Rx;
  const int Ry=m_Ry;
  const int Rz=m_Rz;
  const int Rt=m_Rt;

  const int bulk_size=m_Rst2;
  const int boundary_size=m_Rst_boundary2;
  assert(CommonParameters::Nc() == 3);
  constexpr int Ndf_ferm=12;
  const int Rx2=m_Rx/2;

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

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

  // input/output fermion field
  const vcomplex_t *pw_full=w.ptr_vcomplex();
  vcomplex_t *pv_full=v.ptr_vcomplex();

  // gauge field
  vcomplex_t *Uptr_p, *Uptr_m;
  if(ieo==0){
    Uptr_p=m_Ue.ptr_vcomplex();
    Uptr_m=m_Uo.ptr_vcomplex();
  } else {
    Uptr_p=m_Uo.ptr_vcomplex();
    Uptr_m=m_Ue.ptr_vcomplex();
  }

  // list vectors
  const site_list_t *list_bulk=m_site_list_bulk_ptr[ieo];

#ifdef USE_MPI
  const site_list_t *list_boundary=m_site_list_boundary_ptr[ieo];
  const site_list_pack_t *list_pack_x=m_site_list_pack_x_ptr[ieo];
  const site_list_pack_t *list_pack_y=m_site_list_pack_y_ptr[ieo];
  const site_list_pack_t *list_pack_z=m_site_list_pack_z_ptr[ieo];
  const site_list_pack_t *list_pack_t=m_site_list_pack_t_ptr[ieo];

  // initiate recieving
  comm_buffer.start_receiving();
#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
  {
    timer_calc0->start();
  } // master
#endif // Timer



#define PACK_TO_BUFFER

  if(bools.has_mpi_y){
#ifdef RUN_HOP_YP
    for(int i=isy; i<nsy; ++i){
      const int ss = list_pack_y[i].here;
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
      const int ss = list_pack_y[i].here + Rx2*(Ry-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
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
      const int ss = list_pack_z[i].here;
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
      const int ss = list_pack_z[i].here+Rx2*Ry*(Rz-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
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
      const int ss = list_pack_t[i].here;
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
      const int ss = list_pack_t[i].here + Rx2*Ry*Rz*(Rt-1);
      const vcomplex_t *in = pw_full + ss*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
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
#pragma omp master
  {
#ifdef TimerAfoprDomainwall_5din_eo
    timer_calc0->stop();
#endif // Timer
#ifdef TimerAfoprDomainwall_5din_eo_Comm
    timer_comm->start();
#endif
  }
  comm_buffer.start_sending();
  comm_buffer.test_all();
#endif  // USE_MPI

#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
   timer_calc1->start();
 }
#endif

 for(int s0=is_bulk; s0<ns_bulk; ++s0){

    // get coordinates
    const int s = list_bulk[s0].here;
    const int ieo_work=list_bulk[s0].ieo_work;
    const int sn_xp = list_bulk[s0].neighbor[DirXp];
    const int sn_xm = list_bulk[s0].neighbor[DirXm];
    const int sn_yp = list_bulk[s0].neighbor[DirYp];
    const int sn_ym = list_bulk[s0].neighbor[DirYm];
    const int sn_zp = list_bulk[s0].neighbor[DirZp];
    const int sn_zm = list_bulk[s0].neighbor[DirZm];
    const int sn_tp = list_bulk[s0].neighbor[DirTp];
    const int sn_tm = list_bulk[s0].neighbor[DirTm];

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

    { // source at s
      const vcomplex_t *in=in_base+s*Ndf_ferm5;
      prewilson_at_site(pw0_site, in, Ns, c_array, mq);
    }


    /********************************************************************/
    // + directions
    //
    const vcomplex_t *Uptr=Uptr_p;
    // x direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_XP
    if(sn_xp>-1){
      const vcomplex_t *in=in_base+sn_xp*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_xp
      const vcomplex_t *pw0 = pw0_site;  // at s
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u, Uptr, s, DirX );
      for(int s5=0; s5<Ns; s5++){
#define HOP_DOES_NOT_ACCUM
      if(ieo_work==0){
#include "Wilson_inc/multi_hop_xp_eo_impl1-inc.h"
      }else {
#include "Wilson_inc/multi_hop_xp_impl1-inc.h"
      }
#undef HOP_DOES_NOT_ACCUM
        pw0 += Ndf_ferm;
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif // RUN_HOP_XP
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YP
    if(sn_yp>-1){
      const vcomplex_t *in=in_base+sn_yp*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_yp
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u, Uptr, s, DirY );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_yp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZP
    if(sn_zp>-1){
      const vcomplex_t *in=in_base+sn_zp*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_zp
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u, Uptr, s, DirZ );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_zp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TP
    if(sn_tp>-1){
      const vcomplex_t *in=in_base+sn_tp*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_tp
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u, Uptr, s, DirT );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_tp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif

    /********************************************************************/
    // - directions
    //
    Uptr=Uptr_m;
    // x direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_XM
    if(sn_xm>-1){
      const vcomplex_t *in=in_base+sn_xm*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_xm
      const vcomplex_t *pw0 = pw0_site;  // at s
      vcomplex_t *pv = pv_site;
      if(ieo_work==0){
        GET_GAUGE(u,  Uptr, s, DirX );
        GET_GAUGE(un, Uptr, sn_xm, DirX );
        SHIFT_AND_MERGE_GAUGE(xm, u, un); // (u,un) --> u
        for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_xm_impl1-inc.h"
          pw0 += Ndf_ferm;
          pw  += Ndf_ferm;
          pv  += Ndf_ferm;
        }
      } else {
        GET_GAUGE(u,  Uptr, s, DirX );
        for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_xm_eo_impl1-inc.h"
          pw0 += Ndf_ferm;
          pw  += Ndf_ferm;
          pv  += Ndf_ferm;
        }
      }
    }
#endif
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YM
    if(sn_ym>-1){
      const vcomplex_t *in=in_base+sn_ym*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_ym
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u,  Uptr, sn_ym, DirY );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_ym-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZM
    if(sn_zm>-1){
      const vcomplex_t *in=in_base+sn_zm*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_zm
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u,  Uptr, sn_zm, DirZ );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_zm-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TM
    if(sn_tm>-1){
      const vcomplex_t *in=in_base+sn_tm*Ndf_ferm5;
      prewilson_at_site(pw_site, in, Ns, c_array, mq);
      const vcomplex_t *pw = pw_site;    // at sn_tm
      vcomplex_t *pv = pv_site;
      GET_GAUGE(u,  Uptr, sn_tm, DirT );
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_tm-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif

    /********************************************************************/
    // accumulate
    {
      const vcomplex_t *pv = pv_site;   // hopping
      for(int s5=0; s5<Ns; s5++){
        real_t bb = -0.5*b_array[s5];
        vreal_t vb;
        bcast(vb, bb);
        GET_SPINOR( pv_site, s5 );
        scal(wc0_sp1, vb);
        scal(wc1_sp1, vb);
        scal(wc2_sp1, vb);
        scal(wc0_sp2, vb);
        scal(wc1_sp2, vb);
        scal(wc2_sp2, vb);
        scal(wc0_sp3, vb);
        scal(wc1_sp3, vb);
        scal(wc2_sp3, vb);
        scal(wc0_sp4, vb);
        scal(wc1_sp4, vb);
        scal(wc2_sp4, vb);
        SET_SPINOR( out, s5);
      } //s5
    }
 } // s0 (loop for bulk site)
#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
   timer_calc1->stop();
 }
#endif

 
 #ifdef USE_MPI
 // finalize the communication, all direction
 comm_buffer.wait_receiving();

#ifdef TimerAfoprDomainwall_5din_eo_Comm
#pragma omp master
 {
   //    Communicator::sync();
    timer_comm->stop();
 }  // end of omp master
#endif // TimerAfoprDomaiwall_5din_eo_Comm
#pragma omp barrier

#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
   timer_calc2->start();
 }
#endif
 for(int s0=is_boundary; s0<ns_boundary; ++s0){

   // get coordinates
   const int s = list_boundary[s0].here;
   const int s_packed_y = list_boundary[s0].packed[DirY];
   const int s_packed_z = list_boundary[s0].packed[DirZ];
   const int s_packed_t = list_boundary[s0].packed[DirT];
   const int sn_xp = list_boundary[s0].neighbor[DirXp];
   const int sn_xm = list_boundary[s0].neighbor[DirXm];
   const int sn_yp = list_boundary[s0].neighbor[DirYp];
   const int sn_ym = list_boundary[s0].neighbor[DirYm];
   const int sn_zp = list_boundary[s0].neighbor[DirZp];
   const int sn_zm = list_boundary[s0].neighbor[DirZm];
   const int sn_tp = list_boundary[s0].neighbor[DirTp];
   const int sn_tm = list_boundary[s0].neighbor[DirTm];

#define NEED_TO_READ_BUFFER
    /********************************************************************/
    // + directions
    //
    const vcomplex_t *Uptr=Uptr_p;
    // y direction //////////////////////////////////////////////////////
  #ifdef RUN_HOP_YP
  if(sn_yp>-1){
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
  }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZP
  if(sn_zp>-1){
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
  }
#endif
  // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TP
  if(sn_tp>-1){
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
  }
#endif

    /********************************************************************/
    // - directions
    //
    Uptr=Uptr_m;
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YM
  if(sn_ym>-1){
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
  }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZM
  if(sn_zm>-1){
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
  }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TM
  if(sn_tm>-1){
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
  }
#endif

 } // s0

#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
    timer_calc2->stop();
 }
#endif

 comm_buffer.wait_sending();

#undef NEED_TO_READ_BUFFER

#endif // USE_MPI 
  // final barrier
#pragma omp barrier

}

//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AFIELDspec>::mult_Ddag_eo(AFIELDspec &v,
                                                        const AFIELDspec &w,
                                                        const int ieo)
{

#define SPINPROJ_USE_DAGGER
  typedef AFIELDspec AFIELD;

#ifndef AFOPR_DOMAINWALL_5DIN_EO_FORCE_COMMUNICATE
  assert( !m_bools.has_mpi_x);
#endif

  // reduced, after divided by the simd layout
  const int Rst=m_Rst;
  const int Rx=m_Rx;
  const int Ry=m_Ry;
  const int Rz=m_Rz;
  const int Rt=m_Rt;

  const int bulk_size=m_Rst2;
  const int boundary_size=m_Rst_boundary2;
  assert(CommonParameters::Nc() == 3);
  constexpr int Ndf_ferm=12;
  const int Rx2=m_Rx/2;

  const int Ns=m_Ns;
  const int Ndf_ferm5 = Ndf_ferm*Ns;

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
  const real_t *cm_array = &m_cm_over_b_array[ith][0];
  const real_t *cp_array = &m_cp_over_b_array[ith][0];
  //const real_t *cm_diag_array = &m_cm_over_b_diag_array[ith][0];
  //  const real_t *cp_diag_array = &m_cp_over_b_diag_array[ith][0];
  const real_t *b_array = &m_b_array[ith][0];
  const real_t mq = m_mq;
  const real_t M0 = m_M0;

  // input/output fermion field
  const vcomplex_t *pw_full=w.ptr_vcomplex();
  vcomplex_t *pv_full=v.ptr_vcomplex();

  // gauge field
  vcomplex_t *Uptr_p, *Uptr_m;
  if(ieo==0){
    Uptr_p=m_Ue.ptr_vcomplex();
    Uptr_m=m_Uo.ptr_vcomplex();
  } else {
    Uptr_p=m_Uo.ptr_vcomplex();
    Uptr_m=m_Ue.ptr_vcomplex();
  }

  // list vectors
  const site_list_t *list_bulk=m_site_list_bulk_ptr[ieo];

#ifdef USE_MPI
  const site_list_t *list_boundary=m_site_list_boundary_ptr[ieo];
  const site_list_pack_t *list_pack_x=m_site_list_pack_x_ptr[ieo];
  const site_list_pack_t *list_pack_y=m_site_list_pack_y_ptr[ieo];
  const site_list_pack_t *list_pack_z=m_site_list_pack_z_ptr[ieo];
  const site_list_pack_t *list_pack_t=m_site_list_pack_t_ptr[ieo];

  // initiate recieving
  comm_buffer.start_receiving();
#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
  {
    timer_calc0->start();
  } // master
#endif // Timer



#define PACK_TO_BUFFER

  if(bools.has_mpi_y){
#ifdef RUN_HOP_YP
    for(int i=isy; i<nsy; ++i){
      const int ss = list_pack_y[i].here;
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
      const int ss = list_pack_y[i].here + Rx2*(Ry-1);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
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
      const int ss = list_pack_z[i].here;
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
      const int ss = list_pack_z[i].here+Rx2*Ry*(Rz-1);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
      GET_GAUGE(u,  Uptr, ss, DirZ );
      const vcomplex_t *pw = pw_full + ss*Ndf_ferm5;
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
      const int ss = list_pack_t[i].here;
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
      const int ss = list_pack_t[i].here + Rx2*Ry*Rz*(Rt-1);
      const vcomplex_t *Uptr=Uptr_m;  // only negative hops use the gauge field
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

#pragma omp barrier
#pragma omp master
  {
#ifdef TimerAfoprDomainwall_5din_eo
    timer_calc0->stop();
#endif
#ifdef TimerAfoprDomainwall_5din_eo_Comm
    timer_comm->start();
#endif
  } // master
  comm_buffer.start_sending();
  comm_buffer.test_all();
#endif  // USE_MPI
#ifdef TimerAfoprDomainwall_5din_eo

#pragma omp master
 {
   timer_calc1->start();
 }
#endif
 for(int s0=is_bulk; s0<ns_bulk; ++s0){

    // get coordinates
    const int s = list_bulk[s0].here;
    const int ieo_work=list_bulk[s0].ieo_work;
    const int sn_xp = list_bulk[s0].neighbor[DirXp];
    const int sn_xm = list_bulk[s0].neighbor[DirXm];
    const int sn_yp = list_bulk[s0].neighbor[DirYp];
    const int sn_ym = list_bulk[s0].neighbor[DirYm];
    const int sn_zp = list_bulk[s0].neighbor[DirZp];
    const int sn_zm = list_bulk[s0].neighbor[DirZm];
    const int sn_tp = list_bulk[s0].neighbor[DirTp];
    const int sn_tm = list_bulk[s0].neighbor[DirTm];

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


    /********************************************************************/
    // + directions
    //
    const vcomplex_t *Uptr=Uptr_p;
    // x direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_XP
    if(sn_xp>-1){
      const vcomplex_t *pw  = in_base+sn_xp*Ndf_ferm5;  // at sn_xp
      const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;     // at s
      GET_GAUGE(u, Uptr, s, DirX );
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#define HOP_DOES_NOT_ACCUM
      if(ieo_work==0){
#include "Wilson_inc/multi_hop_xp_eo_impl1-inc.h"
      }else {
#include "Wilson_inc/multi_hop_xp_impl1-inc.h"
      }
#undef HOP_DOES_NOT_ACCUM
        pw0 += Ndf_ferm;
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif // RUN_HOP_XP
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YP
    if(sn_yp>-1){
      GET_GAUGE(u, Uptr, s, DirY );
      const vcomplex_t *pw  = in_base+sn_yp*Ndf_ferm5;  // at sn_yp
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_yp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZP
    if(sn_zp>-1){
      GET_GAUGE(u, Uptr, s, DirZ );
      const vcomplex_t *pw  = in_base+sn_zp*Ndf_ferm5;  // at sn_zp
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_zp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TP
    if(sn_tp>-1){
      GET_GAUGE(u, Uptr, s, DirT );
      const vcomplex_t *pw  = in_base+sn_tp*Ndf_ferm5;  // at sn_tp
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_tp-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif

    /********************************************************************/
    // - directions
    //
    Uptr=Uptr_m;
    // x direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_XM
    if(sn_xm>-1){
      const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;      // at s
      const vcomplex_t *pw  = in_base+sn_xm*Ndf_ferm5;  // at sn_xm
      vcomplex_t *pv = pv_site;
      if(ieo_work==0){
        GET_GAUGE(u,  Uptr, s, DirX );
        GET_GAUGE(un, Uptr, sn_xm, DirX );
        SHIFT_AND_MERGE_GAUGE(xm, u, un); // (u,un) --> u
        for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_xm_impl1-inc.h"
          pw0 += Ndf_ferm;
          pw  += Ndf_ferm;
          pv  += Ndf_ferm;
        }
      } else {
        GET_GAUGE(u,  Uptr, s, DirX );
        for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_xm_eo_impl1-inc.h"
          pw0 += Ndf_ferm;
          pw  += Ndf_ferm;
          pv  += Ndf_ferm;
        }
      }
    }
#endif
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YM
    if(sn_ym>-1){
      GET_GAUGE(u,  Uptr, sn_ym, DirY );
      const vcomplex_t *pw  = in_base+sn_ym*Ndf_ferm5;  // at sn_ym
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_ym-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZM
    if(sn_zm>-1){
      GET_GAUGE(u,  Uptr, sn_zm, DirZ );
      const vcomplex_t *pw  = in_base+sn_zm*Ndf_ferm5;  // at sn_zm
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_zm-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TM
    if(sn_tm>-1){
      GET_GAUGE(u,  Uptr, sn_tm, DirT );
      const vcomplex_t *pw  = in_base+sn_tm*Ndf_ferm5;  // at sn_tm
      vcomplex_t *pv = pv_site;
      for(int s5=0; s5<Ns; s5++){
#include "Wilson_inc/multi_hop_tm-inc.h"
        pw  += Ndf_ferm;
        pv  += Ndf_ferm;
      }
    }
#endif

    /********************************************************************/
    // accumulate
    {
      prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
      const vcomplex_t *pv = pw0_site;
      for(int s5=0; s5<Ns; s5++){
        real_t bb = -0.5*b_array[s5];
        vreal_t vb;
        bcast(vb, bb);
        GET_SPINOR( pv, s5 );
        scal(wc0_sp1, vb);
        scal(wc1_sp1, vb);
        scal(wc2_sp1, vb);
        scal(wc0_sp2, vb);
        scal(wc1_sp2, vb);
        scal(wc2_sp2, vb);
        scal(wc0_sp3, vb);
        scal(wc1_sp3, vb);
        scal(wc2_sp3, vb);
        scal(wc0_sp4, vb);
        scal(wc1_sp4, vb);
        scal(wc2_sp4, vb);
        SET_SPINOR( out, s5);
        //        pv+=Ndf_ferm;
      } //s5
    }
 } // s0 (loop for bulk site)

 #ifdef USE_MPI

#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
   timer_calc1->stop();
 }
#endif

 // finalize the communication, all direction
 comm_buffer.wait_receiving();
#ifdef TimerAfoprDomainwall_5din_eo_Comm
#pragma omp master
  {
    //    Communicator::sync();
    timer_comm->stop();
  }
#endif

#pragma omp barrier

#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
  {
    timer_calc2->start();
  }
#endif

 for(int s0=is_boundary; s0<ns_boundary; ++s0){

   // get coordinates
   const int s = list_boundary[s0].here;
   const int s_packed_y = list_boundary[s0].packed[DirY];
   const int s_packed_z = list_boundary[s0].packed[DirZ];
   const int s_packed_t = list_boundary[s0].packed[DirT];
   const int sn_xp = list_boundary[s0].neighbor[DirXp];
   const int sn_xm = list_boundary[s0].neighbor[DirXm];
   const int sn_yp = list_boundary[s0].neighbor[DirYp];
   const int sn_ym = list_boundary[s0].neighbor[DirYm];
   const int sn_zp = list_boundary[s0].neighbor[DirZp];
   const int sn_zm = list_boundary[s0].neighbor[DirZm];
   const int sn_tp = list_boundary[s0].neighbor[DirTp];
   const int sn_tm = list_boundary[s0].neighbor[DirTm];

#define NEED_TO_READ_BUFFER
#define HOP_DOES_NOT_ACCUM
   /********************************************************************/
    // + directions
    //
    const vcomplex_t *Uptr=Uptr_p;
    // y direction //////////////////////////////////////////////////////
  #ifdef RUN_HOP_YP
  if(sn_yp>-1){
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

  }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZP
  if(sn_zp>-1){
    vcomplex_t *rbuf=rbuf_zp + s_packed_z * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirZ );
    vcomplex_t *pv = pv_site;
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
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
  }
#endif
  // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TP
  if(sn_tp>-1){
    vcomplex_t *rbuf=rbuf_tp + s_packed_t * Ndf_ferm5/2;
    GET_GAUGE(u, Uptr, s, DirT );
    vcomplex_t *pv = pv_site;
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
    prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
    vcomplex_t *out = pv_full + s*Ndf_ferm5;
    const vcomplex_t *in  = pw0_site;
    for(int i=0; i<Ndf_ferm5; i+=4){
      accum_add(out[i],in[i]);
      accum_add(out[i+1],in[i+1]);
      accum_add(out[i+2],in[i+2]);
      accum_add(out[i+3],in[i+3]);
    }
  }
#endif

    /********************************************************************/
    // - directions
    //
    Uptr=Uptr_m;
    // y direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_YM
  if(sn_ym>-1){
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
  }
#endif
    // z direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_ZM
  if(sn_zm>-1){
    vcomplex_t *rbuf=rbuf_zm + s_packed_z * Ndf_ferm5/2;
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
  }
#endif
    // t direction //////////////////////////////////////////////////////
#ifdef RUN_HOP_TM
  if(sn_tm>-1){
    vcomplex_t *rbuf=rbuf_tm + s_packed_t * Ndf_ferm5/2;
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
  }
#endif

 } // s0


#ifdef TimerAfoprDomainwall_5din_eo
#pragma omp master
 {
    timer_calc2->stop();
 }
#endif

 comm_buffer.wait_sending();

#undef NEED_TO_READ_BUFFER
#undef HOP_DOES_NOT_ACCUM

#endif // USE_MPI 
  // final barrier
#pragma omp barrier
#undef SPINPROJ_USE_DAGGER
}


/*
//====================================================================
template<>
void AFopr_Clover_eo<AFIELDspec>::mult_cswinv(AFIELDspec& v,
                                      const AFIELDspec& w, int ieo)
{

  typedef AFIELDspec AFIELD;
  const vcomplex_t *clov_eo[2];
  clov_eo[0]=m_fee_inv.ptr_vcomplex();
  clov_eo[1]=m_foo_inv.ptr_vcomplex();
  const vcomplex_t *clov_base=clov_eo[ieo];
  const vcomplex_t *pw=w.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  static const int NS=4;
  static const int NCOL=3;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rst2);

  for(int s = is; s < ns; ++s){

    size_t offset=s*NS*NCOL;
    const vcomplex_t *src=pw+offset;
    vcomplex_t *out=pv+offset;

    size_t offset_clov=s*36;

    const vcomplex_t *clov=clov_base+offset_clov;

    vcomplex_t wc0_sp1, wc1_sp1, wc2_sp1;
    vcomplex_t wc0_sp2, wc1_sp2, wc2_sp2;
    vcomplex_t wc0_sp3, wc1_sp3, wc2_sp3;
    vcomplex_t wc0_sp4, wc1_sp4, wc2_sp4;
    #include "Wilson_inc/mult_clov-inc.h"

    sstore(out,   wc0_sp1);
    sstore(out+1, wc1_sp1);
    sstore(out+2, wc2_sp1);
    sstore(out+3, wc0_sp2);
    sstore(out+4, wc1_sp2);
    sstore(out+5, wc2_sp2);
    sstore(out+6, wc0_sp3);
    sstore(out+7, wc1_sp3);
    sstore(out+8, wc2_sp3);
    sstore(out+9, wc0_sp4);
    sstore(out+10, wc1_sp4);
    sstore(out+11, wc2_sp4);
  }

}


//====================================================================
template<>
void AFopr_Clover_eo<AFIELDspec>::sub_mult_cswinv(AFIELDspec& v,
                                                  const AFIELDspec& w1, const AFIELDspec& w2,
                                                  int ieo)
{
  // v = w1 - Tinv w2
  typedef AFIELDspec AFIELD;
  const vcomplex_t *clov_eo[2];
  clov_eo[0]=m_fee_inv.ptr_vcomplex();
  clov_eo[1]=m_foo_inv.ptr_vcomplex();
  const vcomplex_t *clov_base=clov_eo[ieo];
  const vcomplex_t *pw1=w1.ptr_vcomplex();
  const vcomplex_t *pw2=w2.ptr_vcomplex();
  vcomplex_t *pv=v.ptr_vcomplex();

  // reduced, after divided by the simd layout
  const int Rst2=m_Rst2;
  static const int NS=4;
  static const int NCOL=3;
  assert(CommonParameters::Nc() == NCOL);
  assert(CommonParameters::Nd() == NS);
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rst2);

  for(int s = is; s < ns; ++s){

    size_t offset=s*NS*NCOL;
    const vcomplex_t *src=pw2+offset;
    const vcomplex_t *tmp=pw1+offset;
    vcomplex_t *out=pv+offset;

    size_t offset_clov=s*36;

    const vcomplex_t *clov=clov_base+offset_clov;

    // set values from w1
    vcomplex_t wc0_sp1=tmp[0];
    vcomplex_t wc1_sp1=tmp[1];
    vcomplex_t wc2_sp1=tmp[2];
    vcomplex_t wc0_sp2=tmp[3];
    vcomplex_t wc1_sp2=tmp[4];
    vcomplex_t wc2_sp2=tmp[5];
    vcomplex_t wc0_sp3=tmp[6];
    vcomplex_t wc1_sp3=tmp[7];
    vcomplex_t wc2_sp3=tmp[8];
    vcomplex_t wc0_sp4=tmp[9];
    vcomplex_t wc1_sp4=tmp[10];
    vcomplex_t wc2_sp4=tmp[11];
#define MULT_CLOV_SUB  // Tinv w1 is subtracted
#include "Wilson_inc/mult_clov-inc.h"
#undef MULT_CLOV_SUB
    sstore(out,   wc0_sp1);
    sstore(out+1, wc1_sp1);
    sstore(out+2, wc2_sp1);
    sstore(out+3, wc0_sp2);
    sstore(out+4, wc1_sp2);
    sstore(out+5, wc2_sp2);
    sstore(out+6, wc0_sp3);
    sstore(out+7, wc1_sp3);
    sstore(out+8, wc2_sp3);
    sstore(out+9, wc0_sp4);
    sstore(out+10, wc1_sp4);
    sstore(out+11, wc2_sp4);
  }
}

 
//====================================================================
template<>
double AFopr_Clover_eo<AFIELDspec>::flop_count()
{
  return ::flop_count(m_repr, m_mode, m_vl, class_name);
}

*/
#endif
