/*
        @file    $Id: matvec_su3_macro-inc.h #$

        @brief    Macros for SU(3) matrix and vector multplication

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

// v_i = conj(u_{ji})*a_j,
// w_i = conj(u_{ji})*b_j,
#define MULT_DAG_2SPIN(v,w,u,a,b) {             \
    vcomplex_t u11=u[0];                        \
    vcomplex_t u12=u[1];                        \
    vcomplex_t u13=u[2];                        \
    set_cmult_conj(v ## 1, u11, a ## 1);        \
    set_cmult_conj(w ## 1, u11, b ## 1);        \
    vcomplex_t u21=u[3];                        \
    set_cmult_conj(v ## 2, u12, a ## 1);        \
    set_cmult_conj(w ## 2, u12, b ## 1);        \
    vcomplex_t u22=u[4];                        \
    set_cmult_conj(v ## 3, u13, a ## 1);        \
    set_cmult_conj(w ## 3, u13, b ## 1);        \
    vcomplex_t u23=u[5];                        \
    accum_cmult_conj(v ## 1, u21, a ## 2);      \
    accum_cmult_conj(w ## 1, u21, b ## 2);      \
    vcomplex_t u31=u[6];                        \
    accum_cmult_conj(v ## 2, u22, a ## 2);      \
    accum_cmult_conj(w ## 2, u22, b ## 2);      \
    vcomplex_t u32=u[7];                        \
    accum_cmult_conj(v ## 3, u23, a ## 2);      \
    accum_cmult_conj(w ## 3, u23, b ## 2);      \
    vcomplex_t u33=u[8];                        \
    accum_cmult_conj(v ## 1, u31, a ## 3);      \
    accum_cmult_conj(w ## 1, u31, b ## 3);      \
    accum_cmult_conj(v ## 2, u32, a ## 3);      \
    accum_cmult_conj(w ## 2, u32, b ## 3);      \
    accum_cmult_conj(v ## 3, u33, a ## 3);      \
    accum_cmult_conj(w ## 3, u33, b ## 3);      \
  }

// v_i_sp1 = conj(u_{ji})*c_j_sp1,
// w_i_sp2 = conj(u_{ji})*c_j_sp2,
#define MULT_DAG_2SPIN_SP(v,u,c) {                              \
    vcomplex_t u11=u[0];                                        \
    vcomplex_t u12=u[1];                                        \
    vcomplex_t u13=u[2];                                        \
    set_cmult_conj(v ## 1 ## _sp1, u11, c ## 1 ## _sp1);        \
    set_cmult_conj(v ## 1 ## _sp2, u11, c ## 1 ## _sp2);        \
    vcomplex_t u21=u[3];                                        \
    set_cmult_conj(v ## 2 ## _sp1, u12, c ## 1 ## _sp1);        \
    set_cmult_conj(v ## 2 ## _sp2, u12, c ## 1 ## _sp2);        \
    vcomplex_t u22=u[4];                                        \
    set_cmult_conj(v ## 3 ## _sp1, u13, c ## 1 ## _sp1);        \
    set_cmult_conj(v ## 3 ## _sp2, u13, c ## 1 ## _sp2);        \
    vcomplex_t u23=u[5];                                        \
    accum_cmult_conj(v ## 1 ## _sp1, u21, c ## 2 ## _sp1);      \
    accum_cmult_conj(v ## 1 ## _sp2, u21, c ## 2 ## _sp2);      \
    vcomplex_t u31=u[6];                                        \
    accum_cmult_conj(v ## 2 ## _sp1, u22, c ## 2 ## _sp1);      \
    accum_cmult_conj(v ## 2 ## _sp2, u22, c ## 2 ## _sp2);      \
    vcomplex_t u32=u[7];                                        \
    accum_cmult_conj(v ## 3 ## _sp1, u23, c ## 2 ## _sp1);      \
    accum_cmult_conj(v ## 3 ## _sp2, u23, c ## 2 ## _sp2);      \
    vcomplex_t u33=u[8];                                        \
    accum_cmult_conj(v ## 1 ## _sp1, u31, c ## 3 ## _sp1);      \
    accum_cmult_conj(v ## 1 ## _sp2, u31, c ## 3 ## _sp2);      \
    accum_cmult_conj(v ## 2 ## _sp1, u32, c ## 3 ## _sp1);      \
    accum_cmult_conj(v ## 2 ## _sp2, u32, c ## 3 ## _sp2);      \
    accum_cmult_conj(v ## 3 ## _sp1, u33, c ## 3 ## _sp1);      \
    accum_cmult_conj(v ## 3 ## _sp2, u33, c ## 3 ## _sp2);      \
  }
