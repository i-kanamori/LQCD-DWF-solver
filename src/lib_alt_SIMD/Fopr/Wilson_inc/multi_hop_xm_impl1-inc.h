/* 
        @file    $Id: multi_hop_xm_impl-inc.h #$

        @brief    Macros for hopping in -x direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_XM(psi) SPIN_PROJ1_XP_(psi) 
#define SPIN_PROJ2_XM(psi) SPIN_PROJ2_XP_(psi) 
#define ACCUM_SPIN1_XM(wc, b, sp) ACCUM_SPIN1_XP_(wc, b, sp)
#define ACCUM_SPIN2_XM(wc, b, sp) ACCUM_SPIN2_XP_(wc, b, sp)
#define ACCUM_SPIN3_XM(wc, b, sp) ACCUM_SPIN3_XP_(wc, b, sp)
#define ACCUM_SPIN4_XM(wc, b, sp) ACCUM_SPIN4_XP_(wc, b, sp)

#else

#define SPIN_PROJ1_XM(psi) SPIN_PROJ1_XM_(psi) 
#define SPIN_PROJ2_XM(psi) SPIN_PROJ2_XM_(psi) 
#define ACCUM_SPIN1_XM(wc, b, sp) ACCUM_SPIN1_XM_(wc, b, sp)
#define ACCUM_SPIN2_XM(wc, b, sp) ACCUM_SPIN2_XM_(wc, b, sp)
#define ACCUM_SPIN3_XM(wc, b, sp) ACCUM_SPIN3_XM_(wc, b, sp)
#define ACCUM_SPIN4_XM(wc, b, sp) ACCUM_SPIN4_XM_(wc, b, sp)

#endif

#ifdef USE_SAP
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_SAP(sp, a)

#elif defined(USE_BLOCK_HOP)
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_BLOCK_HOP(sp, a)

#else
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_(sp, a)
#endif

#define PREPARE_SPINOR_(sp, a)		\
    vcomplex_t sp ## 1=pw0[     0 + a]; \
    vcomplex_t sp ## 2=pw0[  NCOL + a]; \
    vcomplex_t sp ## 3=pw0[2*NCOL + a]; \
    vcomplex_t nsp1=pw[     0 + a];	\
    vcomplex_t nsp2=pw[  NCOL + a];	\
    vcomplex_t nsp3=pw[2*NCOL + a];	\
    shift_and_merge_xm(sp ## 1,nsp1);   \
    vcomplex_t sp ## 4=pw0[3*NCOL + a]; \
    vcomplex_t nsp4=pw[3*NCOL + a];	\
    shift_and_merge_xm(sp ## 2,nsp2);   \
    shift_and_merge_xm(sp ## 3,nsp3);	\
    shift_and_merge_xm(sp ## 4,nsp4);   

#define PREPARE_SPINOR_SAP(sp, a)	\
    vcomplex_t sp ## 1=pw0[     0 + a]; \
    vcomplex_t sp ## 2=pw0[  NCOL + a]; \
    vcomplex_t sp ## 3=pw0[2*NCOL + a]; \
    vcomplex_t sp ## 4=pw0[3*NCOL + a]; \
    vcomplex_t vzero;			\
    set_zero(vzero);			\
    shift_and_merge_xm(sp ## 1,vzero);  \
    shift_and_merge_xm(sp ## 2,vzero);  \
    shift_and_merge_xm(sp ## 3,vzero);  \
    shift_and_merge_xm(sp ## 4,vzero); 

#define PREPARE_SPINOR_BLOCK_HOP(sp, a)	\
  vcomplex_t sp ## 1;			\
  vcomplex_t sp ## 2;			\
  vcomplex_t sp ## 3;			\
  vcomplex_t sp ## 4;			\
  set_zero(sp ## 1);			\
  set_zero(sp ## 2);			\
  set_zero(sp ## 3);			\
  set_zero(sp ## 4);			\
  vcomplex_t nsp1=pw[     0 + a];	\
  vcomplex_t nsp2=pw[  NCOL + a];	\
  vcomplex_t nsp3=pw[2*NCOL + a];	\
  vcomplex_t nsp4=pw[3*NCOL + a];	\
  shift_and_merge_xm(sp ## 1,nsp1);	\
  shift_and_merge_xm(sp ## 2,nsp2);	\
  shift_and_merge_xm(sp ## 3,nsp3);	\
  shift_and_merge_xm(sp ## 4,nsp4);   

#define HOP_XM_A( SET_OR_ACCUM, a ) \
  {				    \
    PREPARE_SPINOR(sp, a);	    \
    SPIN_PROJ1_XM(sp);		    \
    SPIN_PROJ2_XM(sp);		    \
    SET_OR_ACCUM ## _cmult_conj(v0_sp1, u ## a ## 0, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp1, u ## a ## 1, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp1, u ## a ## 2, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v0_sp2, u ## a ## 0, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp2, u ## a ## 1, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp2, u ## a ## 2, sp2);  \
  }

{
  constexpr int NCOL=3;
  assert(NCOL==m_Nc);

  vcomplex_t v0_sp1;
  vcomplex_t v1_sp1;
  vcomplex_t v2_sp1;
  vcomplex_t v0_sp2;
  vcomplex_t v1_sp2;
  vcomplex_t v2_sp2;
  HOP_XM_A(set,   0); // va = u_0a^* in0
  HOP_XM_A(accum, 1); // va += u_1a^* in1
  HOP_XM_A(accum, 2); // va += u_2a^* in2

  DEFINE_AND_SET_SPIN1(wc, pv);
  DEFINE_AND_SET_SPIN2(wc, pv);
  ACCUM_SPIN1_XM(wc, 0, v);
  ACCUM_SPIN1_XM(wc, 1, v);
  ACCUM_SPIN1_XM(wc, 2, v);
  PUT_SPIN1(pv, wc);
  DEFINE_AND_SET_SPIN3(wc, pv);
  ACCUM_SPIN2_XM(wc, 0, v);
  ACCUM_SPIN2_XM(wc, 1, v);
  ACCUM_SPIN2_XM(wc, 2, v);
  PUT_SPIN2(pv, wc);
  DEFINE_AND_SET_SPIN4(wc, pv);
  ACCUM_SPIN3_XM(wc, 0, v);
  ACCUM_SPIN3_XM(wc, 1, v);
  ACCUM_SPIN3_XM(wc, 2, v);
  PUT_SPIN3(pv, wc);
  ACCUM_SPIN4_XM(wc, 0, v);
  ACCUM_SPIN4_XM(wc, 1, v);
  ACCUM_SPIN4_XM(wc, 2, v);
  PUT_SPIN4(pv, wc);

}


#undef HOP_XM_A
#undef SPIN_PROJ1_XM
#undef SPIN_PROJ2_XM
#undef ACCUM_SPIN1_XM
#undef ACCUM_SPIN2_XM
#undef ACCUM_SPIN3_XM
#undef ACCUM_SPIN4_XM
#undef PREPARE_SPINOR
#undef PREPARE_SPINOR_
#undef PREPARE_SPINOR_SAP
#undef PREPARE_SPINOR_BLOCK_HOP

