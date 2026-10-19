/*
        @file    $Id: multi_hop_tm-inc.h #$

        @brief    Macros for hopping in -y direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

// note the gauge field must be shifted in advance, if the shift is needed
#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_TM(psi) SPIN_PROJ1_TP_(psi)
#define SPIN_PROJ2_TM(psi) SPIN_PROJ2_TP_(psi)
#define ACCUM_SPIN1_TM(wc, b, sp) ACCUM_SPIN1_TP_(wc, b, sp)
#define ACCUM_SPIN2_TM(wc, b, sp) ACCUM_SPIN2_TP_(wc, b, sp)
#define ACCUM_SPIN3_TM(wc, b, sp) ACCUM_SPIN3_TP_(wc, b, sp)
#define ACCUM_SPIN4_TM(wc, b, sp) ACCUM_SPIN4_TP_(wc, b, sp)
#define SET_SPIN1_TM(wc, b, sp) SET_SPIN1_TP_(wc, b, sp)
#define SET_SPIN2_TM(wc, b, sp) SET_SPIN2_TP_(wc, b, sp)
#define SET_SPIN3_TM(wc, b, sp) SET_SPIN3_TP_(wc, b, sp)
#define SET_SPIN4_TM(wc, b, sp) SET_SPIN4_TP_(wc, b, sp)

#else

#define SPIN_PROJ1_TM(psi) SPIN_PROJ1_TM_(psi)
#define SPIN_PROJ2_TM(psi) SPIN_PROJ2_TM_(psi)
#define ACCUM_SPIN1_TM(wc, b, sp) ACCUM_SPIN1_TM_(wc, b, sp)
#define ACCUM_SPIN2_TM(wc, b, sp) ACCUM_SPIN2_TM_(wc, b, sp)
#define ACCUM_SPIN3_TM(wc, b, sp) ACCUM_SPIN3_TM_(wc, b, sp)
#define ACCUM_SPIN4_TM(wc, b, sp) ACCUM_SPIN4_TM_(wc, b, sp)
#define SET_SPIN1_TM(wc, b, sp) SET_SPIN1_TM_(wc, b, sp)
#define SET_SPIN2_TM(wc, b, sp) SET_SPIN2_TM_(wc, b, sp)
#define SET_SPIN3_TM(wc, b, sp) SET_SPIN3_TM_(wc, b, sp)
#define SET_SPIN4_TM(wc, b, sp) SET_SPIN4_TM_(wc, b, sp)

#endif

#define HOP_TM_A(SET_OR_ACCUM, a )		\
  {						\
    vcomplex_t sp1=pw[     0 + a];		\
    vcomplex_t sp2=pw[  NCOL + a];		\
    vcomplex_t sp3=pw[2*NCOL + a];		\
    vcomplex_t sp4=pw[3*NCOL + a];		\
    SPIN_PROJ1_TM(sp);				\
    SPIN_PROJ2_TM(sp);				\
    shift_tm(sp1, need_shift);			\
    shift_tm(sp2, need_shift);			\
    SET_OR_ACCUM ## _cmult_conj(v0_sp1, u ## a ## 0, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp1, u ## a ## 1, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp1, u ## a ## 2, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v0_sp2, u ## a ## 0, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp2, u ## a ## 1, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp2, u ## a ## 2, sp2);  \
  }



#define MULT_TM_A( SET_OR_ACCUM, a )		\
  {						\
    vcomplex_t sp1=pw[     0 + a];		\
    vcomplex_t sp2=pw[  NCOL + a];		\
    vcomplex_t sp3=pw[2*NCOL + a];		\
    vcomplex_t sp4=pw[3*NCOL + a];		\
    SPIN_PROJ1_TM(sp);				\
    SPIN_PROJ2_TM(sp);				\
    SET_OR_ACCUM ## _cmult_conj(v0_sp1, u ## a ## 0, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp1, u ## a ## 1, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp1, u ## a ## 2, sp1);  \
    SET_OR_ACCUM ## _cmult_conj(v0_sp2, u ## a ## 0, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v1_sp2, u ## a ## 1, sp2);  \
    SET_OR_ACCUM ## _cmult_conj(v2_sp2, u ## a ## 2, sp2);  \
  }


#define HOP_TM_A_READ_BUFFER(tmp, a)	{	     \
    unpack_tm(tmp ## _sp1, tmp ## _sp2, rbuf, 0, a); \
    shift_tm(tmp ## _sp1);			     \
    shift_tm(tmp ## _sp2);			     \
  }


{
  constexpr int NCOL=3;
  assert(NCOL==m_Nc);

#ifdef AFOPR_HAS_SIMD_T
  //const bool need_simd_shift = has_simd_y && (ry==0);
  //SET_MASK(need_shift, need_simd_shift);
  SET_MASK(need_shift, at_boundary_tm);
#else
  SET_MASK(need_shift, false);
#endif
  
  vcomplex_t v0_sp1;
  vcomplex_t v1_sp1;
  vcomplex_t v2_sp1;
  vcomplex_t v0_sp2;
  vcomplex_t v1_sp2;
  vcomplex_t v2_sp2;

#ifdef NEED_TO_READ_BUFFER
#ifdef AFOPR_HAS_SIMD_T
  MULT_TM_A(set,   0);  // vb_sp? = Ub0 sp?
  MULT_TM_A(accum, 1);  // vb_sp? += Ub1 sp?
  MULT_TM_A(accum, 2);  // vb_sp? += Ub2 sp?
#endif
  HOP_TM_A_READ_BUFFER(v0, 0);
  HOP_TM_A_READ_BUFFER(v1, 1);
  HOP_TM_A_READ_BUFFER(v2, 2);
#else
  HOP_TM_A(set,   0);
  HOP_TM_A(accum, 1);
  HOP_TM_A(accum, 2);
#endif

#ifdef PACK_TO_BUFFER
  pack_tm(pv, v0_sp1, v0_sp2, 0, 0);
  pack_tm(pv, v1_sp1, v1_sp2, 0, 1);
  pack_tm(pv, v2_sp1, v2_sp2, 0, 2);
#else

#ifdef HOP_DOES_NOT_ACCUM
  DEFINE_SPIN1(wc);
  DEFINE_SPIN2(wc);
  DEFINE_SPIN3(wc);
  DEFINE_SPIN4(wc);
  SET_SPIN1_TM(wc, 0, v);
  SET_SPIN1_TM(wc, 1, v);
  SET_SPIN1_TM(wc, 2, v);
  PUT_SPIN1(pv, wc);
  SET_SPIN2_TM(wc, 0, v);
  SET_SPIN2_TM(wc, 1, v);
  SET_SPIN2_TM(wc, 2, v);
  PUT_SPIN2(pv, wc);
  SET_SPIN3_TM(wc, 0, v);
  SET_SPIN3_TM(wc, 1, v);
  SET_SPIN3_TM(wc, 2, v);
  PUT_SPIN3(pv, wc);
  SET_SPIN4_TM(wc, 0, v);
  SET_SPIN4_TM(wc, 1, v);
  SET_SPIN4_TM(wc, 2, v);
  PUT_SPIN4(pv, wc);
#else
  DEFINE_AND_SET_SPIN1(wc, pv);
  DEFINE_AND_SET_SPIN2(wc, pv);
  ACCUM_SPIN1_TM(wc, 0, v);
  ACCUM_SPIN1_TM(wc, 1, v);
  ACCUM_SPIN1_TM(wc, 2, v);
  PUT_SPIN1(pv, wc);
  DEFINE_AND_SET_SPIN3(wc, pv);
  ACCUM_SPIN2_TM(wc, 0, v);
  ACCUM_SPIN2_TM(wc, 1, v);
  ACCUM_SPIN2_TM(wc, 2, v);
  PUT_SPIN2(pv, wc);
  DEFINE_AND_SET_SPIN4(wc, pv);
  ACCUM_SPIN3_TM(wc, 0, v);
  ACCUM_SPIN3_TM(wc, 1, v);
  ACCUM_SPIN3_TM(wc, 2, v);
  PUT_SPIN3(pv, wc);
  ACCUM_SPIN4_TM(wc, 0, v);
  ACCUM_SPIN4_TM(wc, 1, v);
  ACCUM_SPIN4_TM(wc, 2, v);
  PUT_SPIN4(pv, wc);
#endif

#endif
}

#undef MULT_TM_A
#undef HOP_TM_A
#undef HOP_TM_A_READ_BUFFER
#undef SPIN_PROJ1_TM
#undef SPIN_PROJ2_TM
#undef ACCUM_SPIN1_TM
#undef ACCUM_SPIN2_TM
#undef ACCUM_SPIN3_TM
#undef ACCUM_SPIN4_TM
#undef SET_SPIN1_TM
#undef SET_SPIN2_TM
#undef SET_SPIN3_TM
#undef SET_SPIN4_TM
