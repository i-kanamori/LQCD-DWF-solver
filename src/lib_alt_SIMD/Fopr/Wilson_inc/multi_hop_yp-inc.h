/* 
        @file    $Id: multi_hop_yp-inc.h #$

        @brief    Macros for hopping in +y direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/


#ifdef SPINPROJ_USE_DAGGER
#define SPIN_PROJ1_YP(psi) SPIN_PROJ1_YM_(psi)
#define SPIN_PROJ2_YP(psi) SPIN_PROJ2_YM_(psi)
#define ACCUM_SPIN1_YP(wc, b, sp) ACCUM_SPIN1_YM_(wc, b, sp)
#define ACCUM_SPIN2_YP(wc, b, sp) ACCUM_SPIN2_YM_(wc, b, sp)
#define ACCUM_SPIN3_YP(wc, b, sp) ACCUM_SPIN3_YM_(wc, b, sp)
#define ACCUM_SPIN4_YP(wc, b, sp) ACCUM_SPIN4_YM_(wc, b, sp)
#define SET_SPIN1_YP(wc, b, sp) SET_SPIN1_YM_(wc, b, sp)
#define SET_SPIN2_YP(wc, b, sp) SET_SPIN2_YM_(wc, b, sp)
#define SET_SPIN3_YP(wc, b, sp) SET_SPIN3_YM_(wc, b, sp)
#define SET_SPIN4_YP(wc, b, sp) SET_SPIN4_YM_(wc, b, sp)

#else
#define SPIN_PROJ1_YP(psi) SPIN_PROJ1_YP_(psi)
#define SPIN_PROJ2_YP(psi) SPIN_PROJ2_YP_(psi)
#define ACCUM_SPIN1_YP(wc, b, sp) ACCUM_SPIN1_YP_(wc, b, sp)
#define ACCUM_SPIN2_YP(wc, b, sp) ACCUM_SPIN2_YP_(wc, b, sp)
#define ACCUM_SPIN3_YP(wc, b, sp) ACCUM_SPIN3_YP_(wc, b, sp)
#define ACCUM_SPIN4_YP(wc, b, sp) ACCUM_SPIN4_YP_(wc, b, sp)
#define SET_SPIN1_YP(wc, b, sp) SET_SPIN1_YP_(wc, b, sp)
#define SET_SPIN2_YP(wc, b, sp) SET_SPIN2_YP_(wc, b, sp)
#define SET_SPIN3_YP(wc, b, sp) SET_SPIN3_YP_(wc, b, sp)
#define SET_SPIN4_YP(wc, b, sp) SET_SPIN4_YP_(wc, b, sp)

#endif

#define HOP_YP_A( SET_OR_ACCUM, a ) \
  {				    \
    vcomplex_t sp1=pw[     0 + a];  \
    vcomplex_t sp2=pw[  NCOL + a];  \
    vcomplex_t sp3=pw[2*NCOL + a];  \
    vcomplex_t sp4=pw[3*NCOL + a];  \
    SPIN_PROJ1_YP(sp);		    \
    SPIN_PROJ2_YP(sp);		    \
    shift_yp(sp1, need_shift);	    \
    shift_yp(sp2, need_shift);	    \
    SET_OR_ACCUM ## _cmult(v0_sp1, u0 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v1_sp1, u1 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v2_sp1, u2 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v0_sp2, u0 ## a , sp2);  \
    SET_OR_ACCUM ## _cmult(v1_sp2, u1 ## a , sp2);  \
    SET_OR_ACCUM ## _cmult(v2_sp2, u2 ## a , sp2);  \
  }


#define HOP_YP_A_READ_BUFFER(SET_OR_ACCUM, a)	    \
  {						    \
    vcomplex_t sp1=pw[     0 + a];		    \
    vcomplex_t sp2=pw[  NCOL + a];		    \
    vcomplex_t sp3=pw[2*NCOL + a];		    \
    vcomplex_t sp4=pw[3*NCOL + a];		    \
    SPIN_PROJ1_YP(sp);				    \
    SPIN_PROJ2_YP(sp);				    \
    unpack_yp(sp1, sp2, rbuf, 0, a);		    \
    shift_yp(sp1);				    \
    shift_yp(sp2);				    \
    SET_OR_ACCUM ## _cmult(v0_sp1, u0 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v1_sp1, u1 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v2_sp1, u2 ## a , sp1);  \
    SET_OR_ACCUM ## _cmult(v0_sp2, u0 ## a , sp2);  \
    SET_OR_ACCUM ## _cmult(v1_sp2, u1 ## a , sp2);  \
    SET_OR_ACCUM ## _cmult(v2_sp2, u2 ## a , sp2);  \
  }

{
  constexpr int NCOL=3;
  assert(NCOL==m_Nc);

#ifdef AFOPR_HAS_SIMD_Y
  //  const bool need_simd_shift = has_simd_y && (ry+1==Ry);
  //  SET_MASK(need_shift, need_simd_shift);
  SET_MASK(need_shift, at_boundary_yp);
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
  HOP_YP_A_READ_BUFFER(set,   0);
  HOP_YP_A_READ_BUFFER(accum, 1);
  HOP_YP_A_READ_BUFFER(accum, 2);
#else
  HOP_YP_A(set, 0);    // va += Ua0 in0
  HOP_YP_A(accum, 1);  // va += Ua1 in1
  HOP_YP_A(accum, 2);  // va += Ua2 in2
#endif

#ifdef HOP_DOES_NOT_ACCUM
  DEFINE_SPIN2(wc);
  DEFINE_SPIN3(wc);
  DEFINE_SPIN4(wc);
  DEFINE_SPIN1(wc);
  SET_SPIN1_YP(wc, 0, v);
  SET_SPIN1_YP(wc, 1, v);
  SET_SPIN1_YP(wc, 2, v);
  PUT_SPIN1(pv, wc);
  SET_SPIN2_YP(wc, 0, v);
  SET_SPIN2_YP(wc, 1, v);
  SET_SPIN2_YP(wc, 2, v);
  PUT_SPIN2(pv, wc);
  SET_SPIN3_YP(wc, 0, v);
  SET_SPIN3_YP(wc, 1, v);
  SET_SPIN3_YP(wc, 2, v);
  PUT_SPIN3(pv, wc);
  SET_SPIN4_YP(wc, 0, v);
  SET_SPIN4_YP(wc, 1, v);
  SET_SPIN4_YP(wc, 2, v);
  PUT_SPIN4(pv, wc);

#else
  DEFINE_AND_SET_SPIN1(wc, pv);
  DEFINE_AND_SET_SPIN2(wc, pv);
  ACCUM_SPIN1_YP(wc, 0, v);
  ACCUM_SPIN1_YP(wc, 1, v);
  ACCUM_SPIN1_YP(wc, 2, v);
  PUT_SPIN1(pv, wc);
  DEFINE_AND_SET_SPIN3(wc, pv);
  ACCUM_SPIN2_YP(wc, 0, v);
  ACCUM_SPIN2_YP(wc, 1, v);
  ACCUM_SPIN2_YP(wc, 2, v);
  PUT_SPIN2(pv, wc);
  DEFINE_AND_SET_SPIN4(wc, pv);
  ACCUM_SPIN3_YP(wc, 0, v);
  ACCUM_SPIN3_YP(wc, 1, v);
  ACCUM_SPIN3_YP(wc, 2, v);
  PUT_SPIN3(pv, wc);
  ACCUM_SPIN4_YP(wc, 0, v);
  ACCUM_SPIN4_YP(wc, 1, v);
  ACCUM_SPIN4_YP(wc, 2, v);
  PUT_SPIN4(pv, wc);
#endif
}

#undef HOP_YP_A
#undef HOP_YP_A_READ_BUFFER
#undef SPIN_PROJ1_YP
#undef SPIN_PROJ2_YP
#undef ACCUM_SPIN1_YP
#undef ACCUM_SPIN2_YP
#undef ACCUM_SPIN3_YP
#undef ACCUM_SPIN4_YP
#undef SET_SPIN1_YP
#undef SET_SPIN2_YP
#undef SET_SPIN3_YP
#undef SET_SPIN4_YP
