/* 
        @file    $Id: spinproj_macro-inc.h #$

        @brief    Macros for  Spin Projectoin/Reconstruction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

/////////////////////////////////////////////////
// +x
/////////////////////////////////////////////////

// w1 = psi_1+ i psi_4
// w2 = psi_2+ i psi_3
#define GET_SPIN_PROJ1_XP_(psi) ( psi ## 1 + TimesI(psi ## 4) )
#define GET_SPIN_PROJ2_XP_(psi) ( psi ## 2 + TimesI(psi ## 3) )
#define SPIN_PROJ1_XP_(psi) ( accum_TimesI(psi ## 1 , psi ## 4 ) )
#define SPIN_PROJ2_XP_(psi) ( accum_TimesI(psi ## 2 , psi ## 3 ) )


// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = -i w2
// psi_4 = -i w1
#define ACCUM_SPIN1_XP_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_XP_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_XP_(wc, b, sp) ( accum_TimesMinusI(wc ## b ## _sp3, sp ## b ## _sp2) )
#define ACCUM_SPIN4_XP_(wc, b, sp) ( accum_TimesMinusI(wc ## b ## _sp4, sp ## b ## _sp1) )
#define SET_SPIN1_XP_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_XP_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_XP_(wc, b, sp) ( wc ## b ## _sp3 = TimesMinusI(sp ## b ## _sp2) )
#define SET_SPIN4_XP_(wc, b, sp) ( wc ## b ## _sp4 = TimesMinusI(sp ## b ## _sp1) )


#define ACCUM1_XP_(b, sp1)                      \
  {                                             \
    accum_add(wc ## b ## _sp1, sp1);            \
    accum_TimesMinusI(wc ## b ## _sp4, sp1);    \
  }

#define ACCUM2_XP_(b, sp2)                      \
  {                                             \
    accum_add(wc ## b ## _sp2, sp2 );           \
    accum_TimesMinusI(wc ## b ## _sp3, sp2);    \
  }

#define MULT_AND_ACCUM1_XP(b, uba, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult(usp, uba, spa);                   \
    ACCUM1_XP(b, usp);                          \
  }



#define MULT_AND_ACCUM2_XP(b, uba, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult(usp, uba, spa);                   \
    ACCUM2_XP(b, usp);                          \
  }


/////////////////////////////////////////////////
// -x
/////////////////////////////////////////////////

// w1 = psi_1- i psi_4
// w2 = psi_2- i psi_3
#define GET_SPIN_PROJ1_XM_(psi) ( psi ## 1 + TimesMinusI(psi ## 4) )
#define GET_SPIN_PROJ2_XM_(psi) ( psi ## 2 + TimesMinusI(psi ## 3) )
#define SPIN_PROJ1_XM_(psi) ( accum_TimesMinusI( psi ## 1, psi ## 4) )
#define SPIN_PROJ2_XM_(psi) ( accum_TimesMinusI( psi ## 2, psi ## 3) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = i w2
// psi_4 = i w1
#define ACCUM_SPIN1_XM_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_XM_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_XM_(wc, b, sp) ( accum_TimesI(wc ## b ## _sp3, sp ## b ## _sp2) )
#define ACCUM_SPIN4_XM_(_wc, b, sp) ( accum_TimesI(wc ## b ## _sp4, sp ## b ## _sp1) )
#define SET_SPIN1_XM_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_XM_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_XM_(wc, b, sp) ( wc ## b ## _sp3 = TimesI(sp ## b ## _sp2) )
#define SET_SPIN4_XM_(_wc, b, sp) (wc ## b ## _sp4 = TimesI(sp ## b ## _sp1) )


#define ACCUM1_XM_(b,usp)               \
  {                                     \
    accum_add(wc ## b ## _sp1, usp);    \
    accum_TimesI(wc ## b ## _sp4, usp); \
  }

#define ACCUM2_XM_(b,usp)               \
  {                                     \
    accum_add(wc ## b ## _sp2, usp);    \
    accum_TimesI(wc ## b ## _sp3, usp); \
  }


#define MULT_AND_ACCUM1_XM(b, uab, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult_conj(usp, uab, spa);              \
    ACCUM1_XM(b,usp)                            \
  }


#define MULT_AND_ACCUM2_XM(b, uab, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult_conj(usp, uab, spa);              \
    ACCUM2_XM(b,usp)                            \
  }

//#define MULT_AND_ACCUM2_XM_(b,sp) MULT2_AND_ACCUM2_XM(b, b, sp)



/////////////////////////////////////////////////
// +y
/////////////////////////////////////////////////

// spin projectoin
// w1 = psi_1+ psi_4
// w2 = psi_2- psi_3
#define GET_SPIN_PROJ1_YP_(psi) ( (psi ## 1) + (psi ## 4) )
#define GET_SPIN_PROJ2_YP_(psi) ( (psi ## 2) - (psi ## 3) )
#define SPIN_PROJ1_YP_(psi) ( accum_add(psi ## 1, psi ## 4) )
#define SPIN_PROJ2_YP_(psi) ( accum_sub(psi ## 2, psi ## 3) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = - w2
// psi_4 = + w1
#define ACCUM_SPIN1_YP_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_YP_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_YP_(wc, b, sp) ( accum_sub(wc ## b ## _sp3, sp ## b ## _sp2) )
#define ACCUM_SPIN4_YP_(wc, b, sp) ( accum_add(wc ## b ## _sp4, sp ## b ## _sp1) )
#define SET_SPIN1_YP_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_YP_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_YP_(wc, b, sp) { set_zero(wc ## b ## _sp3) ; accum_sub(wc ## b ## _sp3, sp ## b ## _sp2); }
#define SET_SPIN4_YP_(wc, b, sp) ( wc ## b ## _sp4 = sp ## b ## _sp1 )

#define ACCUM1_YP_(b, sp1)              \
  {                                     \
    accum_add(wc ## b ## _sp1, sp1);    \
    accum_add(wc ## b ## _sp4, sp1);    \
  }

#define ACCUM2_YP_(b, sp2)              \
  {                                     \
    accum_add(wc ## b ## _sp2, sp2);    \
    accum_sub(wc ## b ## _sp3, sp2);    \
  }


#define MULT_AND_ACCUM1_YP(b, uba, spa) \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult(usp, uba, spa);           \
    ACCUM1_YP(b, usp);                  \
  }

#define MULT_AND_ACCUM2_YP(b, uba, spa) \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult(usp, uba, spa);           \
    ACCUM2_YP(b, usp);                  \
  }



/////////////////////////////////////////////////
// -y
/////////////////////////////////////////////////

// spin projectoin
// w1 = psi_1- psi_4
// w2 = psi_2+ psi_3
#define GET_SPIN_PROJ1_YM_(psi) ( (psi ## 1) - (psi ## 4) )
#define GET_SPIN_PROJ2_YM_(psi) ( (psi ## 2) + (psi ## 3) )
#define SPIN_PROJ1_YM_(psi) ( accum_sub(psi ## 1, psi ## 4) )
#define SPIN_PROJ2_YM_(psi) ( accum_add(psi ## 2, psi ## 3) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = w2
// psi_4 = -w1
#define ACCUM_SPIN1_YM_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_YM_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_YM_(wc, b, sp) ( accum_add(wc ## b ## _sp3, sp ## b ## _sp2) )
#define ACCUM_SPIN4_YM_(wc, b, sp) ( accum_sub(wc ## b ## _sp4, sp ## b ## _sp1) )
#define SET_SPIN1_YM_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_YM_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_YM_(wc, b, sp) ( wc ## b ## _sp3 = sp ## b ## _sp2 )
#define SET_SPIN4_YM_(wc, b, sp) { set_zero(wc ## b ## _sp4); accum_sub(wc ## b ## _sp4, sp ## b ## _sp1); }

#define ACCUM1_YM_(b,sp1)               \
  {                                     \
    accum_add(wc ## b ## _sp1, sp1);    \
    accum_sub(wc ## b ## _sp4, sp1);    \
  }

#define ACCUM2_YM_(b,sp2)               \
  {                                     \
    accum_add(wc ## b ## _sp2, sp2);    \
    accum_add(wc ## b ## _sp3, sp2);    \
  }

#define MULT_AND_ACCUM1_YM(b,uab,spa)   \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult_conj(usp, uab, spa);      \
    ACCUM1_YM(b,usp)                    \
  }

#define MULT_AND_ACCUM2_YM(b,uab,spa)   \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult_conj(usp, uab, spa);      \
    ACCUM2_YM(b,usp)                    \
  }


/////////////////////////////////////////////////
// +z
/////////////////////////////////////////////////

// spin projectoin
// w1 = psi_1+ i psi_3
// w2 = psi_2- i psi_4
#define GET_SPIN_PROJ1_ZP_(psi) ( (psi ## 1) + TimesI(psi ## 3) )
#define GET_SPIN_PROJ2_ZP_(psi) ( (psi ## 2) + TimesMinusI(psi ## 4) )
#define SPIN_PROJ1_ZP_(psi) ( accum_TimesI(psi ## 1, psi ## 3) )
#define SPIN_PROJ2_ZP_(psi) ( accum_TimesMinusI(psi ## 2, psi ## 4) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = - i w1
// psi_4 = + i w2
#define ACCUM_SPIN1_ZP_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_ZP_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_ZP_(wc, b, sp) ( accum_TimesMinusI(wc ## b ## _sp3, sp ## b ## _sp1) )
#define ACCUM_SPIN4_ZP_(wc, b, sp) ( accum_TimesI(wc ## b ## _sp4, sp ## b ## _sp2) )
#define SET_SPIN1_ZP_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_ZP_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_ZP_(wc, b, sp) ( wc ## b ## _sp3 = TimesMinusI(sp ## b ## _sp1) )
#define SET_SPIN4_ZP_(wc, b, sp) ( wc ## b ## _sp4 = TimesI(sp ## b ## _sp2) )


#define ACCUM1_ZP_(b, sp1)                      \
  {                                             \
    accum_add(wc ## b ## _sp1, sp1);            \
    accum_TimesMinusI(wc ## b ## _sp3, sp1);    \
  }

#define ACCUM2_ZP_(b, sp2)                      \
  {                                             \
    accum_add(wc ## b ## _sp2, sp2);            \
    accum_TimesI(wc ## b ## _sp4, sp2);         \
  }

#define MULT_AND_ACCUM1_ZP(b, uba, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult(usp, uba, spa);                   \
    ACCUM1_ZP(b, usp);                          \
  }

#define MULT_AND_ACCUM2_ZP(b, uba, spa)         \
  {                                             \
    vcomplex_t usp;                             \
    set_cmult(usp, uba, spa);                   \
    ACCUM2_ZP(b, usp);                          \
  }


/////////////////////////////////////////////////
// -z
/////////////////////////////////////////////////

// spin projectoin
// w1 = psi_1- i psi_3
// w2 = psi_2+ i psi_4
#define GET_SPIN_PROJ1_ZM_(psi) ( (psi ## 1) + TimesMinusI(psi ## 3) )
#define GET_SPIN_PROJ2_ZM_(psi) ( (psi ## 2) + TimesI(psi ## 4) )
#define SPIN_PROJ1_ZM_(psi) ( accum_TimesMinusI(psi ## 1, psi ## 3) )
#define SPIN_PROJ2_ZM_(psi) ( accum_TimesI(psi ## 2, psi ## 4) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = i w1
// psi_4 = -i w2
#define ACCUM_SPIN1_ZM_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_ZM_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_ZM_(wc, b, sp) ( accum_TimesI(wc ## b ## _sp3, sp ## b ## _sp1) )
#define ACCUM_SPIN4_ZM_(wc, b, sp) ( accum_TimesMinusI(wc ## b ## _sp4, sp ## b ## _sp2) )
#define SET_SPIN1_ZM_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_ZM_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_ZM_(wc, b, sp) ( wc ## b ## _sp3 = TimesI(sp ## b ## _sp1) )
#define SET_SPIN4_ZM_(wc, b, sp) ( wc ## b ## _sp4 = TimesMinusI(sp ## b ## _sp2) )

#define ACCUM1_ZM_(b,sp1)               \
  {                                     \
    accum_add(wc ## b ## _sp1, sp1);    \
    accum_TimesI(wc ## b ## _sp3, sp1); \
  }

#define ACCUM2_ZM_(b,sp2)                       \
  {                                             \
    accum_add(wc ## b ## _sp2, sp2);            \
    accum_TimesMinusI(wc ## b ## _sp4, sp2);    \
  }

#define MULT_AND_ACCUM1_ZM(b,uab,spa)   \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult_conj(usp, uab, spa);      \
    ACCUM1_ZM(b,usp)                    \
  }

#define MULT_AND_ACCUM2_ZM(b,uab,spa)   \
  {                                     \
    vcomplex_t usp;                     \
    set_cmult_conj(usp, uab, spa);      \
    ACCUM2_ZM(b,usp)                    \
  }

// todo: refactor the treatment of the chiral rep.
#ifdef SPINPROJ_USE_CHIRAL_REP
#include "spinproj_macro_chiral-inc.h"
#else
#include "spinproj_macro_dirac-inc.h"
#endif

/////////////////////////////////////////////////
// general helper
/////////////////////////////////////////////////

#define DEFINE_AND_SET_SPINOR(c, w)     \
  vcomplex_t c ## 1_sp1 =w[0];          \
  vcomplex_t c ## 2_sp1 =w[1];          \
  vcomplex_t c ## 3_sp1 =w[2];          \
  vcomplex_t c ## 1_sp2 =w[3];          \
  vcomplex_t c ## 2_sp2 =w[4];          \
  vcomplex_t c ## 3_sp2 =w[5];          \
  vcomplex_t c ## 1_sp3 =w[6];          \
  vcomplex_t c ## 2_sp3 =w[7];          \
  vcomplex_t c ## 3_sp3 =w[8];          \
  vcomplex_t c ## 1_sp4 =w[9];          \
  vcomplex_t c ## 2_sp4 =w[10];         \
  vcomplex_t c ## 3_sp4 =w[11];


#define DEFINE_AND_SET_SPINOR_C1(c, w)  \
  vcomplex_t c ## 1_sp1 =w[0];          \
  vcomplex_t c ## 1_sp2 =w[3];          \
  vcomplex_t c ## 1_sp3 =w[6];          \
  vcomplex_t c ## 1_sp4 =w[9];

#define DEFINE_AND_SET_SPINOR_C2(c, w)  \
  vcomplex_t c ## 2_sp1 =w[1];          \
  vcomplex_t c ## 2_sp2 =w[4];          \
  vcomplex_t c ## 2_sp3 =w[7];          \
  vcomplex_t c ## 2_sp4 =w[10];         \

#define DEFINE_AND_SET_SPINOR_C3(c, w)  \
  vcomplex_t c ## 3_sp1 =w[2];          \
  vcomplex_t c ## 3_sp2 =w[5];          \
  vcomplex_t c ## 3_sp3 =w[8];          \
  vcomplex_t c ## 3_sp4 =w[11];


#define DEFINE_AND_SET_SPIN1(c, w)      \
  vcomplex_t c ## 0_sp1 =w[0];          \
  vcomplex_t c ## 1_sp1 =w[1];          \
  vcomplex_t c ## 2_sp1 =w[2];

#define DEFINE_AND_SET_SPIN2(c, w)      \
  vcomplex_t c ## 0_sp2 =w[3];          \
  vcomplex_t c ## 1_sp2 =w[4];          \
  vcomplex_t c ## 2_sp2 =w[5];

#define DEFINE_AND_SET_SPIN3(c, w)      \
  vcomplex_t c ## 0_sp3 =w[6];          \
  vcomplex_t c ## 1_sp3 =w[7];          \
  vcomplex_t c ## 2_sp3 =w[8];

#define DEFINE_AND_SET_SPIN4(c, w)      \
  vcomplex_t c ## 0_sp4 =w[9];          \
  vcomplex_t c ## 1_sp4 =w[10];         \
  vcomplex_t c ## 2_sp4 =w[11];

#define DEFINE_SPIN1(c)  vcomplex_t c ## 0_sp1, c ## 1_sp1, c ## 2_sp1
#define DEFINE_SPIN2(c)  vcomplex_t c ## 0_sp2, c ## 1_sp2, c ## 2_sp2
#define DEFINE_SPIN3(c)  vcomplex_t c ## 0_sp3, c ## 1_sp3, c ## 2_sp3
#define DEFINE_SPIN4(c)  vcomplex_t c ## 0_sp4, c ## 1_sp4, c ## 2_sp4

#define PUT_SPIN1(w, c) \
  w[0] = c ## 0_sp1;    \
  w[1] = c ## 1_sp1;    \
  w[2] = c ## 2_sp1;

#define PUT_SPIN2(w, c) \
  w[3] = c ## 0_sp2;    \
  w[4] = c ## 1_sp2;    \
  w[5] = c ## 2_sp2;

#define PUT_SPIN3(w, c) \
  w[6] = c ## 0_sp3;    \
  w[7] = c ## 1_sp3;    \
  w[8] = c ## 2_sp3;

#define PUT_SPIN4(w, c) \
  w[9]  = c ## 0_sp4;   \
  w[10] = c ## 1_sp4;   \
  w[11] = c ## 2_sp4;


#define MULT_CONJ_AND_SET(tmp, b,sp)                    \
    set_cmult_conj(tmp ## b ## _ ## sp, u ## b, sp)

#define MULT_CONJ_AND_ACCUM(tmp, b,sp)                  \
    accum_cmult_conj(tmp ## b ## _ ## sp, u ## b, sp)

#define MULT_CONJ_AND_ACCUM(tmp, b,sp)                  \
    accum_cmult_conj(tmp ## b ## _ ## sp, u ## b, sp)
