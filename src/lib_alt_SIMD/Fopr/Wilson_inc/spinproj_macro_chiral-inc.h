/* 
        @file    $Id: spinproj_macro_chiral-inc.h #$

        @brief    Macros for  Spin Projectoin/Reconstruction (Chiral rep.)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/
/////////////////////////////////////////////////
//
// chiral rep.
//

/////////////////////////////////////////////////
// +t
/////////////////////////////////////////////////
// spin projectoin
// w1 = psi_1+ psi_3
// w2 = psi_2+ psi_4
#define GET_SPIN_PROJ1_TP_(psi) ( (psi ## 1) + (psi ## 3) )
#define GET_SPIN_PROJ2_TP_(psi) ( (psi ## 2) + (psi ## 4) )
#define SPIN_PROJ1_TP_(psi) ( accum_add( (psi ## 1), (psi ## 3) ) )
#define SPIN_PROJ2_TP_(psi) ( accum_add( (psi ## 2), (psi ## 4) ) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = w1
// psi_4 = w2
#define ACCUM_SPIN1_TP_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_TP_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_TP_(wc, b, sp) ( accum_add(wc ## b ## _sp3, sp ## b ## _sp1) )
#define ACCUM_SPIN4_TP_(wc, b, sp) ( accum_add(wc ## b ## _sp4, sp ## b ## _sp2) )
#define SET_SPIN1_TP_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_TP_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_TP_(wc, b, sp) ( wc ## b ## _sp3 = sp ## b ## _sp1 )
#define SET_SPIN4_TP_(wc, b, sp) ( wc ## b ## _sp4 = sp ## b ## _sp2 )

#define ACCUM1_TP_(b, sp1)			\
  {						\
    accum_add(wc ## b ## _sp1, sp1);		\
    accum_add(wc ## b ## _sp3, sp1);		\
  }

#define ACCUM2_TP_(b, sp2)			\
  {						\
    accum_add(wc ## b ## _sp2, sp2);		\
    accum_add(wc ## b ## _sp4, sp2);		\
  }

#define MULT_AND_ACCUM1_TP(b,uba,spa)		\
  {						\
    vcomplex_t usp;				\
    set_cmult(usp, uba, spa);			\
    ACCUM1_TP(b, usp);				\
  }


#define MULT_AND_ACCUM2_TP(b,uba,spa)		\
  {						\
    vcomplex_t usp;				\
    set_cmult(usp, uba, spa);			\
    ACCUM2_TP(b, usp);				\
  }


/////////////////////////////////////////////////
// -t
/////////////////////////////////////////////////

// spin projectoin
// w1 = psi_1- i psi_3
// w2 = psi_2- i psi_4
#define GET_SPIN_PROJ1_TM_(psi) ( (psi ## 1) - (psi ## 3) )
#define GET_SPIN_PROJ2_TM_(psi) ( (psi ## 2) - (psi ## 4) )
#define SPIN_PROJ1_TM_(psi) ( accum_sub( (psi ## 1), (psi ## 3) ) )
#define SPIN_PROJ2_TM_(psi) ( accum_sub( (psi ## 2), (psi ## 4) ) )

// 4 spinor
// psi_1 = w1
// psi_2 = w2
// psi_3 = -w1
// psi_4 = -w2
#define ACCUM_SPIN1_TM_(wc, b, sp) ( accum_add(wc ## b ## _sp1, sp ## b ## _sp1) )
#define ACCUM_SPIN2_TM_(wc, b, sp) ( accum_add(wc ## b ## _sp2, sp ## b ## _sp2) )
#define ACCUM_SPIN3_TM_(wc, b, sp) ( accum_sub(wc ## b ## _sp3, sp ## b ## _sp1) )
#define ACCUM_SPIN4_TM_(wc, b, sp) ( accum_sub(wc ## b ## _sp4, sp ## b ## _sp2) )
#define SET_SPIN1_TM_(wc, b, sp) ( wc ## b ## _sp1 = sp ## b ## _sp1 )
#define SET_SPIN2_TM_(wc, b, sp) ( wc ## b ## _sp2 = sp ## b ## _sp2 )
#define SET_SPIN3_TM_(wc, b, sp) { set_zero( wc ## b ## _sp3 ); accum_sub(wc ## b ## _sp3, sp ## b ## _sp1); }
#define SET_SPIN4_TM_(wc, b, sp) { set_zero( wc ## b ## _sp4 ); accum_sub(wc ## b ## _sp4, sp ## b ## _sp2); }


#define ACCUM1_TM_(b,sp1)		\
  {					\
    accum_add(wc ## b ## _sp1, sp1);	\
    accum_sub(wc ## b ## _sp3, sp1);	\
  }

#define ACCUM2_TM(b,sp2)		\
  {					\
    accum_add(wc ## b ## _sp2, sp2);	\
    accum_sub(wc ## b ## _sp4, sp2);	\
  }

#define MULT_AND_ACCUM1_TM(b,uab,spa)	\
  {					\
    vcomplex_t usp;			\
    set_cmult_conj(usp, uab, spa);	\
    ACCUM1_TM(b,usp)			\
  }


#define MULT_AND_ACCUM2_TM(b,uab,spa)	\
  {					\
    vcomplex_t usp;			\
    set_cmult_conj(usp, uab, spa);	\
    ACCUM2_TM(b,usp)			\
  }

