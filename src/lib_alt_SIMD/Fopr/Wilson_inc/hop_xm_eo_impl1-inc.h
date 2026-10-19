/* 
        @file    $Id: hop_xm_eo_impl1-inc.h #$

        @brief    Macros for hopping in -x direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_XM(psi) SPIN_PROJ1_XP_(psi) 
#define SPIN_PROJ2_XM(psi) SPIN_PROJ2_XP_(psi) 
#define ACCUM1_XM(b,spb) ACCUM1_XP_(b,spb)
#define ACCUM2_XM(b,spb) ACCUM2_XP_(b,spb)

#else

#define SPIN_PROJ1_XM(psi) SPIN_PROJ1_XM_(psi) 
#define SPIN_PROJ2_XM(psi) SPIN_PROJ2_XM_(psi) 
#define ACCUM1_XM(b,spb) ACCUM1_XM_(b,spb)
#define ACCUM2_XM(b,spb) ACCUM2_XM_(b,spb)

#endif
  
#define HOP_XM_A( a )					\
  {							\
    const int offset0=s*NCOL*4+a;			\
    const int offset0_gauge=s*NCOL*NCOL*4+a*NCOL;	\
    const vcomplex_t *u=Uptr+offset0_gauge;		\
    vcomplex_t sp1=pw[offset0];				\
    vcomplex_t sp2=pw[offset0+  NCOL];			\
    vcomplex_t sp3=pw[offset0+2*NCOL];			\
    vcomplex_t sp4=pw[offset0+3*NCOL];			\
    SPIN_PROJ1_XM(sp);					\
    SPIN_PROJ2_XM(sp);					\
    vcomplex_t u0=u[0];					\
    vcomplex_t u1=u[1];					\
    vcomplex_t u2=u[2];					\
    MULT_AND_ACCUM1_XM(0,u0,sp1);			\
    MULT_AND_ACCUM2_XM(0,u0,sp2);			\
    MULT_AND_ACCUM1_XM(1,u1,sp1);			\
    MULT_AND_ACCUM2_XM(1,u1,sp2);			\
    MULT_AND_ACCUM1_XM(2,u2,sp1);			\
    MULT_AND_ACCUM2_XM(2,u2,sp2);			\
  }


{
  constexpr int NCOL=3;
  constexpr int mu=DirX;
  assert(NCOL==m_Nc);

  HOP_XM_A(0);
  HOP_XM_A(1);
  HOP_XM_A(2);

}

#undef HOP_XM_A

#undef SPIN_PROJ1_XM
#undef SPIN_PROJ2_XM
#undef ACCUM1_XM
#undef ACCUM2_XM
