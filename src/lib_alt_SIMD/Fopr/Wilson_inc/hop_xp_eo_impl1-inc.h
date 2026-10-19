/* 
        @file    $Id: hop_xp_eo_impl1-inc.h #$

        @brief    Macros for hopping in +x direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/
  
#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_XP(psi) SPIN_PROJ1_XM_(psi) 
#define SPIN_PROJ2_XP(psi) SPIN_PROJ2_XM_(psi) 
#define ACCUM1_XP(b,spb) ACCUM1_XM_(b,spb)
#define ACCUM2_XP(b,spb) ACCUM2_XM_(b,spb)

#else

#define SPIN_PROJ1_XP(psi) SPIN_PROJ1_XP_(psi) 
#define SPIN_PROJ2_XP(psi) SPIN_PROJ2_XP_(psi) 
#define ACCUM1_XP(b,spb) ACCUM1_XP_(b,spb)
#define ACCUM2_XP(b,spb) ACCUM2_XP_(b,spb)

#endif

#define HOP_XP_A( a )		            \
  {				            \
    const int offset0=s*NCOL*4+a;	    \
    const int offset_gauge=s*NCOL*NCOL*4+a; \
    const vcomplex_t *u=Uptr+offset_gauge;  \
    vcomplex_t sp1=pw[offset0];		    \
    vcomplex_t sp2=pw[offset0+NCOL];	    \
    vcomplex_t sp3=pw[offset0+2*NCOL];	    \
    vcomplex_t sp4=pw[offset0+3*NCOL];	    \
    vcomplex_t u0=u[0];			    \
    vcomplex_t u1=u[NCOL];		    \
    vcomplex_t u2=u[2*NCOL];		    \
    SPIN_PROJ1_XP(sp);			    \
    SPIN_PROJ2_XP(sp);			    \
    MULT_AND_ACCUM1_XP(0,u0,sp1);	    \
    MULT_AND_ACCUM2_XP(0,u0,sp2);	    \
    MULT_AND_ACCUM1_XP(1,u1,sp1);	    \
    MULT_AND_ACCUM2_XP(1,u1,sp2);	    \
    MULT_AND_ACCUM1_XP(2,u2,sp1);	    \
    MULT_AND_ACCUM2_XP(2,u2,sp2);	    \
  }



{
  constexpr int NCOL=3;
  constexpr int mu=DirX;
  assert(NCOL==m_Nc);

  HOP_XP_A(0);
  HOP_XP_A(1);
  HOP_XP_A(2);

}

#undef HOP_XP_A

#undef SPIN_PROJ1_XP
#undef SPIN_PROJ2_XP
#undef ACCUM1_XP
#undef ACCUM2_XP
