/* 
        @file    $Id: hop_xp_impl1-inc.h #$

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

#ifdef USE_SAP
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_SAP(sp, a)

#elif defined(USE_BLOCK_HOP)
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_BLOCK_HOP(sp, a)

#else
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_(sp, a)
#endif

#define PREPARE_SPINOR_(sp, a)		    \
    vcomplex_t sp ## 1=pw[offset0];	    \
    vcomplex_t sp ## 2=pw[offset0+NCOL];    \
    vcomplex_t sp ## 3=pw[offset0+2*NCOL];  \
    vcomplex_t sp ## 4=pw[offset0+3*NCOL];  \
    vcomplex_t nsp1=pw[offsetn];            \
    vcomplex_t nsp2=pw[offsetn+NCOL];	    \
    vcomplex_t nsp3=pw[offsetn+2*NCOL];	    \
    vcomplex_t nsp4=pw[offsetn+3*NCOL];	    \
    shift_and_merge_xp(sp ## 1,nsp1);	    \
    shift_and_merge_xp(sp ## 2,nsp2);	    \
    shift_and_merge_xp(sp ## 3,nsp3);	    \
    shift_and_merge_xp(sp ## 4,nsp4);   

#define PREPARE_SPINOR_SAP(sp, a)	    \
    vcomplex_t sp ## 1=pw[offset0];	    \
    vcomplex_t sp ## 2=pw[offset0+NCOL];    \
    vcomplex_t sp ## 3=pw[offset0+2*NCOL];  \
    vcomplex_t sp ## 4=pw[offset0+3*NCOL];  \
    vcomplex_t vzero;			    \
    set_zero(vzero);			    \
    shift_and_merge_xp(sp ## 1,vzero);	    \
    shift_and_merge_xp(sp ## 2,vzero);	    \
    shift_and_merge_xp(sp ## 3,vzero);	    \
    shift_and_merge_xp(sp ## 4,vzero); 

#define PREPARE_SPINOR_BLOCK_HOP(sp, a)	\
  vcomplex_t sp ## 1;			\
  vcomplex_t sp ## 2;			\
  vcomplex_t sp ## 3;			\
  vcomplex_t sp ## 4;			\
  set_zero(sp ## 1);			\
  set_zero(sp ## 2);			\
  set_zero(sp ## 3);			\
  set_zero(sp ## 4);			\
  vcomplex_t nsp1=pw[offsetn];		\
  vcomplex_t nsp2=pw[offsetn+NCOL];	\
  vcomplex_t nsp3=pw[offsetn+2*NCOL];   \
  vcomplex_t nsp4=pw[offsetn+3*NCOL];   \
  shift_and_merge_xp(sp ## 1,nsp1);	\
  shift_and_merge_xp(sp ## 2,nsp2);	\
  shift_and_merge_xp(sp ## 3,nsp3);	\
  shift_and_merge_xp(sp ## 4,nsp4);   


#define HOP_XP_A( a )		            \
  {				            \
    const int offsetn=sn*NCOL*4+a;	    \
    const int offset0=s*NCOL*4+a;	    \
    const int offset_gauge=s*NCOL*NCOL*4+a; \
    const vcomplex_t *u=Uptr+offset_gauge;  \
    PREPARE_SPINOR(sp, a);		    \
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

  const int sn=sn_xp;

  HOP_XP_A(0);
  HOP_XP_A(1);
  HOP_XP_A(2);

}

#undef HOP_XP_A
#undef SPIN_PROJ1_XP
#undef SPIN_PROJ2_XP
#undef ACCUM1_XP
#undef ACCUM2_XP
#undef PREPARE_SPINOR
#undef PREPARE_SPINOR_
#undef PREPARE_SPINOR_SAP
#undef PREPARE_SPINOR_BLOCK_HOP
