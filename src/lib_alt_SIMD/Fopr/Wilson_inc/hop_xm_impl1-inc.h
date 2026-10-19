/*
        @file    $Id: hop_xm_impl1-inc.h #$

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
  
#ifdef USE_SAP
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_SAP(sp, a)
#define PREPARE_GAUGE(u, a) PREPARE_GAUGE_SAP(u, a)

#elif defined(USE_BLOCK_HOP)
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_BLOCK_HOP(sp, a)
#define PREPARE_GAUGE(u, a) PREPARE_GAUGE_BLOCK_HOP(u, a)

#else
#define PREPARE_SPINOR(sp, a) PREPARE_SPINOR_(sp, a)
#define PREPARE_GAUGE(u, a) PREPARE_GAUGE_(u, a)
#endif

#define PREPARE_SPINOR_(sp, a)              \
  vcomplex_t sp ## 1=pw[offset0];           \
  vcomplex_t sp ## 2=pw[offset0+NCOL];      \
  vcomplex_t sp ## 3=pw[offset0+2*NCOL];    \
  vcomplex_t sp ## 4=pw[offset0+3*NCOL];    \
  vcomplex_t nsp1=pw[offsetn];              \
  vcomplex_t nsp2=pw[offsetn+NCOL];         \
  vcomplex_t nsp3=pw[offsetn+2*NCOL];       \
  vcomplex_t nsp4=pw[offsetn+3*NCOL];       \
  shift_and_merge_xm(sp ## 1,nsp1);         \
  shift_and_merge_xm(sp ## 2,nsp2);         \
  shift_and_merge_xm(sp ## 3,nsp3);         \
  shift_and_merge_xm(sp ## 4,nsp4);

#define PREPARE_SPINOR_SAP(sp, a)           \
  vcomplex_t sp ## 1=pw[offset0];           \
  vcomplex_t sp ## 2=pw[offset0+NCOL];      \
  vcomplex_t sp ## 3=pw[offset0+2*NCOL];    \
  vcomplex_t sp ## 4=pw[offset0+3*NCOL];    \
  vcomplex_t vzero;                         \
  set_zero(vzero);                          \
  shift_and_merge_xm(sp ## 1,vzero);        \
  shift_and_merge_xm(sp ## 2,vzero);        \
  shift_and_merge_xm(sp ## 3,vzero);        \
  shift_and_merge_xm(sp ## 4,vzero);

#define PREPARE_SPINOR_BLOCK_HOP(sp, a) \
  vcomplex_t sp ## 1;                   \
  vcomplex_t sp ## 2;                   \
  vcomplex_t sp ## 3;                   \
  vcomplex_t sp ## 4;                   \
  set_zero(sp ## 1);                    \
  set_zero(sp ## 2);                    \
  set_zero(sp ## 3);                    \
  set_zero(sp ## 4);                    \
  vcomplex_t nsp1=pw[offsetn];          \
  vcomplex_t nsp2=pw[offsetn+NCOL];     \
  vcomplex_t nsp3=pw[offsetn+2*NCOL];   \
  vcomplex_t nsp4=pw[offsetn+3*NCOL];   \
  shift_and_merge_xm(sp ## 1,nsp1);     \
  shift_and_merge_xm(sp ## 2,nsp2);     \
  shift_and_merge_xm(sp ## 3,nsp3);     \
  shift_and_merge_xm(sp ## 4,nsp4);


#define PREPARE_GAUGE_(uu, a)            \
  vcomplex_t uu ## 0=u[0];               \
  vcomplex_t uu ## 1=u[1];               \
  vcomplex_t uu ## 2=u[2];               \
  vcomplex_t uu ## n0=un[0];             \
  vcomplex_t uu ## n1=un[1];             \
  vcomplex_t uu ## n2=un[2];             \
  shift_and_merge_xm(uu ## 0, uu ## n0); \
  shift_and_merge_xm(uu ## 1, uu ## n1); \
  shift_and_merge_xm(uu ## 2, uu ## n2);

#define PREPARE_GAUGE_SAP(uu, a)        \
  vcomplex_t uu ## 0=u[0];              \
  vcomplex_t uu ## 1=u[1];              \
  vcomplex_t uu ## 2=u[2];              \
  vcomplex_t vvzero;                     \
  set_zero(vvzero);                      \
  shift_and_merge_xm(uu ## 0, vvzero);   \
  shift_and_merge_xm(uu ## 1, vvzero);   \
  shift_and_merge_xm(uu ## 2, vvzero);

#define PREPARE_GAUGE_BLOCK_HOP(uu, a)   \
  vcomplex_t uu ## 0;                    \
  vcomplex_t uu ## 1;                    \
  vcomplex_t uu ## 2;                    \
  set_zero(uu ## 0);                     \
  set_zero(uu ## 1);                     \
  set_zero(uu ## 2);                     \
  vcomplex_t uu ## n0=un[0];             \
  vcomplex_t uu ## n1=un[1];             \
  vcomplex_t uu ## n2=un[2];             \
  shift_and_merge_xm(uu ## 0, uu ## n0); \
  shift_and_merge_xm(uu ## 1, uu ## n1); \
  shift_and_merge_xm(uu ## 2, uu ## n2);


#define HOP_XM_A( a )                                   \
  {                                                     \
    const int offsetn=sn*NCOL*4+a;                      \
    const int offsetn_gauge=sn*NCOL*NCOL*4+a*NCOL;      \
    const int offset0=s*NCOL*4+a;                       \
    const int offset0_gauge=s*NCOL*NCOL*4+a*NCOL;       \
    const vcomplex_t *un=Uptr+offsetn_gauge;            \
    const vcomplex_t *u=Uptr+offset0_gauge;             \
    PREPARE_SPINOR(sp,a);                               \
    SPIN_PROJ1_XM(sp);                                  \
    SPIN_PROJ2_XM(sp);                                  \
    PREPARE_GAUGE(u,a);                                 \
    MULT_AND_ACCUM1_XM(0,u0,sp1);                       \
    MULT_AND_ACCUM2_XM(0,u0,sp2);                       \
    MULT_AND_ACCUM1_XM(1,u1,sp1);                       \
    MULT_AND_ACCUM2_XM(1,u1,sp2);                       \
    MULT_AND_ACCUM1_XM(2,u2,sp1);                       \
    MULT_AND_ACCUM2_XM(2,u2,sp2);                       \
  }

{
  constexpr int NCOL=3;
  constexpr int mu=DirX;
  assert(NCOL==m_Nc);

  const int sn=sn_xm;

  HOP_XM_A(0);
  HOP_XM_A(1);
  HOP_XM_A(2);

}

#undef HOP_XM_A
#undef SPIN_PROJ1_XM
#undef SPIN_PROJ2_XM
#undef ACCUM1_XM
#undef ACCUM2_XM
#undef PREPARE_SPINOR
#undef PREPARE_SPINOR_
#undef PREPARE_SPINOR_SAP
#undef PREPARE_SPINOR_BLOCK_HOP
#undef PREPARE_GAUGE
#undef PREPARE_GAUGE_
#undef PREPARE_GAUGE_SAP
#undef PREPARE_GAUGE_BLOCK_HOP

