/* 
        @file    $Id: accum_diag_domainwall-inc.h #$

        @brief    Macros for accumulating diaognal elements

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

// w := b( (-M0+4) w0 - (1/2) hop) + w
//  hop: pv

#define SET_DIAG(sp) {			 \
    constexpr int offset=(sp-1)*3;	 \
    vcomplex_t hop_c0 = pv[ offset  ];	 \
    vcomplex_t hop_c1 = pv[ offset+1];	 \
    vcomplex_t hop_c2 = pv[ offset+2];	 \
    wc0_sp ## sp = ::mult(vh, hop_c0);	 \
    wc1_sp ## sp = ::mult(vh, hop_c1);	 \
    wc2_sp ## sp = ::mult(vh, hop_c2);	 \
}

#define ACCUM_DIAG(sp) {			 \
    constexpr int offset=(sp-1)*3;	 \
    vcomplex_t hop_c0 = pv[ offset  ];	 \
    vcomplex_t hop_c1 = pv[ offset+1];	 \
    vcomplex_t hop_c2 = pv[ offset+2];	 \
    axpy(wc0_sp ## sp, vh, hop_c0);	 \
    axpy(wc1_sp ## sp, vh, hop_c1);	 \
    axpy(wc2_sp ## sp, vh, hop_c2);	 \
}


{

  constexpr int NS=4; // 4 sprinor
  constexpr int NCOL=3;
  assert(NCOL==m_Nc);
  assert(Ndf_ferm==NS*NCOL);

  vreal_t vh;
  bcast(vh, -0.5*bb);

#ifdef USE_SET_DIAG
  SET_DIAG(1); // wc*_sp1
  SET_DIAG(2); // wc*_sp2
  SET_DIAG(3); // wc*_sp3
  SET_DIAG(4); // wc*_sp4
#else
  ACCUM_DIAG(1); // wc*_sp1
  ACCUM_DIAG(2); // wc*_sp2
  ACCUM_DIAG(3); // wc*_sp3
  ACCUM_DIAG(4); // wc*_sp4
#endif
}

#undef SET_DIAG
#undef ACCUM_DIAG

