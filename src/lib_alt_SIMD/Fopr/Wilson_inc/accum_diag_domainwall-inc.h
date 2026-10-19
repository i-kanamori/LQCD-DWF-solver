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

#define SET_ACCUM_DIAG(sp) {		 \
    constexpr int offset=(sp-1)*3;	 \
    vcomplex_t hop_c0 = pv[ offset  ];	 \
    vcomplex_t hop_c1 = pv[ offset+1];	 \
    vcomplex_t hop_c2 = pv[ offset+2];	 \
    scal(wc0_sp ## sp, vm);		 \
    axpy(wc0_sp ## sp, vh, hop_c0);	 \
    scal(wc1_sp ## sp, vm);		 \
    axpy(wc1_sp ## sp, vh, hop_c1);	 \
    scal(wc2_sp ## sp, vm);		 \
    axpy(wc2_sp ## sp, vh, hop_c2);	 \
}


{

  constexpr int NS=4; // 4 sprinor
  constexpr int NCOL=3;
  assert(NCOL==m_Nc);
  assert(Ndf_ferm==NS*NCOL);

  vreal_t vh;
  vreal_t vm;
  bcast(vh, -0.5*bb);
  bcast(vm, (4.0-M0)*bb+1.0);

  SET_ACCUM_DIAG(1); // wc*_sp1
  SET_ACCUM_DIAG(2); // wc*_sp2
  SET_ACCUM_DIAG(3); // wc*_sp3
  SET_ACCUM_DIAG(4); // wc*_sp4
}

#undef SET_ACCUM_DIAG

