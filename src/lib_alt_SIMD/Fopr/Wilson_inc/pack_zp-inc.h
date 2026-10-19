/* 
        @file    $Id: pack_zp-inc.h #$

        @brief    Macros for packing in +z direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/
#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_ZP(psi) SPIN_PROJ1_ZM_(psi) 
#define SPIN_PROJ2_ZP(psi) SPIN_PROJ2_ZM_(psi) 

#else

#define SPIN_PROJ1_ZP(psi) SPIN_PROJ1_ZP_(psi) 
#define SPIN_PROJ2_ZP(psi) SPIN_PROJ2_ZP_(psi) 

#endif
{
  static const int NCOL=3;
  static const int mu=DirZ;
  assert(NCOL==m_Nc);
  assert(m_simd_layout.vlen(mu) <3);

  // spin projection
  const int offset=s*NCOL*4;
  const vcomplex_t *w=pw+offset;
  DEFINE_AND_SET_SPINOR(c,w); // c ## a _sp ## j
  SPIN_PROJ1_ZP(c1_sp);
  SPIN_PROJ2_ZP(c1_sp);
  SPIN_PROJ1_ZP(c2_sp);
  SPIN_PROJ2_ZP(c2_sp);
  SPIN_PROJ1_ZP(c3_sp);
  SPIN_PROJ2_ZP(c3_sp);

  // pack to the buffer
  const int offset_buf=NCOL*s_packed;
  pack_zp(buf_zp,c1_sp1,c1_sp2,offset_buf,0);
  pack_zp(buf_zp,c2_sp1,c2_sp2,offset_buf,1);
  pack_zp(buf_zp,c3_sp1,c3_sp2,offset_buf,2);
}
#undef SPIN_PROJ1_ZP
#undef SPIN_PROJ2_ZP
