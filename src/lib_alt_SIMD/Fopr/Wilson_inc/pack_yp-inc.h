/* 
        @file    $Id: pack_yp-inc.h #$

        @brief    Macros for packing in +y direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/
#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_YP(psi) SPIN_PROJ1_YM_(psi) 
#define SPIN_PROJ2_YP(psi) SPIN_PROJ2_YM_(psi) 

#else

#define SPIN_PROJ1_YP(psi) SPIN_PROJ1_YP_(psi) 
#define SPIN_PROJ2_YP(psi) SPIN_PROJ2_YP_(psi) 

#endif
{
  static const int NCOL=3;
  static const int mu=DirY;
  assert(NCOL==m_Nc);
  assert(m_simd_layout.vlen(mu) <3);

  // spin projection
  const int offset=s*NCOL*4;
  const vcomplex_t *w=pw+offset;
  DEFINE_AND_SET_SPINOR(c,w); // c ## a _sp ## j
  SPIN_PROJ1_YP(c1_sp);
  SPIN_PROJ2_YP(c1_sp);
  SPIN_PROJ1_YP(c2_sp);
  SPIN_PROJ2_YP(c2_sp);
  SPIN_PROJ1_YP(c3_sp);
  SPIN_PROJ2_YP(c3_sp);

  // pack to the buffer
  const int offset_buf=NCOL*s_packed;
  pack_yp(buf_yp,c1_sp1,c1_sp2,offset_buf,0);
  pack_yp(buf_yp,c2_sp1,c2_sp2,offset_buf,1);
  pack_yp(buf_yp,c3_sp1,c3_sp2,offset_buf,2);

}
#undef SPIN_PROJ1_YP
#undef SPIN_PROJ2_YP
