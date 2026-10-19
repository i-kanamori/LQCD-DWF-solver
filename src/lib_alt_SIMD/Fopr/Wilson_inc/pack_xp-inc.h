/* 
        @file    $Id: pack_xp-inc.h #$

        @brief    Macros for packing in +x direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/

#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_XP(psi) SPIN_PROJ1_XM_(psi) 
#define SPIN_PROJ2_XP(psi) SPIN_PROJ2_XM_(psi) 

#else

#define SPIN_PROJ1_XP(psi) SPIN_PROJ1_XP_(psi) 
#define SPIN_PROJ2_XP(psi) SPIN_PROJ2_XP_(psi) 

#endif

{
  static const int NCOL=3;
  assert(NCOL==m_Nc);
  assert(m_simd_layout.vlen(DirX)==1);

  // spin projection
  const int offset=s*NCOL*4;
  const vcomplex_t *w=pw+offset;
  DEFINE_AND_SET_SPINOR(c,w); // c ## a _sp ## j
  SPIN_PROJ1_XP(c1_sp);
  SPIN_PROJ2_XP(c1_sp);
  SPIN_PROJ1_XP(c2_sp);
  SPIN_PROJ2_XP(c2_sp);
  SPIN_PROJ1_XP(c3_sp);
  SPIN_PROJ2_XP(c3_sp);

  // pack to the buffer
  const int offset_buf=NCOL*s_packed;
  pack_xp(buf_xp,c1_sp1,c1_sp2,offset_buf,0);
  pack_xp(buf_xp,c2_sp1,c2_sp2,offset_buf,1);
  pack_xp(buf_xp,c3_sp1,c3_sp2,offset_buf,2);
  
}

#undef SPIN_PROJ1_XP
#undef SPIN_PROJ2_XP
