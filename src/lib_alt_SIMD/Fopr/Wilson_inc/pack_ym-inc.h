/* 
        @file    $Id: pack_ym-inc.h #$

        @brief    Macros for packing in -y direction

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/
#ifdef SPINPROJ_USE_DAGGER

#define SPIN_PROJ1_YM(psi) SPIN_PROJ1_YP_(psi) 
#define SPIN_PROJ2_YM(psi) SPIN_PROJ2_YP_(psi) 

#else

#define SPIN_PROJ1_YM(psi) SPIN_PROJ1_YM_(psi) 
#define SPIN_PROJ2_YM(psi) SPIN_PROJ2_YM_(psi) 

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
  SPIN_PROJ1_YM(c1_sp);
  SPIN_PROJ2_YM(c1_sp);
  SPIN_PROJ1_YM(c2_sp);
  SPIN_PROJ2_YM(c2_sp);
  SPIN_PROJ1_YM(c3_sp);
  SPIN_PROJ2_YM(c3_sp);

  vcomplex_t v1_sp1, v2_sp1, v3_sp1;
  vcomplex_t v1_sp2, v2_sp2, v3_sp2;
  { // multiply gauge field
    const int offset=s*NCOL*NCOL*4 + mu*NCOL*NCOL;
    const vcomplex_t *u=Uptr+offset;
    // v_i_sp1 = conj(u_{ji})*c_j_sp1,
    // v_i_sp2 = conj(u_{ji})*c_j_sp2,
    MULT_DAG_2SPIN_SP(v,u,c);
  }

  // pack to the buffer
  const int offset_buf=NCOL*s_packed;
  pack_ym(buf_ym,v1_sp1,v1_sp2,offset_buf,0);
  pack_ym(buf_ym,v2_sp1,v2_sp2,offset_buf,1);
  pack_ym(buf_ym,v3_sp1,v3_sp2,offset_buf,2);
}

#undef SPIN_PROJ1_YM
#undef SPIN_PROJ2_YM
