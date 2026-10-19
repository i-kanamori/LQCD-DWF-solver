/*!

        @file    vsimd_common_SU3-inc.h

        @brief   arithmetics for simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================

//  vsimd_common_SU3-inc.h

namespace {
  // w = u * v
  inline void mult_Vnn(vcomplexF_t& w, vcomplexF_t *u,  vcomplexF_t *v, int Nc)
  {
    vcomplexF_t u0=u[0];
    vcomplexF_t u1=u[1];
    vcomplexF_t u2=u[2];
    vcomplexF_t v0=v[0];
    vcomplexF_t v1=v[1];
    vcomplexF_t v2=v[2];
    set_cmult(w,u0,v0);
    accum_cmult(w, u1,v1);
    accum_cmult(w, u2,v2);
  }

  inline void mult_Vnn(vcomplexD_t& w, vcomplexD_t *u,  vcomplexD_t *v, int Nc)
  {
    vcomplexD_t u0=u[0];
    vcomplexD_t u1=u[1];
    vcomplexD_t u2=u[2];
    vcomplexD_t v0=v[0];
    vcomplexD_t v1=v[1];
    vcomplexD_t v2=v[2];
    set_cmult(w,u0,v0);
    accum_cmult(w, u1,v1);
    accum_cmult(w, u2,v2);
  }

  // w = conj(u) * v
  inline void mult_Vdn(vcomplexF_t& w, vcomplexF_t *u,  vcomplexF_t *v, int Nc)
  {
    vcomplexF_t u0=u[0];
    vcomplexF_t u1=u[1];
    vcomplexF_t u2=u[2];
    vcomplexF_t v0=v[0];
    vcomplexF_t v1=v[1];
    vcomplexF_t v2=v[2];
    set_cmult_conj(w,u0,v0);
    accum_cmult_conj(w, u1,v1);
    accum_cmult_conj(w, u2,v2);
  }

  inline void mult_Vdn(vcomplexD_t& w, vcomplexD_t *u,  vcomplexD_t *v, int Nc)
  {
    vcomplexD_t u0=u[0];
    vcomplexD_t u1=u[1];
    vcomplexD_t u2=u[2];
    vcomplexD_t v0=v[0];
    vcomplexD_t v1=v[1];
    vcomplexD_t v2=v[2];
    set_cmult_conj(w,u0,v0);
    accum_cmult_conj(w, u1,v1);
    accum_cmult_conj(w, u2,v2);
  }

}
