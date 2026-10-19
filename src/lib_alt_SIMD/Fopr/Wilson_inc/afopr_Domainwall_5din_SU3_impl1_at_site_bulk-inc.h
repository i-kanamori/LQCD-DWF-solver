/*!

        @file    afopr_Domainwall_General_5din_SU3_impl1_at_site_bulk-inc.h

        @brief   Clover Dirac operator with afield (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $

 */
//====================================================================
{
int sss; // working index for macro


/********************************************************************
 //
 // dagger
 //
 *********************************************************************/
#ifdef USE_DAGGER
#define SPINPROJ_USE_DAGGER


/********************************************************************/
// x direction
#ifdef RUN_HOP_XP
#ifdef USE_SAP
{
  const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;     // at s
  GET_GAUGE(u, Uptr, s, DirX );
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xp_impl1-inc.h"
    pw0 += Ndf_ferm;
    pv  += Ndf_ferm;
  }
}
#else // USE_SAP
if(sn_xp>-1){
  const vcomplex_t *pw  = in_base+sn_xp*Ndf_ferm5;  // at sn_xp
  const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;     // at s
  GET_GAUGE(u, Uptr, s, DirX );
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#define HOP_DOES_NOT_ACCUM
#include "multi_hop_xp_impl1-inc.h"
#undef HOP_DOES_NOT_ACCUM
    pw0 += Ndf_ferm;
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
 }
#endif  // else for USE_SPA
#endif  // RUN_HOP_XP

#ifdef RUN_HOP_XM
#ifdef USE_SAP
{
  const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;     // at s
  GET_GAUGE(u,  Uptr, s, DirX );
  SHIFT_AND_ZERO_GAUGE(xm, u);
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xm_impl1-inc.h"
    pw0 += Ndf_ferm;
    pv  += Ndf_ferm;
  }
}
#else // else for USE_SAP
if(sn_xm>-1){
  const vcomplex_t *pw0 = in_base+s*Ndf_ferm5;      // at s
  const vcomplex_t *pw  = in_base+sn_xm*Ndf_ferm5;  // at sn_xm
  GET_GAUGE(u,  Uptr, s, DirX );
  GET_GAUGE(un, Uptr, sn_xm, DirX );
  SHIFT_AND_MERGE_GAUGE(xm, u, un); // (u,un) --> u
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xm_impl1-inc.h"
    pw0 += Ndf_ferm;
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
 }
#endif // USE_SAP
#endif // RUN_HOP_XM


/********************************************************************/
// y direction
if(sn_yp>-1){
#ifdef RUN_HOP_YP
  GET_GAUGE(u, Uptr, s, DirY );
  const vcomplex_t *pw  = in_base+sn_yp*Ndf_ferm5;  // at sn_yp
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_yp-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

if(sn_ym>-1){
#ifdef RUN_HOP_YM
  GET_GAUGE(u,  Uptr, sn_ym, DirY );
  const vcomplex_t *pw  = in_base+sn_ym*Ndf_ferm5;  // at sn_ym
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_ym-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

/********************************************************************/
// z direction
 if(sn_zp>-1){
#ifdef RUN_HOP_ZP
  GET_GAUGE(u, Uptr, s, DirZ );
  const vcomplex_t *pw  = in_base+sn_zp*Ndf_ferm5;  // at sn_zp
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_zp-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

if(sn_zm>-1){
#ifdef RUN_HOP_ZM
  GET_GAUGE(u,  Uptr, sn_zm, DirZ );
  const vcomplex_t *pw  = in_base+sn_zm*Ndf_ferm5;  // at sn_zm
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_zm-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }


/********************************************************************/
// t direction
if(sn_tp>-1){
#ifdef RUN_HOP_TP
  GET_GAUGE(u, Uptr, s, DirT );
  const vcomplex_t *pw  = in_base+sn_tp*Ndf_ferm5;  // at sn_tp
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_tp-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

if(sn_tm>-1){
#ifdef RUN_HOP_TM
  GET_GAUGE(u,  Uptr, sn_tm, DirT );
  const vcomplex_t *pw  = in_base+sn_tm*Ndf_ferm5;  // at sn_tm
  vcomplex_t *pv = pv_site;
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_tm-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
}


/********************************************************************/
// daig. and accumulate
{ // hopping
  prewilson_dag_at_site(pw0_site, pv_site, Ns, cm_array, cp_array, mq);
  const vcomplex_t *pv = pw0_site;

  // diag in Wilson + non_Wilson
  prewilson_dag_at_site(pv_site, in_base+s*Ndf_ferm5, Ns, cm_diag_array, cp_diag_array, mq);

  // output
  for(int s5=0; s5<Ns; s5++){
    real_t bb=b_array[s5];
    GET_SPINOR( pv_site, s5 );
#include "accum_diag_domainwall-inc.h"
    SET_SPINOR( out, s5);
    pv+=Ndf_ferm;
  }
} // hopping
#undef SPINPROJ_USE_DAGGER


/********************************************************************
 //
 // D
 //
 *********************************************************************/
#else  // else for USE_DAGGER
{ // source at s
  const vcomplex_t *in=in_base+s*Ndf_ferm5;
  prewilson_at_site(pw0_site, in, Ns, c_array, mq);
}


/********************************************************************/
// x direction
#ifdef RUN_HOP_XP
#ifdef USE_SAP
{
  const vcomplex_t *pw0 = pw0_site;  // at s
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u, Uptr, s, DirX );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xp_impl1-inc.h"
    pw0 += Ndf_ferm;
    pv  += Ndf_ferm;
  }
}
#else // else for USE_SAP
if(sn_xp>-1){
  const vcomplex_t *in=in_base+sn_xp*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_xp
  const vcomplex_t *pw0 = pw0_site;  // at s
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u, Uptr, s, DirX );
  for(int s5=0; s5<Ns; s5++){
#define HOP_DOES_NOT_ACCUM
#include "multi_hop_xp_impl1-inc.h"
#undef HOP_DOES_NOT_ACCUM
    pw0 += Ndf_ferm;
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
 }
#endif // USE_SAP
#endif // RUN_HOP_XP

#ifdef RUN_HOP_XM
#ifdef USE_SAP
{
  const vcomplex_t *pw0 = pw0_site;  // at s
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u,  Uptr, s, DirX );
  SHIFT_AND_ZERO_GAUGE(xm, u);
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xm_impl1-inc.h"
    pw0 += Ndf_ferm;
    pv  += Ndf_ferm;
  }
}
#else // else for USE_SAP
if(sn_xm>-1){
  const vcomplex_t *in=in_base+sn_xm*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_xm
  const vcomplex_t *pw0 = pw0_site;  // at s
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u,  Uptr, s, DirX );
  GET_GAUGE(un, Uptr, sn_xm, DirX );
  SHIFT_AND_MERGE_GAUGE(xm, u, un); // (u,un) --> u
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_xm_impl1-inc.h"
    pw0 += Ndf_ferm;
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
 }
#endif // USE_SAP
#endif // RUN_HOP_XM


/********************************************************************/
// y direction
if(sn_yp>-1){
#ifdef RUN_HOP_YP
  const vcomplex_t *in=in_base+sn_yp*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_yp
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u, Uptr, s, DirY );
    for(int s5=0; s5<Ns; s5++){
#include "multi_hop_yp-inc.h"
      pw  += Ndf_ferm;
      pv  += Ndf_ferm;
    }
#endif
 }

if(sn_ym>-1){
#ifdef RUN_HOP_YM
  const vcomplex_t *in=in_base+sn_ym*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_ym
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u,  Uptr, sn_ym, DirY );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_ym-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }


/********************************************************************/
// z direction
if(sn_zp>-1){
#ifdef RUN_HOP_ZP
  const vcomplex_t *in=in_base+sn_zp*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_zp
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u, Uptr, s, DirZ );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_zp-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

if(sn_zm>-1){
#ifdef RUN_HOP_ZM
  const vcomplex_t *in=in_base+sn_zm*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_zm
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u,  Uptr, sn_zm, DirZ );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_zm-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }


/********************************************************************/
// t direction
if(sn_tp>-1){
#ifdef RUN_HOP_TP
  const vcomplex_t *in=in_base+sn_tp*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_tp
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u, Uptr, s, DirT );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_tp-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

if(sn_tm>-1){
#ifdef RUN_HOP_TM
  const vcomplex_t *in=in_base+sn_tm*Ndf_ferm5;
  prewilson_at_site(pw_site, in, Ns, c_array, mq);
  const vcomplex_t *pw = pw_site;    // at sn_tm
  vcomplex_t *pv = pv_site;
  GET_GAUGE(u,  Uptr, sn_tm, DirT );
  for(int s5=0; s5<Ns; s5++){
#include "multi_hop_tm-inc.h"
    pw  += Ndf_ferm;
    pv  += Ndf_ferm;
  }
#endif
 }

/********************************************************************/
// diag and accumulate
{
  prewilson_at_site(pw0_site, in_base+s*Ndf_ferm5, Ns, c_diag_array, mq); // non-Wilson + diag in Wilson
  const vcomplex_t *pv = pv_site;   // hopping
  for(int s5=0; s5<Ns; s5++){
    real_t bb=b_array[s5];
    GET_SPINOR( pw0_site, s5 );
#include "accum_diag_domainwall-inc.h"
    SET_SPINOR( out, s5);
    pv+=Ndf_ferm;
  }
}

#endif  // USE_DAGGER
}
