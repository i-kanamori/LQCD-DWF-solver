/*!
      @file    afopr_Domainwall_5din_eo_mpreal.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#ifdef USE_MPREAL

#include "mpreal.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

// vector length
#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

using mpfr::mpreal;
typedef mpreal real_t;

#include "lib_alt_QXS/inline/vsimd_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_mpreal-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_qxs_parts-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"

// template file
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo-tmpl.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  AFopr<AField<mpreal, QXS> > *create_object_with_params1(
    const Parameters& params)
  {
    return new AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >(params);
  }


  bool init1 = AFopr<AField<mpreal, QXS> >::Factory_params::Register(
    "Domainwall_5din_eo", create_object_with_params1);
  // temporal for transition of the name
  init1 &= AFopr<AField<mpreal, QXS> >::Factory_params::Register(
    "Domainwall_General_5din_eo", create_object_with_params1);
}
#endif

template<>
const std::string AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >
::class_name = "AFopr_Domainwall_5din_eo<AField<mpreal,QXS> >";


template<>
void AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >::setup_channels()
{
  // use double for sending data
  typedef double real_t;

  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  m_buf_xp.resize(m_Nbdsize[0]);
  m_buf_xm.resize(m_Nbdsize[0]);
  m_buf_yp.resize(m_Nbdsize[1]);
  m_buf_ym.resize(m_Nbdsize[1]);
  m_buf_zp.resize(m_Nbdsize[2]);
  m_buf_zm.resize(m_Nbdsize[2]);
  m_buf_tp.resize(m_Nbdsize[3]);
  m_buf_tm.resize(m_Nbdsize[3]);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif
    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }
  }
}



//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >::mult_D_eo(AField<mpreal, QXS> & v,
                                                 const AField<mpreal, QXS>& w,
                                                 const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AField<mpreal, QXS> *>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector
  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_clear(vp, m_Ns, m_Nsize);

  BridgeQXS::mult_domainwall_5din_eo_5dir_dirac(yp, wp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, do_comm);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = &m_buf_xp[0];
    real_t *buf1_xm = &m_buf_xm[0];
    real_t *buf1_yp = &m_buf_yp[0];
    real_t *buf1_ym = &m_buf_ym[0];
    real_t *buf1_zp = &m_buf_zp[0];
    real_t *buf1_zm = &m_buf_zm[0];
    real_t *buf1_tp = &m_buf_tp[0];
    real_t *buf1_tm = &m_buf_tm[0];

    double *ch_xp = (double*)chsend_dn[0].ptr();
    double *ch_xm = (double*)chsend_up[0].ptr();
    double *ch_yp = (double*)chsend_dn[1].ptr();
    double *ch_ym = (double*)chsend_up[1].ptr();
    double *ch_zp = (double*)chsend_dn[2].ptr();
    double *ch_zm = (double*)chsend_up[2].ptr();
    double *ch_tp = (double*)chsend_dn[3].ptr();
    double *ch_tm = (double*)chsend_up[3].ptr();

    BridgeQXS::mult_domainwall_5din_eo_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, yp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      &m_Leo[0], m_Nsize, do_comm, ieo);

#pragma omp barrier
    {
      int ith, nth, is, ns;
      if(do_comm[0]) {
        set_threadtask(ith, nth, is, ns, m_buf_xp.size());
        for(int i=is; i<ns; ++i) {
          ch_xp[i] = to_double(buf1_xp[i]);
          ch_xm[i] = to_double(buf1_xm[i]);
        }
      }
      if(do_comm[1]) {
        set_threadtask(ith, nth, is, ns, m_buf_yp.size());
        for(int i=is; i<ns; ++i) {
          ch_yp[i] = to_double(buf1_yp[i]);
          ch_ym[i] = to_double(buf1_ym[i]);
        }
      }
      if(do_comm[2]) {
        set_threadtask(ith, nth, is, ns, m_buf_zp.size());
        for(int i=is; i<ns; ++i) {
          ch_zp[i] = to_double(buf1_zp[i]);
          ch_zm[i] = to_double(buf1_zm[i]);
        }
      }
      if(do_comm[3]) {
        set_threadtask(ith, nth, is, ns, m_buf_tp.size());
        for(int i=is; i<ns; ++i) {
          ch_tp[i] = to_double(buf1_tp[i]);
          ch_tm[i] = to_double(buf1_tm[i]);
        }
      }
    }

    #pragma omp barrier


    if (ith == 0) chset_send.start();
  }
#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_hopb_dirac(vp, up, yp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                &m_Leo[0], m_Nsize, do_comm, ieo);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = &m_buf_xp[0];
    real_t *buf2_xm = &m_buf_xm[0];
    real_t *buf2_yp = &m_buf_yp[0];
    real_t *buf2_ym = &m_buf_ym[0];
    real_t *buf2_zp = &m_buf_zp[0];
    real_t *buf2_zm = &m_buf_zm[0];
    real_t *buf2_tp = &m_buf_tp[0];
    real_t *buf2_tm = &m_buf_tm[0];

    double *ch_xp = (double*)chrecv_up[0].ptr();
    double *ch_xm = (double*)chrecv_dn[0].ptr();
    double *ch_yp = (double*)chrecv_up[1].ptr();
    double *ch_ym = (double*)chrecv_dn[1].ptr();
    double *ch_zp = (double*)chrecv_up[2].ptr();
    double *ch_zm = (double*)chrecv_dn[2].ptr();
    double *ch_tp = (double*)chrecv_up[3].ptr();
    double *ch_tm = (double*)chrecv_dn[3].ptr();

    {
      int ith, nth, is, ns;
      if(do_comm[0]) {
        set_threadtask(ith, nth, is, ns, m_buf_xp.size());
        for(int i=is; i<ns; ++i) {
          buf2_xp[i] = ch_xp[i];
          buf2_xm[i] = ch_xm[i];
        }
      }
      if(do_comm[1]) {
        set_threadtask(ith, nth, is, ns, m_buf_yp.size());
        for(int i=is; i<ns; ++i) {
          buf2_yp[i] = ch_yp[i];
          buf2_ym[i] = ch_ym[i];
        }
      }
      if(do_comm[2]) {
        set_threadtask(ith, nth, is, ns, m_buf_zp.size());
        for(int i=is; i<ns; ++i) {
          buf2_zp[i] = ch_zp[i];
          buf2_zm[i] = ch_zm[i];
        }
      }
      if(do_comm[3]) {
        set_threadtask(ith, nth, is, ns, m_buf_tp.size());
        for(int i=is; i<ns; ++i) {
          buf2_tp[i] = ch_tp[i];
          buf2_tm[i] = ch_tm[i];
          //          printf("i=%d  recv, ch_tp=%14.23e ch_tm=%14.23e\n", i, ch_tp[i], ch_tm[i]);
          //          fflush(stdout);
        }
      }
    }

#pragma omp barrier

    BridgeQXS::mult_domainwall_5din_eo_hop2_dirac(vp, up, yp,
                                                  buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                                  buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                                  m_mq, m_M0, m_Ns, &m_boundary[0],
                                                  &m_Leo[0], m_Nsize, do_comm, ieo);
    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}



//====================================================================
template<>
void AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >::mult_Ddag_eo(AField<mpreal, QXS>& v,
                                                    const AField<mpreal, QXS>& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AField<mpreal, QXS> *>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector
  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

  BridgeQXS::mult_domainwall_5din_eo_clear(yp, m_Ns, m_Nsize);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = &m_buf_xp[0];
    real_t *buf1_xm = &m_buf_xm[0];
    real_t *buf1_yp = &m_buf_yp[0];
    real_t *buf1_ym = &m_buf_ym[0];
    real_t *buf1_zp = &m_buf_zp[0];
    real_t *buf1_zm = &m_buf_zm[0];
    real_t *buf1_tp = &m_buf_tp[0];
    real_t *buf1_tm = &m_buf_tm[0];

    BridgeQXS::mult_domainwall_5din_eo_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, vp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      &m_Leo[0], m_Nsize, do_comm, ieo);


#pragma omp barrier
    {
      double *ch_xp = (double*)chsend_dn[0].ptr();
      double *ch_xm = (double*)chsend_up[0].ptr();
      double *ch_yp = (double*)chsend_dn[1].ptr();
      double *ch_ym = (double*)chsend_up[1].ptr();
      double *ch_zp = (double*)chsend_dn[2].ptr();
      double *ch_zm = (double*)chsend_up[2].ptr();
      double *ch_tp = (double*)chsend_dn[3].ptr();
      double *ch_tm = (double*)chsend_up[3].ptr();

      int ith, nth, is, ns;
      if(do_comm[0]) {
        set_threadtask(ith, nth, is, ns, m_buf_xp.size());
        for(int i=is; i<ns; ++i) {
          ch_xp[i] = to_double(buf1_xp[i]);
          ch_xm[i] = to_double(buf1_xm[i]);
        }
      }
      if(do_comm[1]) {
        set_threadtask(ith, nth, is, ns, m_buf_yp.size());
        for(int i=is; i<ns; ++i) {
          ch_yp[i] = to_double(buf1_yp[i]);
          ch_ym[i] = to_double(buf1_ym[i]);
        }
      }
      if(do_comm[2]) {
        set_threadtask(ith, nth, is, ns, m_buf_zp.size());
        for(int i=is; i<ns; ++i) {
          ch_zp[i] = to_double(buf1_zp[i]);
          ch_zm[i] = to_double(buf1_zm[i]);
        }
      }
      if(do_comm[3]) {
        set_threadtask(ith, nth, is, ns, m_buf_tp.size());
        for(int i=is; i<ns; ++i) {
          ch_tp[i] = to_double(buf1_tp[i]);
          ch_tm[i] = to_double(buf1_tm[i]);
        }
      }
    }
#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_domainwall_5din_eo_hopb_dirac(yp, up, vp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                &m_Leo[0], m_Nsize, do_comm, ieo);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = &m_buf_xp[0];
    real_t *buf2_xm = &m_buf_xm[0];
    real_t *buf2_yp = &m_buf_yp[0];
    real_t *buf2_ym = &m_buf_ym[0];
    real_t *buf2_zp = &m_buf_zp[0];
    real_t *buf2_zm = &m_buf_zm[0];
    real_t *buf2_tp = &m_buf_tp[0];
    real_t *buf2_tm = &m_buf_tm[0];

    {
      double *ch_xp = (double*)chrecv_up[0].ptr();
      double *ch_xm = (double*)chrecv_dn[0].ptr();
      double *ch_yp = (double*)chrecv_up[1].ptr();
      double *ch_ym = (double*)chrecv_dn[1].ptr();
      double *ch_zp = (double*)chrecv_up[2].ptr();
      double *ch_zm = (double*)chrecv_dn[2].ptr();
      double *ch_tp = (double*)chrecv_up[3].ptr();
      double *ch_tm = (double*)chrecv_dn[3].ptr();

      int ith, nth, is, ns;
      if(do_comm[0]) {
        set_threadtask(ith, nth, is, ns, m_buf_xp.size());
        for(int i=is; i<ns; ++i) {
          buf2_xp[i] = ch_xp[i];
          buf2_xm[i] = ch_xm[i];
        }
      }
      if(do_comm[1]) {
        set_threadtask(ith, nth, is, ns, m_buf_yp.size());
        for(int i=is; i<ns; ++i) {
          buf2_yp[i] = ch_yp[i];
          buf2_ym[i] = ch_ym[i];
        }
      }
      if(do_comm[2]) {
        set_threadtask(ith, nth, is, ns, m_buf_zp.size());
        for(int i=is; i<ns; ++i) {
          buf2_zp[i] = ch_zp[i];
          buf2_zm[i] = ch_zm[i];
        }
      }
      if(do_comm[3]) {
        set_threadtask(ith, nth, is, ns, m_buf_tp.size());
        for(int i=is; i<ns; ++i) {
          buf2_tp[i] = ch_tp[i];
          buf2_tm[i] = ch_tm[i];
        }
      }
    }

#pragma omp barrier

    BridgeQXS::mult_domainwall_5din_eo_hop2_dirac(yp, up, vp,
                                                  buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                                  buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                                  m_mq, m_M0, m_Ns, &m_boundary[0],
                                                  &m_Leo[0], m_Nsize, do_comm, ieo);
#pragma omp barrier


    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_5dirdag_dirac(vp, yp,
                                                   m_mq, m_M0, m_Ns, &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);

#pragma omp barrier
}


// explicit instanciation
template class AFopr_Domainwall_5din_eo<AField<mpreal, QXS> >;


#endif
//============================================================END=====
