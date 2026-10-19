/*!
      @file    bridgeQXS_Domainwall.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_DOMAINWALL_INCLUDED
#define BRIDGEQXS_DOMAINWALL_INCLUDED

#ifdef USE_MPREAL
#include "bridge_mpreal.h"
using mpfr::mpreal;
#endif

namespace BridgeQXS {
  // real_t = double

  void mult_domainwall_5din_5dir_dirac(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(double *vp, double *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_clear(double *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    double *vp, double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    double *vp, double *up, double *wp,
    double *buf2_xp, double *buf2_xm,
    double *buf2_yp, double *buf2_ym,
    double *buf2_zp, double *buf2_zm,
    double *buf2_tp, double *buf2_tm,
    double mq, double M0, int Ns, int *bc,
    int *Nsize, int *do_comm);


  void mult_domainwall_5din_eo_5dir_dirac(double *yp, double *wp,
                                          double mq, double M0, int Ns, int *bc,
                                          double *b, double *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(double *vp, double *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(double *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    double *vp, double *yp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    double *vp, double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    double *vp, double *up, double *wp,
    double *buf2_xp, double *buf2_xm,
    double *buf2_yp, double *buf2_ym,
    double *buf2_zp, double *buf2_zm,
    double *buf2_tp, double *buf2_tm,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    double *vp, double *up, double *wp,
    double *yp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  // real_t = float

  void mult_domainwall_5din_5dir_dirac(
    float *vp, float *yp, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    float *vp, float *yp, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(float *vp, float *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_clear(float *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    float *vp, float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    float *vp, float *up, float *wp,
    float *buf2_xp, float *buf2_xm,
    float *buf2_yp, float *buf2_ym,
    float *buf2_zp, float *buf2_zm,
    float *buf2_tp, float *buf2_tm,
    float mq, float M0, int Ns, int *bc,
    int *Nsize, int *do_comm);


  void mult_domainwall_5din_eo_5dir_dirac(float *yp, float *wp,
                                          float mq, float M0, int Ns, int *bc,
                                          float *b, float *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(float *vp, float *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(float *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    float *vp, float *yp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    float *vp, float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    float *vp, float *up, float *wp,
    float *buf2_xp, float *buf2_xm,
    float *buf2_yp, float *buf2_ym,
    float *buf2_zp, float *buf2_zm,
    float *buf2_tp, float *buf2_tm,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    float *vp, float *up, float *wp,
    float *yp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

    // real_t
  void mult_domainwall_5din_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

#ifdef USE_MPREAL
  // real_t = mpreal

  void mult_domainwall_5din_5dir_dirac(
    mpreal *vp, mpreal *yp, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    mpreal *vp, mpreal *yp, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(mpreal *vp, mpreal *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_clear(mpreal *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    mpreal *vp, mpreal *up, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    mpreal *buf1_xp, mpreal *buf1_xm,
    mpreal *buf1_yp, mpreal *buf1_ym,
    mpreal *buf1_zp, mpreal *buf1_zm,
    mpreal *buf1_tp, mpreal *buf1_tm,
    mpreal *up, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    mpreal *vp, mpreal *up, mpreal *wp,
    mpreal *buf2_xp, mpreal *buf2_xm,
    mpreal *buf2_yp, mpreal *buf2_ym,
    mpreal *buf2_zp, mpreal *buf2_zm,
    mpreal *buf2_tp, mpreal *buf2_tm,
    mpreal mq, mpreal M0, int Ns, int *bc,
    int *Nsize, int *do_comm);


  void mult_domainwall_5din_eo_5dir_dirac(mpreal *yp, mpreal *wp,
                                          mpreal mq, mpreal M0, int Ns, int *bc,
                                          mpreal *b, mpreal *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(mpreal *vp, mpreal *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(mpreal *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    mpreal *vp, mpreal *yp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    mpreal *vp, mpreal *up, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    mpreal *buf1_xp, mpreal *buf1_xm,
    mpreal *buf1_yp, mpreal *buf1_ym,
    mpreal *buf1_zp, mpreal *buf1_zm,
    mpreal *buf1_tp, mpreal *buf1_tm,
    mpreal *up, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    mpreal *buf1_xp, mpreal *buf1_xm,
    mpreal *buf1_yp, mpreal *buf1_ym,
    mpreal *buf1_zp, mpreal *buf1_zm,
    mpreal *buf1_tp, mpreal *buf1_tm,
    mpreal *up, mpreal *wp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    mpreal *vp, mpreal *up, mpreal *wp,
    mpreal *buf2_xp, mpreal *buf2_xm,
    mpreal *buf2_yp, mpreal *buf2_ym,
    mpreal *buf2_zp, mpreal *buf2_zm,
    mpreal *buf2_tp, mpreal *buf2_tm,
    mpreal mq, mpreal M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    mpreal *vp, mpreal *up, mpreal *wp,
    mpreal *yp,
    mpreal mq, mpreal M0, int Ns, int *bc,
    mpreal *b, mpreal *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);


#endif  // mpreal

}
#endif
