/*!
      @file    bridgeQXS_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_WILSON_INCLUDED
#define BRIDGEQXS_WILSON_INCLUDED

#ifdef USE_MPREAL
#include "mpreal.h"
using mpfr::mpreal;
#endif

namespace BridgeQXS {
  // real_t = double

  void mult_wilson_bulk_dirac(double *v2, double *up, double *v1,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double *up, double *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(double *v2, double *up, double *v1,
                           double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(double *v2, double *up, double *v1,
                                 double *xp,
                                 double kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double *up, double *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(double *v2, double *up, double *v1,
                              double *xp,
                              double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(double *v2, double *v1, int *Nsize);

  // real_t = float

  void mult_wilson_bulk_dirac(float *v2, float *up, float *v1,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float *up, float *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(float *v2, float *up, float *v1,
                           float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(float *v2, float *up, float *v1,
                                 float *xp,
                                 float kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float *up, float *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(float *v2, float *up, float *v1,
                              float *xp,
                              float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(float *v2, float *v1, int *Nsize);


#ifdef USE_MPREAL
  // real_t = mpreal

  void mult_wilson_bulk_dirac(mpreal *v2, mpreal *up, mpreal *v1,
                              mpreal kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(mpreal *buf_xp, mpreal *buf_xm,
                           mpreal *buf_yp, mpreal *buf_ym,
                           mpreal *buf_zp, mpreal *buf_zm,
                           mpreal *buf_tp, mpreal *buf_tm,
                           mpreal *up, mpreal *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(mpreal *v2, mpreal *up, mpreal *v1,
                           mpreal *buf_xp, mpreal *buf_xm,
                           mpreal *buf_yp, mpreal *buf_ym,
                           mpreal *buf_zp, mpreal *buf_zm,
                           mpreal *buf_tp, mpreal *buf_tm,
                           mpreal kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(mpreal *v2, mpreal *up, mpreal *v1,
                                 mpreal *xp,
                                 mpreal kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(mpreal *buf_xp, mpreal *buf_xm,
                              mpreal *buf_yp, mpreal *buf_ym,
                              mpreal *buf_zp, mpreal *buf_zm,
                              mpreal *buf_tp, mpreal *buf_tm,
                              mpreal *up, mpreal *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(mpreal *v2, mpreal *up, mpreal *v1,
                              mpreal *xp,
                              mpreal *buf_xp, mpreal *buf_xm,
                              mpreal *buf_yp, mpreal *buf_ym,
                              mpreal *buf_zp, mpreal *buf_zm,
                              mpreal *buf_tp, mpreal *buf_tm,
                              mpreal kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(mpreal *v2, mpreal *v1, int *Nsize);

#endif // mpreal
}
#endif
