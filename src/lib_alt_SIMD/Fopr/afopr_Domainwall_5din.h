/*!
 
        @file    afopr_Domainwall_5din.h

        @brief   Domainwall_5din Dirac operator with afield (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */
//====================================================================

//#define DEBUG_DOMAINWALL_SIMD

#ifndef AFOPR_DOMAINWALL_5DIN_SIMD_H
#define AFOPR_DOMAINWALL_5DIN_SIMD_H

#include "lib_alt_SIMD/alt_simd_def.h"

//#define TimerAfoprDomainwall_5din_SimdComm
//#define TimerAfoprDomainwall_5din_Simd

// force to copy boundary
#ifdef AFOPR_DOMAINWALL_5DIN_SIMD_FORCE_COMMUNICATE
#undef AFOPR_DOMAINWALL_5DIN_SIMD_FORCE_COMMUNICATE
#endif
#ifdef AFOPR_DOMAINWALL_5DIN_SIMD_FORCE_COMMUNICATE_YZT
#undef AFOPR_DOMAINWALL_5DIN_SIMD_FORCE_COMMUNICATE_YZT
#endif
#ifdef FORCE_COMMUNICATE
//#define AFOPR_DOMAINWALL_5DIN_FORCE_COMMUNICATE
#define AFOPR_DOMAINWALL_5DIN_SIMD_FORCE_COMMUNICATE_YZT
#endif


// the master thread does not compute during the communication
#ifdef AFOPR_DOMAINWALL_5DIN_SIMD_EXCLUSIVE_THREAD_FOR_COMM
#undef AFOPR_DOMAINWALL_5DIN_SIMD_EXCLUSIVE_THREAD_FOR_COMM
#endif
#ifdef EXCLUSIVE_THREAD_FOR_COMM
#define AFOPR_DOMAINWALL_5DIN_SIMD_EXCLUSIVE_THREAD_FOR_COMM
#endif
//#define AFOPR_DOMAINWALL_5DIN_SIMD_EXCLUSIVE_THREAD_FOR_COMM

//#define DEBUG_DOMAINWALL_5DIN_SIMD

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"
#include "lib_alt_SIMD/Field/aindex_lex.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#ifdef TimerAfoprDomainwall_5din_Simd
//#include "lib/Tools/timer.h"
#include "simple_timer.h"
#endif

#include "lib/Communicator/communicator_impl.h"

#include "mpi_bools.h"
#include "site_list.h"
#include "comm_buffer.h"

class Field;

//====================================================================
template<typename AFIELD>
class AFopr_Domainwall_5din : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

  using vreal_t   =ValType_real<real_t>;
  using vcomplex_t=ValType_complex<real_t>;
  using complex_t    =typename vcomplex_t::scalar_type;

  double etime_comm;
  double etime0, etime1, etime2;


#ifdef TimerAfoprDomainwall_5dinComm
  unique_ptr<SimpleTimer> timer_comm;
#endif
#ifdef TimerAfoprDomainwall_5din
  unique_ptr<SimpleTimer> timer_calc0;
  unique_ptr<SimpleTimer> timer_calc1;
  unique_ptr<SimpleTimer> timer_calc2;
#endif


 protected:
  // parameters common to overlap fermion
  real_t m_mq;            //!< quark mass
  real_t m_M0;            //!< domain-wall height
  real_t m_kappa;         //!< kappa for domain-wall height
  int    m_Ns;            //!< size of fifth-dimension
  std::vector<real_t> m_b;
  std::vector<real_t> m_c;

  aligned_vector<aligned_vector<real_t> > m_b_array;
  aligned_vector<aligned_vector<real_t> > m_c_over_b_array;
  aligned_vector<aligned_vector<real_t> > m_cp_over_b_array;
  aligned_vector<aligned_vector<real_t> > m_cm_over_b_array;
  aligned_vector<aligned_vector<real_t> > m_c_over_b_diag_array;
  aligned_vector<aligned_vector<real_t> > m_cp_over_b_diag_array;
  aligned_vector<aligned_vector<real_t> > m_cm_over_b_diag_array;

  aligned_vector<aligned_vector<real_t> > m_coeff_m1_array;


  int m_Nc, m_Nd;

  int m_Nst, m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Rst, m_Rx, m_Ry, m_Rz, m_Rt;   // reduced, after divided by the simd layout
  int m_Rst_boundary;

  int Ndim;
  int m_Ndf_ferm, m_Ndf_gauge, m_Ndf_ferm2;
  int m_Ndf_ferm5;
  int m_NinF;

  int m_nodeid; // for debug

  static constexpr int m_Ndim=4;

  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;        //!< original gauge config.
  AFIELD m_U;           //!< copied gauge config. with boundary conditions.
  vcomplex_t *m_Uptr;   //!< pointer to data member in m_U
  simd_layout_complex<real_t> m_simd_layout; //!< carries information of the simd usage

  std::string  m_mode;  //!< mult mode


  AFIELD m_w1, m_v1, m_v2;    //!< working spinor field.

  AIndex_lex<real_t,AFIELD::IMPL> m_index;

  aligned_vector<site_list_t> m_site_list_bulk;
  aligned_vector<site_list_t> m_site_list_boundary;
  aligned_vector<site_list_pack_t> m_site_list_pack_x;
  aligned_vector<site_list_pack_t> m_site_list_pack_y;
  aligned_vector<site_list_pack_t> m_site_list_pack_z;
  aligned_vector<site_list_pack_t> m_site_list_pack_t;

  aligned_vector<aligned_vector<vcomplex_t> > m_w_at_site_array;
  aligned_vector<aligned_vector<vcomplex_t> > m_w0_at_site_array;
  aligned_vector<aligned_vector<vcomplex_t> > m_v_at_site_array;
  aligned_vector<aligned_vector<vcomplex_t> > m_v0_at_site_array;

  mpi_bools m_bools;

#ifdef DEBUG_DOMAINWALL_5DIN_SIMD
  static constexpr bool debug=true;
#else
  static constexpr bool debug=false;
#endif

  // buffer: body
#ifdef USE_MPI
  CommBuffer<vcomplex_t> comm_buffer;
#endif

  // buffer: alias
  vcomplex_t *buf_xp=nullptr;
  vcomplex_t *rbuf_xp=nullptr;
  vcomplex_t *buf_xm=nullptr;
  vcomplex_t *rbuf_xm=nullptr;

  vcomplex_t *buf_yp=nullptr;
  vcomplex_t *rbuf_yp=nullptr;
  vcomplex_t *buf_ym=nullptr;
  vcomplex_t *rbuf_ym=nullptr;

  vcomplex_t *buf_zp=nullptr;
  vcomplex_t *rbuf_zp=nullptr;
  vcomplex_t *buf_zm=nullptr;
  vcomplex_t *rbuf_zm=nullptr;

  vcomplex_t *buf_tp=nullptr;
  vcomplex_t *rbuf_tp=nullptr;
  vcomplex_t *buf_tm=nullptr;
  vcomplex_t *rbuf_tm=nullptr;


public:
  //! constructor.
  AFopr_Domainwall_5din(){
    vout.paranoiac("AFopr_Domainwall_5din(): calling init()...\n");
    init();
  }

  AFopr_Domainwall_5din(const Parameters& params){
    vout.paranoiac("AFopr_Domainwall_5din(param): calling init()...\n");
    init();
    vout.paranoiac("AFopr_Domainwall_5din(param): calling set_parameters(params)\n");
    set_parameters(params);
  };
  //! destructor.
  ~AFopr_Domainwall_5din(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const double mq, const double M0,
                      const int Ns, const std::vector<int> bc,
                      const double b, const double c);

  //! set parameters of kernel operaotr.
  void set_kernel_parameters(const Parameters& params);

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<double> b,
                        const std::vector<double> c);

  //! this class needs convert of fermion field.
  bool needs_convert(){ return true; }

  void convert(AFIELD&, const Field&);

  void reverse(Field&,  const AFIELD&);

  //! setting gauge configuration.
  void set_config(Field* u);

  void set_config(unique_ptr<Field_G>& U) {set_config(U.get()); }

  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! returns the reference to gauge configuration (AField).
  inline AFIELD *get_confA(void){return &m_U;};

  //! size of domainwall_general term.
  int get_domainwall_general_size();

  void set_boundary_config(AFIELD& U, const int mu);


  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const { return m_mode; }

  std::string get_repr() const { return m_repr; }

  //! mult with specified mode.
  //! Possible modes are: "D", "Ddag", "DdagD",
  //!   "Prec", "D_prec", "Ddag_prec", "DdagD_prec".
  void mult(AFIELD& v, const AFIELD& w, std::string mode);

  //! mult
  void mult(AFIELD& v, const AFIELD& w);

  //! mult
  void mult_dag(AFIELD& v, const AFIELD& w);

  //! mult_dag with specified mode.
  void mult_dag(AFIELD& v, const AFIELD& w, std::string mode);

  void DdagD(AFIELD&, const AFIELD&);
  void DDdag(AFIELD&, const AFIELD&);
  void     D(AFIELD&, const AFIELD&);
  void  Ddag(AFIELD&, const AFIELD&);
  void     H(AFIELD&, const AFIELD&);
  void  Hdag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);
  void mult_gm5R(AFIELD&, const AFIELD&);
  void mult_R(AFIELD&, const AFIELD&);

  //  preconditioner
  void DdagD_prec(AFIELD&, const AFIELD&);
  void     D_prec(AFIELD&, const AFIELD&);
  void  Ddag_prec(AFIELD&, const AFIELD&);
  void Prec(AFIELD&, const AFIELD&);
  void Precdag(AFIELD&, const AFIELD&);

  void L_inv(AFIELD&, const AFIELD&);
  void U_inv(AFIELD&, const AFIELD&);
  void Ldag_inv(AFIELD&, const AFIELD&);
  void Udag_inv(AFIELD&, const AFIELD&);

  //! returns inner size parameter.
  int field_nin() { return m_NinF; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! this returns the number of floating point number operations.
  double flop_count() { return flop_count(m_mode); }

  //! flop-count for specified mode.
  double flop_count(std::string mode);

  //! timer
  void report_timing(){
#ifdef TimerAfoprDomainwall_5din
    vout.general(m_vl, "reproting timing for afopr_Domainwall_5din\n");
    timer_calc0->report();
    timer_calc1->report();
    timer_calc2->report();
    timer_comm->report();
#else
    vout.general("no timer is active: check if TimerAfoprDomainwall_5din is defined in %s.\n",__FILE__);
#endif
  }
 private:
  //! initial setup.
  void init();

  //! initial setup for communication.
  void init_comm_buffer();

  //! initial setup for site list.
  void init_site_list();

  //! final tidy-up.
  void tidyup();

  //! setting boundary condition.
  void set_boundary();

  //  void mult_gm5(real_t*, real_t*);

  //  void clear(real_t*);
  //  void aypx(real_t, real_t*, real_t*);
  //  void gm5_aypx(real_t, real_t*, real_t*);

  //  void clear(real_t*, const int, const int);
  //  void aypx(real_t, real_t*, real_t*, const int, const int);
  //  void gm5_aypx(real_t, real_t*, real_t*, const int, const int);

  // operation on a site: the pointers should point the value on the site
  //  void aypx_site(real_t, vcomplex_t*, const vcomplex_t*, const int);
  //  void scal_site(real_t, vcomplex_t*, const int);
  //  void clear_site(vcomplex_t*, const int);


  //static const int Nc=3;

  // for block operation
  aligned_vector<size_t> m_index_e;
  aligned_vector<size_t> m_index_o;
  aligned_vector<size_t> m_index_xp;
  aligned_vector<size_t> m_index_xm;
  aligned_vector<size_t> m_index_yp;
  aligned_vector<size_t> m_index_ym;
  aligned_vector<size_t> m_index_zp;
  aligned_vector<size_t> m_index_zm;
  aligned_vector<size_t> m_index_tp;
  aligned_vector<size_t> m_index_tm;

#ifdef USE_FACTORY
private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_5din(params); }

public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Domainwall_5din",
                                              create_object_with_params);
    return init1;
  }
#endif


};


#endif  // AFOPR_DOMAINWALL_5DIN_SIMD_H
//============================================================END=====
