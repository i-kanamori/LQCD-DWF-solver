/*!

        @file    afopr_Domainwall_5din_eo-tmpl.h

        @brief   template file for Domainwall_5din_eo Dirac operator with afield (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */

//====================================================================
#ifndef AFOPR_DOMAINWALL_5DIN_EO_SIMD_TMPL_H_
#define AFOPR_DOMAINWALL_5DIN_EO_SIMD_TMPL_H_


namespace {
#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_EXCLUSIVE_THREAD_FOR_COMM
  inline void set_threadtask_without_master(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    set_threadtask_without_master_afield(ith,nth,is,ns,size);
  }
#else
  inline void set_threadtask_without_master(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    set_threadtask(ith,nth,is,ns,size);
  }
#endif

} // anonymous


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init(const Parameters &param)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);
  init();
  set_parameters(param);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  vout.detailed(m_vl, "afopr_domainwall_general: init, started\n");

  m_repr = "Chiral";  // only the Chiral repr is available.

  m_Nc = CommonParameters::Nc();
  m_Nd = CommonParameters::Nd();

  assert(m_Ndim == CommonParameters::Ndim());
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();
  m_Nx2  = m_Nx/2;
  m_Nst2 = m_Nst/2;

  m_NinF=0; // to be set later after fixing Ns
  m_Ndf_ferm=4*m_Nc;
  m_Ndf_ferm2=m_Ndf_ferm/2;
  m_Ndf_gauge=m_Nd*m_Nc*m_Nc;

  vout.detailed(m_vl, "AFopr_Domainwall_5din_eo: resetting m_U\n");
  m_U.reset(2*m_Nc*m_Nc*m_Ndim, m_Nst, 1);
  vout.detailed(m_vl, "AFopr_Domainwall_5din_eo: resetting m_U, done\n");
  m_Uptr= &(m_U.ptr_vcomplex()[0]);
  m_index.init(m_Nx, m_Ny, m_Nz, m_Nt);

  // even-odd gauge field
  m_Ue.reset(2*m_Nc*m_Nc*m_Ndim, m_Nst2, 1);
  m_Uo.reset(2*m_Nc*m_Nc*m_Ndim, m_Nst2, 1);

  // lattice size in simd vectors
  m_Rx=m_Nx/m_simd_layout.vlen(DirX);
  m_Ry=m_Ny/m_simd_layout.vlen(DirY);
  m_Rz=m_Nz/m_simd_layout.vlen(DirZ);
  m_Rt=m_Nt/m_simd_layout.vlen(DirT);
  m_Rst=m_Rx*m_Ry*m_Rz*m_Rt;
  m_Rx2 = m_Rx/2;
  m_Rst2 = m_Rst/2;

  m_eo_offset=0;

  m_bools.has_mpi_x=false;
  m_bools.has_mpi_y=false;
  m_bools.has_mpi_z=false;
  m_bools.has_mpi_t=false;

#ifdef USE_MPI
  {
    size_t tmp
      =(size_t)m_Nx*Communicator::ipe(DirX)
      +(size_t)m_Ny*Communicator::ipe(DirY)
      +(size_t)m_Nz*Communicator::ipe(DirZ)
      +(size_t)m_Nt*Communicator::ipe(DirT);
    m_eo_offset = tmp %2;
  }

  if(CommonParameters::NPEx()>1){
    m_bools.has_mpi_x=true;
  }
  if(CommonParameters::NPEy()>1){
    m_bools.has_mpi_y=true;
  }
  if(CommonParameters::NPEz()>1){
    m_bools.has_mpi_z=true;
  }
  if(CommonParameters::NPEt()>1){
    m_bools.has_mpi_t=true;
  }

#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_FORCE_COMMUNICATE
  m_bools.has_mpi_x=true;
  m_bools.has_mpi_y=true;
  m_bools.has_mpi_z=true;
  m_bools.has_mpi_t=true;
  vout.general(m_vl, "WARNING!!!!! force to mpi communicate (AFOPR_DOMAINWALL_GENERAL_SIMD_FORCE_COMMUNICATE is defined)\n");
#endif

#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_FORCE_COMMUNICATE_YZT
  m_bools.has_mpi_y=true;
  m_bools.has_mpi_z=true;
  m_bools.has_mpi_t=true;
  vout.general(m_vl, "WARNING!!!!! force to mpi communicate (AFOPR_DOMAINWALL_GENERAL_SIMD_FORCE_COMMUNICATE_YZT is defined)\n");
#endif

#endif // USE_MPI

  print_mpi_bools(m_bools,m_vl);
  print_simd(m_vl);

  // condition check
  if(m_Nx % (2*m_simd_layout.vlen(DirX)) != 0){
    vout.crucial(m_vl, "%s: Nx must be mulriple of 2*simd_layout.vlen(DirX).  Nx=%d, sisimd_layout.vlen(DirX)=%d\n",
                 class_name.c_str(), m_Nx, m_simd_layout.vlen(DirX));
    exit(EXIT_FAILURE);
  }


  init_site_list();

#if defined(AFOPR_DOMAINWALL_5DIN_EO_SIMD_COMMOVERLAP) && (AFOPR_DOMAINWALL_GENERAL_SIMD_COMMOVERLAP == 2)
  vout.general(m_vl, "AFOPR_DOMAINWALL_5DIN_EO_SIMD_COMMOVERLAP=%d (ALL)\n", AFOPR_DOMAINWALL_5DIN_EO_SIMD_COMMOVERLAP);
#endif // (AFOPR_DOMAINWALL_5DIN_EO_COMMOVERLAP) && (AFOPR_DOMAINWALL_5DIN_EO_COMMOVERLAP == 2)

#if defined(AFOPR_DOMAINWALL_5DIN_EO_SIMD_COMMOVERLAP) && (AFOPR_DOMAINWALL_GENERAL_SIMD_COMMOVERLAP == 1)
#error "AFOPR_DOMAINWALL_5DIN_EO_SIMD_COMMOVERLAP=1" is not supported
#endif

#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_USE_TILING
#warning AFOPR_DOMAINWALL_5DIN_EO_SIMD_USE_TILING is set but not supported
#endif


#ifdef TimerAfoprDomainwall_5din_eo_Comm
  timer_comm.reset(new SimpleTimer("  communicatoin (btw send-wait)           "));
#endif
#ifdef TimerAfoprDomainwall_5din_eo
  timer_calc0.reset(new SimpleTimer("interval 0 (packing)                      "));
  timer_calc1.reset(new SimpleTimer("interval 1 (bulk, [1st thread, 1st node]) "));
  timer_calc2.reset(new SimpleTimer("interval 2 (boundary [1st node])          "));
#endif

  etime_comm=0.0;
  etime0=0.0;
  etime1=0.0;
  etime2=0.0;

#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_PREFETCH
  vout.general(m_vl, "afopr_Domainwall_5din_eo: prefetch is ON (AFOPR_DOMAINWALL_5DIN_EO_SIMD_PREFETCH is defined)\n");
#else
  vout.general(m_vl, "afopr_Domainwall_5din_eo: prefetch is OFF (AFOPR_DOMAINWALL_5DIN_EO_SIMD_PREFETCH is not defined)\n");
#endif

#ifdef AFOPR_DOMAINWALL_5DIN_EO_SIMD_EXCLUSIVE_THREAD_FOR_COMM
  vout.general(m_vl, "afopr_Domainwall_5din_eo: the master thread does not compute during communication (AFOPR_DOMAINWALL_5DIN_EO_SIMD_EXCLUSIVE_THREAD_FOR_COMM is defined)\n");
#else
  vout.general(m_vl, "afopr_Domainwall_5din_eo: the master thread computes during communication (AFOPR_DOMAINWALL_5DIN_EO_SIMD_EXCLUSIVE_THREAD_FOR_COMM is not defined)\n");
#endif

#ifdef AFOPR_FIRST_TOUCH_BUF
  vout.general(m_vl, "afopr_Domainwall_5din_eo: 1st touch for the communication buffer is on: (AFOPR_FIRST_TOUCH_BUF is defined)\n");
#else
  vout.general(m_vl, "afopr_Domainwall_5din_eo: 1st touch for the communication buffer is off: (AFOPR_FIRST_TOUCH_BUF is not defined)\n");
#endif


  vout.detailed(m_vl, "allocating buffer for communication, done\n");
}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::tidyup()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);
#ifdef TimerAfoprDomainwall_5din_eo_Comm
  vout.general(m_vl, "%s\n", class_name.c_str());
  timer_comm->report(m_vl);
#endif
#ifdef TimerAfoprDomainwall_5din_eo
  vout.general(m_vl, "%s\n", class_name.c_str());
  timer_calc0->report(m_vl);
  timer_calc1->report(m_vl);
  timer_calc2->report(m_vl);
  timer_calc0->print_resolution(m_vl);
#endif

}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init_site_list()
{
  // assumes:
  //   bulk:     some directions may have commnucation but are swichted off
  //   bundary:  only communcation requireing directoins are withced on

  // prefetch distance for boundary
  constexpr int dist0=1;
  constexpr int dist1=4;

  // initilize the utility
  assert(m_Rx%2==0);
  std::vector<int> Rlattice={m_Rx, m_Ry, m_Rz, m_Rt};
  site_list_utils list(Rlattice, m_bools);
  list.set_prefetch_dist(dist0, dist1);

  // set site index for the bulk
  list.set_site_eo_bulk(m_site_list_bulk, m_eo_offset);
  m_site_list_bulk_ptr[0]=&(m_site_list_bulk[0][0]);
  m_site_list_bulk_ptr[1]=&(m_site_list_bulk[1][0]);

  // set site index for the boundary
  m_Rst_boundary2=0;
#ifdef USE_MPI
  list.set_site_eo_boundary(m_site_list_boundary,
                            m_site_list_pack_x,
                            m_site_list_pack_y,
                            m_site_list_pack_z,
                            m_site_list_pack_t,
                            m_site_list_bulk,
                            m_eo_offset);
  m_Rst_boundary2=m_site_list_boundary[0].size();

  // alias
  if(m_Rst_boundary2>0){
    m_site_list_boundary_ptr[0]=&(m_site_list_boundary[0][0]);
    m_site_list_boundary_ptr[1]=&(m_site_list_boundary[1][0]);
  }

  if(m_site_list_pack_x[0].size()>0){
    m_site_list_pack_x_ptr[0] = &(m_site_list_pack_x[0][0]);
    m_site_list_pack_x_ptr[1] = &(m_site_list_pack_x[1][0]);
  } else {
    m_site_list_pack_x_ptr[0] = 0;
    m_site_list_pack_x_ptr[1] = 0;
  }

  if(m_site_list_pack_y[0].size()>0){
    m_site_list_pack_y_ptr[0] = &(m_site_list_pack_y[0][0]);
    m_site_list_pack_y_ptr[1] = &(m_site_list_pack_y[1][0]);
  } else {
    m_site_list_pack_y_ptr[0] = 0;
    m_site_list_pack_y_ptr[1] = 0;
  }

  if(m_site_list_pack_z[0].size()>0){
    m_site_list_pack_z_ptr[0] = &(m_site_list_pack_z[0][0]);
    m_site_list_pack_z_ptr[1] = &(m_site_list_pack_z[1][0]);
  } else {
    m_site_list_pack_z_ptr[0] = 0;
    m_site_list_pack_z_ptr[1] = 0;
  }

  if(m_site_list_pack_t[0].size()>0){
    m_site_list_pack_t_ptr[0] = &(m_site_list_pack_t[0][0]);
    m_site_list_pack_t_ptr[1] = &(m_site_list_pack_t[1][0]);
  } else {
    m_site_list_pack_t_ptr[0] = 0;
    m_site_list_pack_t_ptr[1] = 0;
  }
#endif

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init_comm_buffer()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

#ifdef USE_MPI
  vout.general(m_vl, "afopr_Domainwall_5din_eo: MPI is on\n");

  // buffer for send/receive:   buffers are made of simd variables
  size_t Nbx= (size_t)m_Ndf_ferm2*m_Ny*m_Nz*m_Nt*m_Ns/vcomplex_t::Nsimd();
  size_t Nby= (size_t)m_Ndf_ferm2*m_Nx*m_Nz*m_Nt*m_Ns/vcomplex_t::Nsimd();
  size_t Nbz= (size_t)m_Ndf_ferm2*m_Nx*m_Ny*m_Nt*m_Ns/vcomplex_t::Nsimd();
  size_t Nbt= (size_t)m_Ndf_ferm2*m_Nx*m_Ny*m_Nz*m_Ns/vcomplex_t::Nsimd();

  if(m_bools.has_mpi_x){
    vout.detailed(m_vl, "allocating buffer for communication: Nbx=%d\n",Nbx);
    comm_buffer.append(&rbuf_xp, DirXp, &buf_xp, DirXm, Nbx, DirXp);
    comm_buffer.append(&rbuf_xm, DirXm, &buf_xm, DirXp, Nbx, DirXm);
  }
  if(m_bools.has_mpi_y){
    vout.detailed(m_vl, "allocating buffer for communication: Nby=%d\n",Nby);
    comm_buffer.append(&rbuf_yp, DirYp, &buf_yp, DirYm, Nby, DirYp);
    comm_buffer.append(&rbuf_ym, DirYm, &buf_ym, DirYp, Nby, DirYm);
  }
  if(m_bools.has_mpi_z){
    vout.detailed(m_vl, "allocating buffer for communication: Nbz=%d\n",Nbz);
    comm_buffer.append(&rbuf_zp, DirZp, &buf_zp, DirZm, Nbz, DirZp);
    comm_buffer.append(&rbuf_zm, DirZm, &buf_zm, DirZp, Nbz, DirZm);
  }
  if(m_bools.has_mpi_t){
    vout.detailed(m_vl, "allocating buffer for communication: Nbt=%d\n",Nbt);
    comm_buffer.append(&rbuf_tp, DirTp, &buf_tp, DirTm, Nbt, DirTp);
    comm_buffer.append(&rbuf_tm, DirTm, &buf_tm, DirTp, Nbt, DirTm);
    vout.general("buf_tp=%p, rbuf_tp=%p\n", buf_tp, rbuf_tp);
  }
  vout.detailed(m_vl, "allocating buffer for communication, done\n");

#else
  vout.general(m_vl, "afopr_Domwainwall_5din_eo: MPI is off\n");
#endif //USE_MPI
}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string str_gmset_type;
  double mq, M0;
  int Ns;
  std::vector<int> bc;
  double b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", str_gmset_type);
  if(!err_optional ){
    vout.detailed(m_vl, "%s: gamma_matrix_type = %s\n", class_name.c_str(), str_gmset_type.c_str());
    if(str_gmset_type != "Chiral"){
      vout.crucial(m_vl, "Error at %s: only Chiral rep. is ready, sorry (given: %s)\n", class_name.c_str(), str_gmset_type.c_str());
      exit(EXIT_FAILURE);
    }
  }
  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c);

}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(
                                const double mq,
                                const double M0,
                                const int Ns,
                                const vector<int> bc,
                                const double b,
                                const double c)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());

  m_M0 = real_t(M0);
  m_mq = real_t(mq);
  m_Ns = Ns;
  m_kappa = real_t(1.0/(8.0 - 2.0 * M0));

  m_Ndf_ferm5 = m_Ndf_ferm * m_Ns;
  m_NinF = 2*m_Ndf_ferm5;


  assert(bc.size() == m_Ndim);
  if(m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
  }

  if(m_b.size() != m_Ns){
    m_b.resize(m_Ns);
    m_c.resize(m_Ns);
  }
  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = real_t(b);
    m_c[is] = real_t(c);
  }

  //- print input parameters
  vout.general(m_vl, "  mq   = %8.4f\n",m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n",m_M0);
  vout.general(m_vl, "  Ns   = %4d\n",m_Ns);
  for(int mu = 0; mu < m_Ndim; ++mu){
    vout.general(m_vl, "  boundary[%d] = %2d\n",mu,m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
                                                  m_b[0], m_c[0]);


  // working 5d vectors.
  m_w1.reset(m_NinF, m_Nst2, 1);
  m_v1.reset(m_NinF, m_Nst2, 1);
  m_v2.reset(m_NinF, m_Nst2, 1);

  int nthread=ThreadManager_OpenMP::get_num_threads_available();
  m_w_at_site_array.resize(nthread);
  m_w0_at_site_array.resize(nthread);
  m_v_at_site_array.resize(nthread);
  m_v0_at_site_array.resize(nthread);

  m_b_array.resize(nthread);
  m_c_array.resize(nthread);
  m_dpinv_array.resize(nthread);
  m_dm_array.resize(nthread);
  m_e_array.resize(nthread);
  m_f_array.resize(nthread);
  m_c_over_b_array.resize(nthread);
  m_c_over_b_diag_array.resize(nthread);
  m_cp_over_b_array.resize(nthread);
  m_cm_over_b_array.resize(nthread);
  m_cp_over_b_diag_array.resize(nthread);
  m_cm_over_b_diag_array.resize(nthread);
  //  m_coeff_m1_array.resize(nthread);

#pragma omp parallel
  { // preapre working area (w/ first touch)
    int ithread=ThreadManager_OpenMP::get_thread_id();
    // local copy of coefficient b, c
    m_b_array[ithread].resize(m_Ns);
    m_c_array[ithread].resize(m_Ns);
    m_dpinv_array[ithread].resize(m_Ns);
    m_dm_array[ithread].resize(m_Ns);
    m_e_array[ithread].resize(m_Ns);
    m_f_array[ithread].resize(m_Ns);
    m_c_over_b_array[ithread].resize(m_Ns);
    m_c_over_b_diag_array[ithread].resize(m_Ns);
    m_cp_over_b_array[ithread].resize(m_Ns);
    m_cm_over_b_array[ithread].resize(m_Ns);
    m_cp_over_b_diag_array[ithread].resize(m_Ns);
    m_cm_over_b_diag_array[ithread].resize(m_Ns);
    //    m_coeff_m1_array[ithread].resize(m_Ns);

    m_w_at_site_array[ithread].resize(m_Ndf_ferm5);
    m_w0_at_site_array[ithread].resize(m_Ndf_ferm5);
    m_v_at_site_array[ithread].resize(m_Ndf_ferm5);
    m_v0_at_site_array[ithread].resize(m_Ndf_ferm5);

    vcomplex_t *pw =&m_w_at_site_array[ithread][0];
    vcomplex_t *pw0=&m_w0_at_site_array[ithread][0];
    vcomplex_t *pv =&m_v_at_site_array[ithread][0];
    vcomplex_t *pv0=&m_v0_at_site_array[ithread][0];

    vcomplex_t vzero;
    set_zero(vzero);
    for(int i=0; i<m_Ndf_ferm5; i++){
      pw[i]=vzero;
      pw0[i]=vzero;
      pv[i]=vzero;
      pv0[i]=vzero;
    }

    m_dm_array[ithread].resize(m_Ns);
    m_e_array[ithread].resize(m_Ns-1);
    m_f_array[ithread].resize(m_Ns-1);

    for(int s5=0; s5<Ns; s5++){
      m_b_array[ithread][s5] = m_b[s5];
      m_c_array[ithread][s5] = m_c[s5];
    }
    set_precond_parameters(ithread);

    real_t mm=(4.0-m_M0);
    for(int s5=0; s5<Ns; s5++){
      //      m_b_array[ithread][s5] = m_b[s5];
      //      m_c_array[ithread][s5] = m_b[s5];
      //      m_dpinv_array[ithread][s5] = m_c[s5]/m_b[s5];
      //      m_dm_array[ithread][s5] = m_c[ (s5+1) % m_Ns] / m_b[s5];
      //      m_e_array[ithread][s5] = m_c[ (s5-1 +m_Ns)% m_Ns] / m_b[s5];
      //      m_f_array[ithread][s5] = m_c[ (s5-1 +m_Ns)% m_Ns] / m_b[s5];
      m_c_over_b_array[ithread][s5] = m_c[s5]/m_b[s5];
      m_cp_over_b_array[ithread][s5] = m_c[ (s5+1) % m_Ns] / m_b[s5];
      m_cm_over_b_array[ithread][s5] = m_c[ (s5-1 +m_Ns)% m_Ns] / m_b[s5];
      m_c_over_b_diag_array[ithread][s5] = (m_c[s5]*mm-1.0) /(mm*m_b[s5] + 1.0);
      m_cp_over_b_diag_array[ithread][s5] = (m_c[ (s5+1) % m_Ns]*mm-1.0) / (mm*m_b[s5] +1.0);
      m_cm_over_b_diag_array[ithread][s5] = (m_c[ (s5-1+m_Ns) % m_Ns]*mm-1.0) / (mm*m_b[s5] + 1.0);

      //      m_coeff_m1_array[ithread][s5]=-1.0;
    }
  }
  init_comm_buffer();
}


template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_coefficients(const std::vector<double> b,
                                                        const std::vector<double> c){
  for(int ithread=0; ithread<m_b_array.size(); ++ithread){
    for(int s5=0; s5<m_Ns; ++s5){
      m_b_array[ithread][s5] = b[s5];
      m_c_array[ithread][s5] = c[s5];
    }
    set_precond_parameters(ithread);
  }
}


template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_precond_parameters(int ithread){

  real_t *b=&m_b_array[ithread][0];
  real_t *c=&m_c_array[ithread][0];

  real_t *dpinv=&m_dpinv_array[ithread][0];
  real_t *dm=&m_dm_array[ithread][0];
  real_t *e=&m_e_array[ithread][0];
  real_t *f=&m_f_array[ithread][0];

  std::vector<real_t> dp(m_Ns);

  for(int is=0; is < m_Ns; ++is){
    dp[is] = 1.0 + b[is]*(4.0 - m_M0);
    dm[is] = 1.0 - c[is]*(4.0 - m_M0);
  }

  e[0] = m_mq * dm[m_Ns-1]/dp[0];
  f[0] = m_mq * dm[0];
  for(int is=1; is < m_Ns-1; ++is){
    e[is] = e[is-1]*dm[is-1]/dp[is];
    f[is] = f[is-1]*dm[is]/dp[is-1];
  }

  real_t g=e[m_Ns-2]*dm[m_Ns-2];
  if(ithread==0){
    m_g=g;
  }

  for(int is=0; is < m_Ns-1; ++is){
    dpinv[is] = 1.0/dp[is];
  }
  dpinv[m_Ns-1] = 1.0/(dp[m_Ns-1]+g);

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config(Field* u)
{
  int nth = ThreadManager::get_num_threads();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_omp(Field* u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_impl(Field* u)
{
#pragma omp master
  {
    m_conf = u;
  }

   // the compiler confuses with the member function convert so "::" is needed
   // convert(m_index, m_U, *u);
   ::convert(m_index, m_U, *u);
#pragma omp barrier
   set_boundary();
   AIndex_eo<real_t, AFIELD::IMPL> index_eo;
#pragma omp barrier
   index_eo.split(m_Ue, m_Uo, m_U);

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_boundary()
{
  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1) {
      real_t bc=(real_t) m_boundary[mu];
      assert(m_Ndim==4);
      int rsize[4] = {m_Rx, m_Ry, m_Rz, m_Rt};
      ::set_boundary_config(m_U, mu, bc, rsize);
    }
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::convert(AFIELD &v, const Field &w)
{
   // the compiler confuses with the member function convert so "::" is needed
   ::convert(m_index, v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::reverse(Field &v, const AFIELD &w)
{
   ::revert(m_index, v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD &v,
                                                 const AFIELD &w,
                                                 std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if(mode=="Deo"){
    D_eo(v, w, 0);
  }else if(mode=="Doe"){
    D_eo(v, w, 1);
  }else if(mode=="Dee"){
    D_ee(v, w, 0);
  }else if(mode=="Doo"){
    D_ee(v, w, 1);
  }else if(mode=="Dee_inv"){
    D_ee_inv(v, w, 0);
  }else if(mode=="Doo_inv"){
    D_ee_inv(v, w, 1);
  }else{
    std::cout << "mode undeifined in AFopr_Domainwall_eo.\n";
    abort();
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD &v,
                                  const AFIELD &w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.  (mult, mode=%s)\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD &v,
                                                     const AFIELD &w,
                                                     std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if(m_mode=="Deo"){
    Ddag_eo(v, w, 0);
  }else if(m_mode=="Doe"){
    Ddag_eo(v, w, 1);
  }else if(m_mode=="Dee"){
    Ddag_ee(v, w, 0);
  }else if(m_mode=="Doo"){
    Ddag_ee(v, w, 1);
  }else if(m_mode=="Dee_inv"){
    Ddag_ee_inv(v, w, 0);
  }else if(m_mode=="Doo_inv"){
    Ddag_ee_inv(v, w, 1);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.  (mult_dag, mode=%s)\n", class_name.c_str(),mode.c_str());
    exit(EXIT_FAILURE);
  }

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD &v,
                                                     const AFIELD &w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.  (mult_dag, mode=%s)\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D(AFIELD &v,
                               const AFIELD &w)
{
  if(debug) {
    vout.general("afopr_Domainwall_5din_eo: D, started\n");
  }

  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);

#pragma omp barrier

  //  vout.crucial("explicit instance is not yet ready for AFopr_Domainwall_5din_eo\n");
  //  exit(EXIT_FAILURE);


}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag(AFIELD &v,
                                                 const AFIELD &w)
{
#pragma omp barrier

  Udag_inv(m_v1, w);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::DdagD(AFIELD &v, const AFIELD &w)
{
  //  vout.crucial(m_vl, "%s: sorry, DdagD is not ready.\n", class_name.c_str());
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#pragma omp barrier

  Udag_inv(v, m_v1);
  Ldag_inv(m_v2, v);
  Ddag_eo(v, m_v2, 1);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
  Ddag_eo(m_v2, v, 0);

  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::DDdag(AFIELD &v, const AFIELD &w)
{
  vout.crucial(m_vl, "%s: sorry, DDdag is not ready.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}




//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::H(AFIELD &v,
                               const AFIELD &w)
{
  if(debug){
    vout.paranoiac("H for AFopr_Domainwall_5din_eo\n");
  }
  vout.crucial(m_vl, "%s: sorry, H is not ready.\n", class_name.c_str());
  exit(EXIT_FAILURE);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
  if(debug) {
    vout.paranoiac("Hdag for AFopr_Domainwall_5din_eo\n");
  }
  vout.crucial(m_vl, "%s: sorry, Hdag is not ready.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}




//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_eo(AFIELD& v, const AFIELD& w,
                                            const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  mult_D_eo(v, w, ieo);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_eo(AFIELD& v, const AFIELD& w,
                                               const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  mult_Ddag_eo(v, w, ieo);
}


template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee_inv(AFIELD& v, const AFIELD& w,
                                                const int ieo)
{
#pragma omp barrier
  L_inv(m_v1, w);
  U_inv(v, m_v1);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee(AFIELD& v, const AFIELD& w,
                                            const int ieo){
  mult_D_ee(v,w,ieo);
}


template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee_inv(AFIELD& v, const AFIELD& w,
                                               const int ieo)
{
#pragma omp barrier
  Udag_inv(m_v1, w);
  Ldag_inv(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee(AFIELD& v, const AFIELD& w,
                                               const int ieo)
{
  mult_Ddag_ee(v,w,ieo);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::L_inv(AFIELD& v, const AFIELD& w)
{
  mult_L_inv(v,w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ldag_inv(AFIELD& v, const AFIELD& w)
{
  mult_Ldag_inv(v,w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::U_inv(AFIELD& v, const AFIELD& w)
{
  mult_U_inv(v,w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Udag_inv(AFIELD& v, const AFIELD& w)
{
  mult_Udag_inv(v,w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_D_eo(AFIELD& v, const AFIELD& w,
                                                 const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  vout.crucial(m_vl, "%s: mult_D_eo is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_Ddag_eo(AFIELD& v, const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  vout.crucial(m_vl, "%s: mult_Ddag_eo is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_D_ee(AFIELD& v, const AFIELD& w,
                                                 const int ieo)
{
  vout.crucial(m_vl, "%s: mult_D_ee is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_Ddag_ee(AFIELD& v, const AFIELD& w,
                                                    const int ieo)
{
  vout.crucial(m_vl, "%s: mult_Ddag_ee is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_L_inv(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_L_inv is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_Ldag_inv(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_Ldag_inv is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_U_inv(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_U_inv is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_Udag_inv(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_Udag_inv is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din_eo<AFIELD>::flop_count(std::string mode)
{
  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  if (m_repr == "Dirac") {
    vout.crucial(m_vl, "flop_count for Dirac rep. is not yet ready\n");
    return 0;
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      //      + 8 * m_Nc * m_Nc * m_Nd * m_Nd
                                    ); // <- clover term
  } else if (m_repr == "Chiral") {
    //    vout.general(m_vl, "hoge: flop_count\n");
    flop_site = static_cast<double>(
                                    m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))  // Wilson
                                    + m_Nc * m_Nd * 2 * 6   // algebra in 5th direction
                                    );
    flop_site *= m_Ns;
  } else {
    //    vout.crucial(m_vl, "%s: input repr is undefined.\n",
    //                 class_name.c_str());
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DdagD_prec")) flop *= 2.0;

  flop=0;
  return flop;
}


// end of codes

//====================================================================

// explicit instanciation.


//============================================================END=====

#endif
