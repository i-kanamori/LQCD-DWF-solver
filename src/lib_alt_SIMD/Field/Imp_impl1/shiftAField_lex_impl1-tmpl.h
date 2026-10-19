/*!

        @file    shiftAField_lex_impl1-tmpl.h

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */


namespace {
#ifdef SHIFT_EXCLUSIVE_THREAD_FOR_COMM
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
}

template<typename AFIELD>
const std::string ShiftAField_lex<AFIELD>::class_name =
                                         "ShiftAField_lex<AFIELD> [impl1]";
//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::init(int Nx, int Ny, int Nz, int Nt, int Nin)
{
  int Ndim = CommonParameters::Ndim();
  std::vector<int> bc(Ndim);
  for(int mu = 0; mu < Ndim; ++mu) bc[mu] = 1;

  init(Nx, Ny, Nz, Nt, Nin, bc);

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::init(int Nx, int Ny, int Nz, int Nt, int Nin, const std::vector<int>& bc)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: being constructed.\n", class_name.c_str());

  m_Nin = Nin;
  vout.general(m_vl, "  Nin = %d\n", m_Nin);

  m_Nx = Nx;
  m_Ny = Ny;
  m_Nz = Nz;
  m_Nt = Nt;
  m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;

  m_Ndim = CommonParameters::Ndim();


  m_Nxv  = m_Nx/VLEN2;
  m_Nstv = m_Nvol/VLEN2;
  vout.general(m_vl, "  Nxv = %d  Nstv = %d\n", m_Nxv, m_Nstv);

  if(bc.size() != m_Ndim){
    vout.crucial(m_vl, "%s: incorrect size of boundary condition\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_boundary.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#ifdef SHIFT_FORCE_COMMNUICATE
  int req_comm = 1;  // set 1 if communication forced any time
#else
  int req_comm = 0;
#endif
  do_comm_any = 0;
  do_comm[0] = 0;
  vout.general("  do_comm[%d] = %d\n", 0, do_comm[0]);
  for(int mu = 1; mu < m_Ndim; ++mu){
    do_comm[mu]=1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

#ifdef USE_MPI
  // buffer for send/receive:   buffers are made of simd variables
  std::vector<size_t> buf_size(4);
  buf_size[0] = 0;
  buf_size[1] = do_comm[1]*(size_t)m_Nin*m_Nx*m_Nz*m_Nt/vcomplex_t::Nsimd()/2;
  buf_size[2] = do_comm[2]*(size_t)m_Nin*m_Nx*m_Ny*m_Nt/vcomplex_t::Nsimd()/2;
  buf_size[3] = do_comm[3]*(size_t)m_Nin*m_Nx*m_Ny*m_Nz/vcomplex_t::Nsimd()/2;

  setup_comm(buf_size);
#endif

  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());

}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::setup_comm(const std::vector<size_t> &buf_size)
{
#ifdef USE_MPI

#pragma omp master
  {
    // make buffer
    size_t bsize_total=0;
    for(int mu=0; mu<m_Ndim; ++mu){
      bsize_total+=buf_size[mu];
      buf_tag[2*mu  ]=0;
      buf_tag[2*mu+1]=0;
    }
    buf_up_body.resize(bsize_total);
    buf_dn_body.resize(bsize_total);
    rbuf_up_body.resize(bsize_total);
    rbuf_dn_body.resize(bsize_total);

    // pointer to the buffer
    size_t offset=0;
    for(int mu=0; mu<m_Ndim; ++mu){
      if(buf_size[mu]>0){
        buf_tag[2*mu  ]=&(buf_dn_body[offset]);
        buf_tag[2*mu+1]=&(buf_up_body[offset]);
        rbuf_tag[2*mu  ]=&(rbuf_up_body[offset]);
        rbuf_tag[2*mu+1]=&(rbuf_dn_body[offset]);
        buf_up[mu]=&(buf_up_body[offset]);
        buf_dn[mu]=&(buf_dn_body[offset]);
        rbuf_up[mu]=&(rbuf_up_body[offset]);
        rbuf_dn[mu]=&(rbuf_dn_body[offset]);
      }
      offset+=buf_size[mu];
    }

    // list of neighbor nodes
    int myrank=Communicator::self();

    ipe_neighbor.resize(2*m_Ndim);
    aligned_vector<int> coord(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      Communicator::grid_coord(&coord[0], Communicator::self());

      int coord_keep=coord[mu];
      // upward
      coord[mu] = (coord_keep + 1) % Communicator::npe(mu);
      Communicator::grid_rank(&ipe_neighbor[2*mu]  ,&coord[0]);

      // downward
      coord[mu] = (coord_keep - 1 + Communicator::npe(mu)) % Communicator::npe(mu);
      Communicator::grid_rank(&ipe_neighbor[2*mu+1]  ,&coord[0]);
    }

    m_request_send.resize(2*m_Ndim);
    m_request_recv.resize(2*m_Ndim);
    m_status.resize(2*m_Ndim);

    // start persistent communication
    for(int mu=1; mu< m_Ndim; ++mu){ // do not communicate in x-direction
      size_t size=buf_size[mu]*vcomplex_t::Nsimd()*2;  // buf_size is counted as an array of vcomplex_t
      int proc_p = ipe_neighbor[2*mu];
      int proc_m = ipe_neighbor[2*mu+1];

      { // recieve from up (send to down)
        int tag=2*mu;
        int recv_from=proc_p;
        int send_to=proc_m;
        void *buf=reinterpret_cast<void*>(buf_tag[tag]);
        void *rbuf=reinterpret_cast<void*>(rbuf_tag[tag]);

        int retv;
        retv = MPI_Send_init(buf, sizeof(real_t)*size, MPI_BYTE, send_to, tag, Communicator_impl::world(), &m_request_send[tag]);
        if (retv != MPI_SUCCESS) {
          printf("MPI_Send_init tag=%d, failed with %d in %s\n", tag, retv, __FILE__);
          fflush(0);
          exit(EXIT_FAILURE);
        }

        retv = MPI_Recv_init(rbuf, sizeof(real_t)*size, MPI_BYTE, recv_from, tag, Communicator_impl::world(), &m_request_recv[tag]);
        if (retv != MPI_SUCCESS) {
          printf("MPI_Recv_init tag=%d, failed with %d in %s\n", tag, retv, __FILE__);
          fflush(0);
          exit(EXIT_FAILURE);
        }
      }

      { // recieve from down (send to up)
        int tag=2*mu+1;
        int recv_from=proc_m;
        int send_to=proc_p;
        void *buf=reinterpret_cast<void*>(buf_tag[tag]);
        void *rbuf=reinterpret_cast<void*>(rbuf_tag[tag]);

        int retv;
        retv = MPI_Send_init(buf, sizeof(real_t)*size, MPI_BYTE, send_to, tag, Communicator_impl::world(), &m_request_send[tag]);
        if (retv != MPI_SUCCESS) {
          printf("MPI_Send_init tag=%d, failed with %d in %s\n", tag, retv, __FILE__);
          fflush(0);
          exit(EXIT_FAILURE);
        }

        retv = MPI_Recv_init(rbuf, sizeof(real_t)*size, MPI_BYTE, recv_from, tag, Communicator_impl::world(), &m_request_recv[tag]);
        if (retv != MPI_SUCCESS) {
          printf("MPI_Recv_init tag=%d, failed with %d in %s\n", tag, retv, __FILE__);
          fflush(0);
          exit(EXIT_FAILURE);
        }
      }

    } // mu loop

  } // omp master
#endif
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::setup_channels()
{
  /*
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    int Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu] = Communicator::send_init(Nvsize, mu, -1);
    chsend_up[mu] = Communicator::send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1);
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu]->ptr();
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu]->ptr();
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1, buf_dn);
#endif
    if(do_comm[mu] == 1){
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

  }
  */
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const AFIELD& w,
                                                         const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  Element_type::type elem_type=w.field_element_type();
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  assert(elem_type==Element_type::COMPLEX);
  int ANin=m_Nin;
  if(elem_type==Element_type::COMPLEX){
    ANin/=2;
    for(int ex = 0; ex < Nex; ++ex){
      vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex);
      const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex);

      if(mu == 0) {
        up_x(vp, wp, ANin, ANin, ANin, index);
      }else if(mu == 1) {
        up_y(vp, wp, ANin, ANin, ANin, index);
      }else if(mu == 2) {
        up_z(vp, wp, ANin, ANin, ANin, index);
      }else if(mu == 3) {
      	up_t(vp, wp, ANin, ANin, ANin, index);
      }else{
        vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                     class_name.c_str());
        exit(EXIT_FAILURE);
      }
    }
  } else { // real
    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

    for(int ex = 0; ex < Nex; ++ex){
      vreal_t *vp = v.ptr_vreal()+index.vector_idx(0, ANin, 0, ex);
      const vreal_t *wp = w.ptr_vreal()+index.vector_idx(0, ANin, 0, ex);
      //...
    }
  } // end real
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const int ex1,
                                 const AFIELD& w, const int ex2,
                                                  const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);


  Element_type::type elem_type=w.field_element_type();
  assert(elem_type==Element_type::COMPLEX);
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  int ANin=m_Nin;
  if(elem_type==Element_type::COMPLEX){
    ANin/=2;
    vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex1);
    const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex2);
    if(mu == 0) {
      up_x(vp, wp, ANin, ANin, ANin, index);
    }else if(mu == 1) {
      up_y(vp, wp, ANin, ANin, ANin, index);
    }else if(mu == 2) {
      up_z(vp, wp, ANin, ANin, ANin, index);
    }else if(mu == 3) {
      up_t(vp, wp, ANin, ANin, ANin, index);
    }else{
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  } else { // real
    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

    vreal_t *vp = v.ptr_vreal()+index.vector_idx(0, ANin, 0, ex1);
    const vreal_t *wp = w.ptr_vreal()+index.vector_idx(0, ANin, 0, ex2);
    //...
  }


}



//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const int in_start1, const int ex1,
                                       const AFIELD& w, const int in_start2, const int ex2,
                                       const int size, const int mu)
{
  // assumes complex
  const int ANin1=v.nin()/2;
  const int ANin2=w.nin()/2;

  assert(in_start1+size<=ANin1);
  assert(ex1 < v.nex());
  assert(in_start2+size<=ANin2);
  assert(ex2 < w.nex());

  Element_type::type elem_type=w.field_element_type();
  assert(elem_type==Element_type::COMPLEX);
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  if(elem_type==Element_type::COMPLEX){
    vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin1, 0, ex1)+in_start1;
    const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin2, 0, ex2)+in_start2;
    if(mu == 0) {
      up_x(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 1) {
      up_y(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 2) {
      up_z(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 3) {
      up_t(vp, wp, ANin1, ANin2, size, index);
    }else{
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  } else { // real
    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

  }


}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const AFIELD& w,
                                                         const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  Element_type::type elem_type=w.field_element_type();
  assert(elem_type==Element_type::COMPLEX);
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  int ANin=m_Nin;
  if(elem_type==Element_type::COMPLEX){
    ANin/=2;
    for(int ex = 0; ex < Nex; ++ex){
      vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex);
      const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex);

      if (mu == 0) {
        dn_x(vp, wp, ANin, ANin, ANin, index);
      } else if (mu == 1) {
        dn_y(vp, wp, ANin, ANin, ANin, index);
      } else if (mu == 2) {
        dn_z(vp, wp, ANin, ANin, ANin, index);
      } else if (mu == 3) {
        dn_t(vp, wp, ANin, ANin, ANin, index);
      } else {
        vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                     class_name.c_str());
        exit(EXIT_FAILURE);
      }
    } // ex-loop

  } else { // real

    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

    for(int ex = 0; ex < Nex; ++ex){
      vreal_t *vp = v.ptr_vreal()+index.vector_idx(0, ANin, 0, ex);
      const vreal_t *wp = w.ptr_vreal()+index.vector_idx(0, ANin, 0, ex);
      //...
    }
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const int ex1,
                                const AFIELD& w, const int ex2,
                                                 const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);


  Element_type::type elem_type=w.field_element_type();
  assert(elem_type==Element_type::COMPLEX);
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  int ANin=m_Nin;
  if(elem_type==Element_type::COMPLEX){
    ANin/=2;
    vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex1);
    const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin, 0, ex2);
    if (mu == 0) {
      dn_x(vp, wp, ANin, ANin, ANin, index);
    } else if (mu == 1) {
      dn_y(vp, wp, ANin, ANin, ANin, index);
    } else if (mu == 2) {
      dn_z(vp, wp, ANin, ANin, ANin, index);
    } else if (mu == 3) {
      dn_t(vp, wp, ANin, ANin, ANin, index);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  } else {
    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const int in_start1, const int ex1,
                                      const AFIELD& w, const int in_start2, const int ex2,
                                      const int size, const int mu)
{
  // assume complex
  const int ANin1=v.nin_complex();
  const int ANin2=w.nin_complex();

  assert(in_start1+size<=ANin1);
  assert(ex1 < v.nex());
  assert(in_start2+size<=ANin2);
  assert(ex2 < w.nex());

  Element_type::type elem_type=w.field_element_type();
  assert(elem_type==Element_type::COMPLEX);
  index_t index(m_Nx, m_Ny, m_Nz, m_Nt, elem_type);
  if(elem_type==Element_type::COMPLEX){
    vcomplex_t *vp = v.ptr_vcomplex()+index.vector_idx(0, ANin1, 0, ex1)+in_start1;
    const vcomplex_t *wp = w.ptr_vcomplex()+index.vector_idx(0, ANin2, 0, ex2)+in_start2;
    if(mu == 0) {
      dn_x(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 1) {
      dn_y(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 2) {
      dn_z(vp, wp, ANin1, ANin2, size, index);
    }else if(mu == 3) {
      dn_t(vp, wp, ANin1, ANin2, size, index);
    }else{
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  } else { // real
    vout.crucial(m_vl, "%s: shift for real field is not ready\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

  }


}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_x(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(CommonParameters::NPEx()==1);
  assert(size<=ANin1);
  assert(size<=ANin2);

  real_t bc2 = real_t(m_boundary[0]);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
  vreal_t boundary1;
  vreal_t boundary2;
  bcast(boundary1, 1.0);
  set_boundary_xp(boundary2, bc2);

#pragma omp barrier
  for(int site=is; site<ns; ++site){
    vreal_t boundary=boundary1;
    int ix=site%Rx;
    if(ix==Rx-1) {
      boundary=boundary2;
    }
    //    int iyzt=site/Rx;
    //    int site_n= (ix+1) % Rx + Rx*iyzt;
    int site_n = (site-ix) + (ix+1) % Rx;
    for(int nin=0; nin<size; ++nin){
      int i1 =site*ANin1+nin;
      int i2 =site  *ANin2+nin;
      int i2n=site_n*ANin2+nin;
      vcomplex_t val1=wp[i2];
      vcomplex_t val2=wp[i2n];
      shift_and_merge_xp(val1,val2);
      apply_boundary(val1, boundary);
      vp[i1]=val1;
    }
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_x(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(CommonParameters::NPEx()==1);
  assert(size<=ANin1);
  assert(size<=ANin2);

  //  real_t bc2 = 1.0;
  //  if (Communicator::ipe(0) == 0) bc2 = real_t(m_boundary[0]);
  real_t bc2 = real_t(m_boundary[0]);
  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
  vreal_t boundary1;
  vreal_t boundary2;
  bcast(boundary1, 1.0);
  set_boundary_xm(boundary2, bc2);

#pragma omp barrier
  for(int site=is; site<ns; ++site){
    vreal_t boundary=boundary1;
    int ix=site%Rx;
    if(ix==0) {
      boundary=boundary2;
    }
    //int iyzt=site/Rx;
    //    int site_n= (ix-1+Rx) % Rx + Rx*iyzt;
    int site_n = (site-ix) + (ix-1+Rx) % Rx;
    for(int nin=0; nin<size; ++nin){
      int i1 = site  *ANin1+nin;
      int i2  = site  *ANin2+nin;
      int i2n = site_n*ANin2+nin;
      vcomplex_t val1=wp[i2];
      vcomplex_t val2=wp[i2n];
      shift_and_merge_xm(val1,val2);
      apply_boundary(val1, boundary);
      vp[i1]=val1;
    }
  }
#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_y(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);

  const int mu=YDIR;
  const int tag=2*mu;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == Communicator::npe(mu)-1) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_yp(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier

#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
#else
  static const bool need_comm=false;
#endif

  if(need_comm) {
    #ifdef USE_MPI
    // prepare buffer
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*1 *Rz*Rt); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole

    // s_packed= rx + Rx*(rz + Rz*rt) , tmp=rz + Rz*rt
    // site    = rx + Rx*(0 + Ry*(rz + Rz*rt) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/Rx;
      int site=s_packed + Rx*(Ry-1)*tmp;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
      }
    }
#pragma omp barrier

#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    // site    = rx + Rx*(iy   + Ry*(rz + Rz*rt) )
    // site_n  = rx + Rx*(iy+1 + Ry*(rz + Rz*rt) )  [ iy\neq Ry-1]
    for(int site=is; site<ns; ++site){
      int iy=(site/Rx) % Ry;
      if(iy==Ry-1) {
        continue;
      }
      int site_n= site+Rx;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    // s_packed= rx + Rx*(rz + Rz*rt) , tmp=rz + Rz*rt
    // site    = rx + Rx*(Ry-1 + Ry*(rz + Rz*rt) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/Rx;
      int site=s_packed + Rx*(Ry-1)*(tmp+1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed
    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int iy=(site/Rx) % Ry;
      int site_n= site+Rx;
      if(iy==Ry-1) {
        boundary=boundary2;
        site_n-=Rx*Ry;
      }
      //      int site_n= site+Rx*(-iy + ((iy+1)%Ry));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_y(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);

  const int mu=YDIR;
  const int tag=2*mu+1;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == 0) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_ym(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier

#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
  //    const bool need_comm=true;
#else
  static const bool need_comm=false;
#endif
  if(need_comm) {
#ifdef USE_MPI
    // prepare buffer
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*1 *Rz*Rt); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole

    // s_packed= rx + Rx*(rz + Rz*rt) , tmp=rz + Rz*rt
    // site    = rx + Rx*(Ry-1 + Ry*(rz + Rz*rt) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/Rx;
      int site=s_packed + Rx*(Ry-1)*(tmp+1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
      }
    }

#pragma omp barrier
#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    // site    = rx + Rx*(iy   + Ry*(rz + Rz*rt) )
    // site_n  = rx + Rx*(iy-1 + Ry*(rz + Rz*rt) )  [ iy\neq 0]
    for(int site=is; site<ns; ++site){
      int iy=(site/Rx) % Ry;
      if(iy==0) {
        continue;
      }
      int site_n= site-Rx;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    // s_packed= rx + Rx*(rz + Rz*rt) , tmp=rz + Rz*rt
    // site    = rx + Rx*(0 + Ry*(rz + Rz*rt) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/Rx;
      int site=s_packed + Rx*(Ry-1)*tmp;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed

    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int iy=(site/Rx) % Ry;
      int site_n= site-Rx;
      if(iy==0) {
        boundary=boundary2;
        site_n+=Rx*Ry;
      }
      //      int site_n= site+Rx*(-iy + ((iy-1+Ry)%Ry));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_z(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);
  const int mu=ZDIR;
  const int tag=2*mu;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == Communicator::npe(mu)-1) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_zp(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier

#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
#else
  static const bool need_comm=false;
#endif

  if(need_comm) {
    #ifdef USE_MPI
    // prepare buffer
    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*Ry*1*Rt); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/(Rx*Ry);  // =t
      int site=s_packed + Rx*Ry*(Rz-1)*tmp;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
      }
    }
#pragma omp barrier
#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    for(int site=is; site<ns; ++site){
      int iz=(site/(Rx*Ry)) % Rz;
      if(iz==Rz-1) {
        continue;
      }
      int site_n= site+Rx*Ry;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/(Rx*Ry);
      int site=s_packed + Rx*Ry*(Rz-1)*(tmp+1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed
    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int iz=(site/(Rx*Ry)) % Rz;
      int site_n= site+Rx*Ry;
      if(iz==Rz-1) {
        boundary=boundary2;
        site_n-=Rx*Ry*Rz;
      }
      //      int site_n= site+Rx*Ry*(-iz + ((iz+1)%Rz));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_z(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);
  const int mu=ZDIR;
  const int tag=2*mu+1;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == 0) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_zm(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier

#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
#else
  static const bool need_comm=false;
#endif


  if(need_comm) {
    #ifdef USE_MPI
    // prepare buffer
    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*Ry*1*Rt); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/(Rx*Ry);  // =t
      int site=s_packed + Rx*Ry*(Rz-1)*(tmp+1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
        //buf[i]=wp[j];
      }
    }
#pragma omp barrier
#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    for(int site=is; site<ns; ++site){
      int iz=(site/(Rx*Ry)) % Rz;
      if(iz==0) {
        continue;
      }
      int site_n= site-Rx*Ry;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int tmp=s_packed/(Rx*Ry);
      int site=s_packed + Rx*Ry*(Rz-1)*tmp;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed
    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int iz=(site/(Rx*Ry)) % Rz;
      int site_n= site-Rx*Ry;
      if(iz==0) {
        boundary=boundary2;
        site_n+=Rx*Ry*Rz;
      }
      //      int site_n= site+Rx*Ry*(-iz + ((iz-1+Rz)%Rz));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_t(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);
  const int mu=TDIR;
  const int tag=2*mu;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == Communicator::npe(mu)-1) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_yp(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier
#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
#else
  static const bool need_comm=false;
#endif

  if(need_comm) {
    #ifdef USE_MPI
    // prepare buffer
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*Ry*Rz*1); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole

    // s_packed= rx + Rx*(ry + Ry*rz)
    // site    = rx + Rx*(ry + Ry*(rz + Rz*0) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int site=s_packed;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
      }
    }
#pragma omp barrier

#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    // site    = rx + Rx*(ry + Ry*(rz + Rz*it) )
    // site_n  = rx + Rx*(ry + Ry*(rz + Rz*(it+1) )  [ it\neq Rt-1]
    for(int site=is; site<ns; ++site){
      int it=site/(Rx*Ry*Rz);
      if(it==Rt-1) {
        continue;
      }
      int site_n= site+Rx*Ry*Rz;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    // s_packed= rx + Rx*(ry + Ry*rz)
    // site    = rx + Rx*(ry + Ry*(rz + Rz*(Rt-1)) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int site=s_packed + Rx*Ry*Rz*(Rt-1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed
    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int it=site/(Rx*Ry*Rz);
      int site_n= site+Rx*Ry*Rz;
      if(it==Rt-1) {
        boundary=boundary2;
        site_n -= Rx*Ry*Rz*Rt;
      }
      //      int site_n= site+Rx*Ry*Rz * (-it + ((it+1)%Rt));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier

}

//====================================================================
//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_t(vcomplex_t *vp, const vcomplex_t *wp,
                                   const int ANin1, const int ANin2,
                                   const int size, const index_t &index)
{
  assert(size<=ANin1);
  assert(size<=ANin2);
  assert(2*size<=m_Nin);
  const int mu=TDIR;
  const int tag=2*mu+1;

  real_t bc2 = 1.0;
  if (Communicator::ipe(mu) == 0) { bc2 = real_t(m_boundary[mu]); }
  vreal_t boundary2;
  set_boundary_ym(boundary2, bc2);

  const int Rx=index.rlattice(XDIR);
  const int Ry=index.rlattice(YDIR);
  const int Rz=index.rlattice(ZDIR);
  const int Rt=index.rlattice(TDIR);

#pragma omp barrier

#ifdef USE_MPI
  const bool need_comm=( do_comm[mu] == 1 );
  //    const bool need_comm=true;
#else
  static const bool need_comm=false;
#endif
  if(need_comm) {
#ifdef USE_MPI
    // prepare buffer
    vcomplex_t *buf=buf_tag[tag];
    vcomplex_t *rbuf=rbuf_tag[tag];

    int ith, nth, is, ns,is1,ns1;
    set_threadtask(ith, nth, is1, ns1, Rx*Ry*Rz*1); // for boundary
    set_threadtask_without_master(ith, nth, is,  ns,  Rx*Ry*Rz*Rt); // for whole

    // s_packed= rx + Rx*(ry + Ry*rz)
    // site    = rx + Rx*(ry + Ry*(rz + Rz*(Rt-1)) )

    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int site=s_packed + Rx*Ry*Rz*(Rt-1);
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin2+nin;
        sstore(buf+i, wp[j]);
      }
    }

#pragma omp barrier
#pragma omp master
    {
      MPI_Start(&(m_request_recv[tag]));
      MPI_Start(&(m_request_send[tag]));
    }

    // bulk
    // site    = rx + Rx*(ry + Ry*(rz + Rz*it) )
    // site_n  = rx + Rx*(ry + Ry*(rz + Rz*(it-1)) )  [ it\neq 0]
    for(int site=is; site<ns; ++site){
      int it=site/(Rx*Ry*Rz);
      if(it==0) {
        continue;
      }
      int site_n= site-Rx*Ry*Rz;
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        vp[i1]=val;
      }
    }

#pragma omp master
    {
      MPI_Wait(&(m_request_send[tag]),&(m_status[tag]));
      MPI_Wait(&(m_request_recv[tag]),&(m_status[tag]));
    }
#pragma omp barrier

    // s_packed= rx + Rx*(rz + Rz*rt) , tmp=rz + Rz*rt
    // site    = rx + Rx*(ry + Ry*(rz + Rz*0) )
    for(int s_packed=is1; s_packed<ns1; ++s_packed){
      int site=s_packed;
      for(int nin=0; nin<size; ++nin){
        int i=s_packed*size+nin;
        int j=site*ANin1+nin;
        vcomplex_t val=rbuf[i];
        apply_boundary(val, boundary2);
        vp[j] = val;
      }
    }
#endif

  } else { // no communication is needed

    int ith, nth, is, ns;
    vreal_t boundary1;
    bcast(boundary1, 1.0);
    set_threadtask(ith, nth, is, ns, Rx*Ry*Rz*Rt);
    for(int site=is; site<ns; ++site){
      vreal_t boundary=boundary1;
      int it=site/(Rx*Ry*Rz);
      int site_n= site-Rx*Ry*Rz;
      if(it==0) {
        boundary=boundary2;
        site_n+=Rx*Ry*Rz*Rt;
      }
      //      int site_n= site+Rx*Ry*Rz(-it + ((it-1+Rt)%Rt));
      for(int nin=0; nin<size; ++nin){
        int i1=site*ANin1+nin;
        int i2=site_n*ANin2+nin;
        vcomplex_t val=wp[i2];
        apply_boundary(val, boundary);
        vp[i1]=val;
      }
    }
  }
#pragma omp barrier
}

//============================================================END=====
