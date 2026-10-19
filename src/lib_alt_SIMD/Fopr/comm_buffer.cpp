/*!

        @file    $Id: comm_buffer.cpp #$

        @brief   auxially class to initialize communication buffer

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */
#include <assert.h>

#include "comm_buffer.h"
#include "lib_alt_SIMD/alt_simd_def.h"
#include "lib_alt_SIMD/inline/vector_arith-inc.h"
#include "lib/ResourceManager/threadManager.h"
#ifdef USE_MPI
#include "lib/Communicator/MPI/communicator_mpi.h"
#endif


namespace{
  void  set_zero(complexD_t &zero){
    zero=complexD_t(0.0,0.0);
  };
  void  set_zero(complexF_t &zero){
    zero=complexF_t(0.0,0.0);
  };
}

template<typename ValType>
void CommBuffer<ValType>::init(const int num_hop) {

  if(num_hop>2 || num_hop<1){
    vout.crucial(m_vl, "%s: bad num_hop\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if(num_hop>0){
    init_nn(); // nearest neighbor
  }
  if(num_hop>1){
    init_nn(); // next nearest neighbor
  }
  vout.general(m_vl, "%s: initialization, done for num_hup=%d\n", class_name.c_str(), num_hop);
}

template<typename ValType>
void CommBuffer<ValType>::init_nn() {

#ifdef USE_MPI
  ThreadManager_OpenMP::assert_single_thread(class_name);

  // array of buffers
  int myrank=Communicator::self();
  int Ndim=CommonParameters::Ndim();
  m_ipe_neighbor.resize(2*Ndim);
  aligned_vector<int> coord(Ndim);

  m_rbuf_array.resize(2*Ndim);
  m_sbuf_array.resize(2*Ndim);
  m_rrequest.resize(2*Ndim);
  m_srequest.resize(2*Ndim);
  m_request_id.resize(2*Ndim);

  for (int i = 0; i < Ndim; ++i) {
    Communicator::grid_coord(&coord[0], Communicator::self());

    int coord_keep=coord[i];
    // upward
    coord[i] = (coord_keep + 1) % Communicator::npe(i);
    Communicator::grid_rank(&m_ipe_neighbor[2*i]  ,&coord[0]);
    m_request_id[2*i]=-1;

    // downward
    coord[i] = (coord_keep - 1 + Communicator::npe(i)) % Communicator::npe(i);
    Communicator::grid_rank(&m_ipe_neighbor[2*i+1]  ,&coord[0]);
    m_request_id[2*i+1]=-1;

#ifdef DEBUG
    {
      int nodeid=Communicator::nodeid();
      printf("hoge: nideid=%d: i=%d, npe(i)=%d, ipe(i)=%d, rank_up=%d, rank_down=%d\n",nodeid,i,
             Communicator::npe(i),
             Communicator::ipe(i),
             m_ipe_neighbor[2*i], m_ipe_neighbor[2*i+1]);
    }
#endif
  }
#endif // USE_MPI

}

template<typename ValType>
void CommBuffer<ValType>::init_nnn() {
  vout.crucial(m_vl, "%s: next nearest neighbering is not ready\n", class_name.c_str());
    exit(EXIT_FAILURE);
}

template<typename ValType>
void CommBuffer<ValType>::append(val_t **rbuf, const int rdir,
                                 val_t **sbuf, const int sdir,
                                 const int size, const int tag){
#ifdef USE_MPI
  ThreadManager_OpenMP::assert_single_thread(class_name);

#ifdef DEBUG
  vout.general("%s:  append, num=%d, size=%d\n", class_name.c_str(), m_num, size);
#endif

  m_rbuf_array[m_num].resize(size);
  m_sbuf_array[m_num].resize(size);
  *rbuf=&(m_rbuf_array[m_num][0]);
  *sbuf=&(m_sbuf_array[m_num][0]);

  // first touch
#ifdef AFOPR_FIRST_TOUCH_BUF
  val_t vzero;
  set_zero(vzero);
  for(int i=0; i<size; i++) { (*rbuf)[i]=vzero; }
  for(int i=0; i<size; i++) { (*sbuf)[i]=vzero; }
#endif

  int to   = m_ipe_neighbor[sdir];
  int from = m_ipe_neighbor[rdir];

   MPI_Recv_init(*rbuf, sizeof(val_t)*size, MPI_BYTE, from, tag, Communicator_impl::world(), &m_rrequest[m_num]);
   MPI_Send_init(*sbuf, sizeof(val_t)*size, MPI_BYTE, to,   tag, Communicator_impl::world(), &m_srequest[m_num]);
   m_request_id[rdir]=m_num;
   m_num++;

#endif
}



// explicit instanciation.
template class CommBuffer<vcomplexF_t>;
template<>
const std::string CommBuffer<vcomplexF_t>::class_name = "CommBuffer<vcomplexF_t>";

template class CommBuffer<vcomplexD_t>;
template<>
const std::string CommBuffer<vcomplexD_t>::class_name = "CommBuffer<vcomplexD_t>";

template class CommBuffer<complexF_t>;
template<>
const std::string CommBuffer<complexF_t>::class_name = "CommBuffer<complexF_t>";

template class CommBuffer<complexD_t>;
template<>
const std::string CommBuffer<complexD_t>::class_name = "CommBuffer<complexD_t>";
