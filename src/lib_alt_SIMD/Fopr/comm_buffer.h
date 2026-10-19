#ifndef COMM_BUFFER_H
#define COMM_BUFFER_H
/*!

        @file    $Id: comm_buffer.h #$

        @brief   auxially class to initialize communication buffer

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */

#ifdef USE_MPI
#include <mpi.h>
#endif


#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
#include "lib_alt_SIMD/alt_simd_def.h"

using Bridge::vout;

template<typename ValType>
class CommBuffer{
 public:
  using val_t = ValType;
  using buffer_t = aligned_vector<ValType>;
  static const std::string class_name;

  CommBuffer() : m_num(0), m_vl(CommonParameters::Vlevel()){
    init(1);  // 8 is for 4 dim x (forward, backward)
  }

  //! initilize: [num_hop = 1: nearest neighbor, 2: include next nearest]
  void init(const int num_hop);

  //! verbose level
  void set_vl(Bridge::VerboseLevel vl) { m_vl = vl; }

  //! allocate buffer and initiate persistent communication
  void append(val_t **rbuf, const int rdir,
              val_t **sbuf, const int sdir,
              const int size, const int tag);



  //! start communication (receive)
  void start_receiving(){
#ifdef USE_MPI
#pragma omp master
  {
    MPI_Startall(m_num, &m_rrequest[0]);
  }
#endif
  }

  //! start communication (receive)
  void start_receiving(int rdir){
#ifdef USE_MPI
#pragma omp master
  {
    int id=m_request_id[rdir];
    if(id>=0){
      MPI_Start(&m_rrequest[id]);
    }
  }
#endif
  }


  //! start communication (send)
  void start_sending(){
#ifdef USE_MPI
#pragma omp master
  {
    MPI_Startall(m_num, &m_srequest[0]);
  }
#endif
  }

  //! start communication (send)
  void start_sending(int rdir){
#ifdef USE_MPI
#pragma omp master
  {
    int id=m_request_id[rdir];
    if(id>=0){
      MPI_Start(&m_srequest[id]);
    }
  }
#endif
  }


  //! call MPI_Testall
  void test_all(){
#ifdef USE_MPI
#pragma omp master
    {
    int complete1, complete2;
    MPI_Testall(m_num, &m_rrequest[0], &complete1, MPI_STATUSES_IGNORE);
    MPI_Testall(m_num, &m_srequest[0], &complete2, MPI_STATUSES_IGNORE);
    }
#endif
  }

  void test_send(){
#ifdef USE_MPI
#pragma omp master
    {
      int complete1;
      MPI_Testall(m_num, &m_srequest[0], &complete1, MPI_STATUSES_IGNORE);

    }
#endif
  }

  //! call MPI_Testall
  void test_all(int rdir){
#ifdef USE_MPI
#pragma omp master
    {
      int id=m_request_id[rdir];
      if(id >= 0){
        int complete1, complete2;
        MPI_Test(&m_rrequest[id], &complete1, MPI_STATUSES_IGNORE);
        MPI_Test(&m_srequest[id], &complete2, MPI_STATUSES_IGNORE);
      }
    }
#endif
  }

  //! call MPI_Testall
  void test_send(int rdir){
#ifdef USE_MPI
#pragma omp master
    {
      int id=m_request_id[rdir];
      if(id >= 0){
        int complete1, complete2;
        MPI_Test(&m_srequest[id], &complete1, MPI_STATUSES_IGNORE);
      }
    }
#endif
  }


  //! MPI_Waitall
  void wait_receiving(){
#ifdef USE_MPI
#pragma omp master
    {
      MPI_Waitall(m_num, &m_rrequest[0], MPI_STATUSES_IGNORE);
    }
#endif
  }

  //! MPI_Waitall
  void wait_receiving(int rdir){
#ifdef USE_MPI
#pragma omp master
    {
      int id=m_request_id[rdir];
      if(id >= 0){
        MPI_Wait(&m_rrequest[id], MPI_STATUSES_IGNORE);
      }

    }
#endif
  }


  //! MPI_Waitall
  void wait_sending(){
#ifdef USE_MPI
#pragma omp master
    {
      MPI_Waitall(m_num, &m_srequest[0], MPI_STATUSES_IGNORE);
    }
#endif
  }

  //! MPI_Waitall
  void wait_sending(int rdir){
#ifdef USE_MPI
#pragma omp master
    {
      int id=m_request_id[rdir];
      if(id >= 0){
        MPI_Wait(&m_srequest[id], MPI_STATUSES_IGNORE);
      }
    }
#endif
  }


  //#ifdef USE_MPI
  //  void get_rrequest(MPI_Request *req, int n) {req=&m_rrequest[n]; }
  //  void get_srequest(MPI_Request *req, int n) {req=&m_srequest[n]; }
  //#endif

protected:
  int m_num;
  Bridge::VerboseLevel m_vl;

 private:
  CommBuffer(const CommBuffer&);

  void init_nn();  // nearest neighbor
  void init_nnn(); // next nearest neighbor

  //! mpi grid
  std::vector<int> m_ipe_neighbor;

  //! body of the buffers
  std::vector<buffer_t> m_rbuf_array;
  std::vector<buffer_t> m_sbuf_array;

#ifdef USE_MPI
  //! MPI request
  aligned_vector<MPI_Request> m_rrequest;
  aligned_vector<MPI_Request> m_srequest;

  //! map to request id from rdir
  std::vector<int> m_request_id;
#endif

};

#endif
