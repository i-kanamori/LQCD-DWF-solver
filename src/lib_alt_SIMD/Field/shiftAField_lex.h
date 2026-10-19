/*!
        @file    $Id:: shiftAField_lex.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2017-12-21 10:39:00 #$

        @version $LastChangedRevision: 1687 $
*/

#ifndef SHIFTAFIELD_LEX_SIMD_INCLUDED
#define SHIFTAFIELD_LEX_SIMD_INCLUDED

#include "lib_alt_SIMD/alt_simd_def.h"

#ifdef SHIFT_FORCE_COMMNUICATE
#undef SHIFT_FORCE_COMMNUICATE
#endif
#ifdef FORCE_COMMNNICATE
#define SHIFT_FORCE_COMMNUICATE
#endif

#ifdef SHIFT_EXCLUSIVE_THREAD_FOR_COMM
#undef SHIFT_EXCLUSIVE_THREAD_FOR_COMM
#endif
#ifdef EXCLUSIVE_THREAD_FOR_COMM
#define SHIFT_EXCLUSIVE_THREAD_FOR_COMM
#endif

#include <vector>

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "aindex_lex.h"

//! Methods to shift a field in the lexical site index.

/*!
   This template class is an alternative version of ShiftField
   in the Bridge core libirary.
   This implementation assumes the complex array which is
   integrated in the inner-most SIMD vector.
   Differently from the Field version, the Nin size must be
   specified at the construction.
   Note that the Nin/2 is placed as the size of second degree
   of freedom if the field array.
   The boundary condition is also given at the construction.
                                   [25 Dec 2017 H.Matsufuru]
 */
template<typename AFIELD>
class ShiftAField_lex {

 public:
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;
  typedef typename AFIELD::vreal_t vreal_t;
  static const Impl IMPL = AFIELD::IMPL;
  typedef AIndex_lex<real_t, IMPL> index_t;

  static const std::string class_name;

 private:
  int m_Nin;             //!< internal degree of freedom.
  int m_Nvol;
  int m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nxv, m_Nstv;

  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  //std::vector<int> m_Nbdsize;
  //  std::vector<Channel*> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  //  ChannelSet chset_send, chset_recv;

#ifdef USE_MPI

    // buffer: body
  //   TODO: replace with the chanel
  aligned_vector<vcomplex_t> buf_up_body;  // for sending
  aligned_vector<vcomplex_t> rbuf_up_body; // for receiving
  aligned_vector<vcomplex_t> buf_dn_body;  // for sending
  aligned_vector<vcomplex_t> rbuf_dn_body; // for receiving

  vcomplex_t *buf_up[4], *rbuf_up[4], *buf_dn[4], *rbuf_dn[4];
  vcomplex_t *buf_tag[8], *rbuf_tag[8];

  // nearest neighbor nodes
  std::vector<int> ipe_neighbor;

  std::vector<MPI_Request> m_request_send;  //!< a collection of MPI request.
  std::vector<MPI_Request> m_request_recv;  //!< a collection of MPI request.
  std::vector<MPI_Status> m_status;  //!< a collection of MPI status.
#endif

 public:
  ShiftAField_lex(int nin)
    {
      init(CommonParameters::Nx(),
           CommonParameters::Ny(),
           CommonParameters::Nz(),
           CommonParameters::Nt(), nin);
    }

  ShiftAField_lex(int nin, const std::vector<int>& bc)
    {
      init(CommonParameters::Nx(),
           CommonParameters::Ny(),
           CommonParameters::Nz(),
           CommonParameters::Nt(), nin, bc);
    }

  ShiftAField_lex(int Nx, int Ny, int Nz, int Nt, int nin)
    { init(Nx,Ny,Nz,Nt,nin); }

  ShiftAField_lex(int Nx, int Ny, int Nz, int Nt, int nin, const std::vector<int>& bc)
    { init(Nx,Ny,Nz,Nt,nin, bc); }

 private:
  // non-copyable
  ShiftAField_lex(const ShiftAField_lex&);
  ShiftAField_lex& operator=(const ShiftAField_lex&);

 public:
  void forward(AFIELD&, const AFIELD&, const int mu);
  void backward(AFIELD&, const AFIELD&, const int mu);

  void forward( AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);
  void backward(AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);

  void forward(AFIELD&, const int in_start1, const int ex1,
               const AFIELD&, const int in_start2, const int ex2, const int size,
               const int mu);
  void backward(AFIELD&, const int in_start1, const int ex1,
                const AFIELD&, const int in_start2, const int ex2, const int size,
                const int mu);

  void init(int Nx, int Ny, int Nz, int Nt, int Nin);
  void init(int Nx, int Ny, int Nz, int Nt, int Nin, const std::vector<int>& bc);

 private:

  //! setup channels for communication.
  void setup_channels();

  void setup_comm(const std::vector<size_t> &buf_size);


  void up_x(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void up_y(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void up_z(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void up_t(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void dn_x(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void dn_y(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void dn_z(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);
  void dn_t(vcomplex_t *, const vcomplex_t *, const int ANin1 ,const int ANin2, const int size, const index_t &index);

};

#endif
