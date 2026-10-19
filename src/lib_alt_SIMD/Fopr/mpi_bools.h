/*!
        @file    $Id:: mpi_bool.h#$

        @brief

        @author  <Issaku Kanamori> kanamori@hiroshima-u.ac.jp
                 $LastChangedBy:  $

        @date    $LastChangedDate::  #$

        @version $LastChangedRevision:  $
*/
#ifndef AFOPR_MPI_BOOLS
#define AFOPR_MPI_BOOLS

struct mpi_bools{
  bool has_mpi_x;
  bool has_mpi_y;
  bool has_mpi_z;
  bool has_mpi_t;
};

namespace{
  void print_mpi_bools(const mpi_bools &bools, const Bridge::VerboseLevel &m_vl){
#ifdef USE_MPI
    if(bools.has_mpi_x){
      vout.general(m_vl, "has_mpi_x: yes\n");
    } else {
      vout.general(m_vl, "has_mpi_x: no\n");
    }
    if(bools.has_mpi_y){
      vout.general(m_vl, "has_mpi_y: yes\n");
    } else {
      vout.general(m_vl, "has_mpi_y: no\n");
    }
    if(bools.has_mpi_z){
      vout.general(m_vl, "has_mpi_z: yes\n");
    } else {
      vout.general(m_vl, "has_mpi_z: no\n");
    }
    if(bools.has_mpi_t){
      vout.general(m_vl, "has_mpi_t: yes\n");
    } else {
      vout.general(m_vl, "has_mpi_t: no\n");
    }
#else
    vout.general(m_vl, "has_mpi: USE_MPI is off\n");
#endif
  }
}
#endif
