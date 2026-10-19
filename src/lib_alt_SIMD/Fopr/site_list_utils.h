#ifndef SITE_LIST_UTIL_H
#define SITE_LIST_UTIL_H
/*!

        @file    $Id: site_list_util.h #$

        @brief   site index utility class

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */
#include  "Parameters/commonParameters.h"
#include  "IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_SIMD/inline/valtype.h"
#include "site_list.h"
#include "mpi_bools.h"

/*
  It was originally meant for reduced lattice with simd variables,
  but does not use any information of the simd size or local lattice
  size from CommonParameters.
  Therefore it works for the normal lattice as well as coarse grid
  lattice provided the target lattice size is given.

 */

class site_list_utils{

public:
  using vector_site_t = aligned_vector<site_list_t>;
  using vector_pack_t = aligned_vector<site_list_pack_t>;
  const std::string class_name="site_list_utils";

  site_list_utils() : m_dist0(0), m_dist1(0), m_vl(CommonParameters::Vlevel())
  {  }

  site_list_utils(const std::vector<int> &Rlattice, const mpi_bools &bools) : m_dist0(0), m_dist1(0), m_vl(CommonParameters::Vlevel())  {
    init(Rlattice, bools);
  }

  void init(const std::vector<int> &Rlattice, const mpi_bools &bools);

  void init(const std::vector<int> &Rlattice, const mpi_bools &bools,
            const int prefetch_dist0, const int prefetch_dist1){
    init(Rlattice, bools);
    set_prefetch_dist(prefetch_dist0, prefetch_dist1);
  }

  void set_vl(Bridge::VerboseLevel vl) { m_vl = vl; }

  void set_prefetch_dist(const int prefetch_dist0, const int prefetch_dist1){
    m_dist0=prefetch_dist0;
    m_dist1=prefetch_dist1;
  }

  void set_site_lex_bulk(vector_site_t &site_list_bulk) const;

  void set_site_lex_boundary(vector_site_t &site_list_boundary,
                             vector_pack_t &pack_list_x,
                             vector_pack_t &pack_list_y,
                             vector_pack_t &pack_list_z,
                             vector_pack_t &pack_list_t,
                             const vector_site_t &site_list_bulk) const;

  void set_site_eo_bulk(vector_site_t *site_list_bulk_ptr, const int eo_offset) const;

  void set_site_eo_boundary(vector_site_t *site_list_bounary,
                            vector_pack_t *pack_list_x,
                            vector_pack_t *pack_list_y,
                            vector_pack_t *pack_list_z,
                            vector_pack_t *pack_list_t,
                            const vector_site_t *site_list_bulk,
                            const int eo_offset) const;

private:
  int m_Rx, m_Ry, m_Rz, m_Rt;
  int m_Rst, m_Rst_bulk, m_Rst_boundary;
  mpi_bools m_bools;
  int m_dist0;
  int m_dist1;

protected:
  Bridge::VerboseLevel m_vl;

private:
  site_list_utils(const site_list_utils&);

};

#endif
