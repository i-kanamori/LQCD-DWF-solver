/*!
        @file    aindex_lex_impl1.h

        @brief   index implementation for impl1

        @author  <Issaku Kanamori> kanamori@hiroshima-u.ac.jp
                 $LastChangedBy:  $

        @date    $LastChangedDate::  #$

        @version $LastChangedRevision:  $
*/
//====================================================================

#ifndef AINDEX_LEX_SIMD_IMPL1_INC_H
#define AINDEX_LEX_SIMD_IMPL1_INC_H

#include <assert.h>

//====================================================================
template<class REALTYPE>
class AIndex_lex<REALTYPE, SIMD>{
 protected:
  size_t m_Nx, m_Ny, m_Nz, m_Nt;
  size_t m_Nc, m_Nd, m_Nvol, m_Ndim, m_Rvol;
  Element_type::type m_element_type;
  int VLEN;
  std::vector<size_t> m_ostride;
  std::vector<size_t> m_istride;
  std::vector<size_t> m_rlattice;
  std::vector<int> m_simd_layout;
public:
  using real_t=REALTYPE;

public:

  AIndex_lex()
  {
    init(CommonParameters::Nx(),  CommonParameters::Ny(),  CommonParameters::Nz(),  CommonParameters::Nt());
  }

  AIndex_lex(Element_type::type elem)
  {
    init(CommonParameters::Nx(),  CommonParameters::Ny(),  CommonParameters::Nz(),  CommonParameters::Nt(), elem);
  }

  AIndex_lex(int Nx, int Ny, int Nz, int Nt)
  {
    init(Nx,Ny,Nz,Nt);
  }

  AIndex_lex(int Nx, int Ny, int Nz, int Nt, Element_type::type elem)
  {
    init(Nx,Ny,Nz,Nt,elem);
  }

  int site(const int x, const int y, const int z, const int t) const
  {  return m_Nx * (m_Ny * (m_Nz * t + z) + y) + x;  }

  int Nx() const { return m_Nx; }
  int Ny() const { return m_Ny; }
  int Nz() const { return m_Nz; }
  int Nt() const { return m_Nt; }

private:

public:

  void init(int Nx, int Ny, int Nz, int Nt, Element_type::type elem=Element_type::COMPLEX) {
    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
    m_Nvol = Nx*Ny*Nz*Nt;
    m_element_type=elem;
    m_Ndim = CommonParameters::Ndim();


    std::vector<int> llattice(m_Ndim);
    assert(m_Ndim==4);
    llattice[XDIR] = m_Nx;
    llattice[YDIR] = m_Ny;
    llattice[ZDIR] = m_Nz;
    llattice[TDIR] = m_Nt;
    vout.paranoiac("index: m_Nx= %d\n",m_Nx);
    vout.paranoiac("index: m_Ny= %d\n",m_Ny);
    vout.paranoiac("index: m_Nz= %d\n",m_Nz);
    vout.paranoiac("index: m_Nt= %d\n",m_Nt);

    m_ostride.resize(m_Ndim);
    m_istride.resize(m_Ndim);
    m_simd_layout.resize(m_Ndim);
    m_rlattice.resize(m_Ndim);

    if(m_element_type == Element_type::COMPLEX){      // complex field
      simd_layout_complex<real_t> simd_layout;
      VLEN=ValType_complex<real_t>::Nsimd();
      for(int mu=0; mu <m_Ndim; ++mu){
        m_simd_layout[mu]=simd_layout.vlen(mu);
      }
    } else {                          // real field
      simd_layout_real<real_t> simd_layout;
      VLEN=ValType_real<real_t>::Nsimd();
      for(int mu=0; mu <m_Ndim; ++mu){
        m_simd_layout[mu]=simd_layout.vlen(mu);
      }
    }

    m_Rvol = m_Nvol/VLEN;
    for(int mu=0; mu <m_Ndim; ++mu){
      m_rlattice[mu]=llattice[mu]/m_simd_layout[mu];
    }

    m_ostride[0]=1;
    for(int mu=1; mu<m_Ndim; ++mu){
      m_ostride[mu]=m_ostride[mu-1]*m_rlattice[mu-1];
    }
    vout.paranoiac("index: x: m_ostride=%d, rlattice=%d\n",m_ostride[0], m_rlattice[0]);
    vout.paranoiac("index: y: m_ostride=%d, rlattice=%d\n",m_ostride[1], m_rlattice[1]);
    vout.paranoiac("index: z: m_ostride=%d, rlattice=%d\n",m_ostride[2], m_rlattice[2]);
    vout.paranoiac("index: t: m_ostride=%d, rlattice=%d\n",m_ostride[3], m_rlattice[3]);


    for(int mu=0; mu<m_Ndim; ++mu){
      if(m_simd_layout[mu]==1){
        m_istride[mu]=1;
      } else {
        m_istride[mu]=0;
      }
    }

#pragma omp barrier
  }


  inline size_t ostride(const int mu) const {return m_ostride[mu];}
  inline size_t istride(const int mu) const {return m_istride[mu];}

  size_t vector_site(const int& x, const int& y, const int& z, const int& t) const
  {
    // site index to specify vector variable
    //   it does not count internal/outer d.o.f
    // see oIndex() in cartesian/Cartesian_base.h in GRID
    size_t idx = x / m_simd_layout[0];
    idx+=m_ostride[1]*y;
    idx+=m_ostride[2]*z;
    idx+=m_ostride[3]*t;
    return idx;
  }

  size_t rcoord_to_vector_site(const int& rx, const int& ry, const int& rz, const int& rt) const
  {
    // site index to specify vector variable (from reduced coordinate)
    size_t idx=m_ostride[0]*(rx);
    idx+=m_ostride[1]*(ry);
    idx+=m_ostride[2]*(rz);
    idx+=m_ostride[3]*(rt);
    return idx;
  }

  size_t iIndex(const int& x, const int& y, const int& z, const int& t) const
  {
    // index to specify the lane in vector variable
    return (x % m_simd_layout[0]);
  }


  size_t vector_idx(const int& in, const int& Nin,
                    const int& x, const int& y, const int& z, const int& t,
                    const int& ex) const
  {
    // index to access as an array of vector variables
    return in +  Nin*(vector_site(x,y,z,t) + m_Rvol*ex);
  }

  size_t vector_idx(const int& in, const int& Nin,
                    const int& vsite,
                    const int& ex) const
  {
    // index to access as an array of vector variables
    return in +  Nin*(vsite + m_Rvol*ex);
  }

  size_t scalar_idx(const int& in, const int& Nin,
                    const int& x, const int& y, const int& z, const int& t,
                    const int& ex) const
  {
    // index to access as an array of scalar variables
    return iIndex(x,y,z,t) + VLEN*vector_idx(in, Nin, x, y, z, t, ex);
  }

  size_t scalar_idx(const int& in, const int& Nin,
                    const int& site,
                    const int& ex) const
  {
    size_t vector_site = site/VLEN;
    size_t ix = site%VLEN;
    size_t vector_idx = in+Nin*(vector_site + m_Rvol*ex);
    return ix + VLEN*vector_idx;
  }

  size_t idx_SPr(const int ic, const int id, const int site, const int ex) const {
    int in = 2*(ic + id*m_Nc);
    return scalar_idx(in, 2*m_Nc*m_Nd, site, ex);
  }

  size_t idx_SPi(const int ic, const int id, const int site, const int ex) const {
    int in = 2*(ic + id*m_Nc) + 1;
    return scalar_idx(in, 2*m_Nc*m_Nd, site, ex);
  }

  void vsite_to_rcoord(const size_t &vector_ist, int &rx, int &ry, int &rz, int &rt) const
  {
    // vector_site (=site index to specify vector variable)
    // to reduced coordinate (=coordiante for vector vaiables)
    size_t idx=vector_ist;
    rx = idx % m_rlattice[0];
    idx = idx / m_rlattice[0];

    ry = idx % m_rlattice[1];
    idx = idx / m_rlattice[1];

    rz = idx % m_rlattice[2];
    idx = idx / m_rlattice[2];

    rt = idx % m_rlattice[3];

  }


  void vsite_to_coord(const size_t &vector_ist, const size_t &iidx, int &x, int &y, int &z, int &t) const
  {
    int rx;
    vsite_to_rcoord(vector_ist, rx, y, z, t);
    x=iidx+VLEN*rx;
  }

  inline void vsite_to_site(const size_t &vector_ist, std::vector<size_t> &site_array) const
  {
    const size_t idx=vector_ist*VLEN;
    for(int i=0; i<VLEN; ++i){
      site_array[i]=idx+i;
    }
  }

  inline int rlattice(int mu) const {return m_rlattice[mu];}
};

#endif // include guard
