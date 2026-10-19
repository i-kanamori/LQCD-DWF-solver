/*!

        @file    aindex_eo_impl1.h

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */


#ifndef AINDEX_LEX_EO_SIMD_IMPL1_H
#define AINDEX_LEX_EO_SIMD_IMPL1_H


template<class REALTYPE>
class AIndex_eo<REALTYPE, SIMD>{
 protected:
  size_t m_Nx, m_Ny, m_Nz, m_Nt;
  size_t m_Nc, m_Nd, m_Nvol, m_Ndim, m_Rvol;
  size_t m_Nx2;
  size_t m_Nvol2, m_Rvol2;
  Element_type::type m_element_type;
  int VLEN;
  int m_pieo;
  std::vector<size_t> m_ostride;
  std::vector<size_t> m_istride;
  std::vector<size_t> m_rlattice;
  std::vector<int> m_simd_layout;
 public:
  using real_t=REALTYPE;


public:

  AIndex_eo()
  {
    init(CommonParameters::Nx(),  CommonParameters::Ny(),  CommonParameters::Nz(),  CommonParameters::Nt());
  };

  AIndex_eo(Element_type::type elem)
  {
    init(CommonParameters::Nx(),  CommonParameters::Ny(),  CommonParameters::Nz(),  CommonParameters::Nt(), elem);
  }

  AIndex_eo(int Nx, int Ny, int Nz, int Nt)
  {
    init(Nx,Ny,Nz,Nt);
  }

  AIndex_eo(int Nx, int Ny, int Nz, int Nt, Element_type::type elem)
  {
    init(Nx,Ny,Nz,Nt,elem);
  }

  void init(){
    init(CommonParameters::Nx(),  CommonParameters::Ny(),  CommonParameters::Nz(),  CommonParameters::Nt());
  }

  int site(const int x, const int y, const int z, const int t) const
  {
    vout.crucial("site is not tested : %s\n", __FILE__);
    exit(EXIT_FAILURE);

  }

  int Nx() const { return m_Nx; }
  int Ny() const { return m_Ny; }
  int Nz() const { return m_Nz; }
  int Nt() const { return m_Nt; }


 private:



 public:

  void init(int Nx, int Ny, int Nz, int Nt, Element_type::type elem=Element_type::COMPLEX) {

    assert(Nx/2 %2 ==0);

    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
    m_Nvol = Nx*Ny*Nz*Nt;

    m_Nvol2 = m_Nvol/2;
    m_Nx2 = Nx/2;
    m_Ndim = CommonParameters::Ndim();
    m_element_type=elem;

    std::vector<int> llattice(m_Ndim);
    assert(m_Ndim==4);
    llattice[XDIR] = m_Nx2;
    llattice[YDIR] = m_Ny;
    llattice[ZDIR] = m_Nz;
    llattice[TDIR] = m_Nt;
    vout.paranoiac("index: m_Nx= %d\n",m_Nx);
    vout.paranoiac("index: m_Ny= %d\n",m_Ny);
    vout.paranoiac("index: m_Nz= %d\n",m_Nz);
    vout.paranoiac("index: m_Nt= %d\n",m_Nt);

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
    m_Rvol2 = m_Rvol/2;
    m_Ndim = CommonParameters::Ndim();

    vout.paranoiac("index: m_Nx= %d\n",m_Nx);
    vout.paranoiac("index: m_Ny= %d\n",m_Ny);
    vout.paranoiac("index: m_Nz= %d\n",m_Nz);
    vout.paranoiac("index: m_Nt= %d\n",m_Nt);

    m_ostride.resize(m_Ndim);
    m_istride.resize(m_Ndim);
    m_ostride[0]=1;

    for(int mu=0; mu <m_Ndim; ++mu){
      m_rlattice[mu]=llattice[mu]/m_simd_layout[mu];
      if(m_rlattice[mu] == 0 ){
        vout.crucial("Index_eo_alt: too small local lattice: cannot fit the simd varaibles\n");
        vout.crucial("mu=%d, llattice_eo=%d, simd=%d, sizeof(real_t)=%d\n", mu, llattice[mu], m_simd_layout[mu], sizeof(real_t));
        abort();
      }
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

    // even-odd-ness of the origin of the local lattice
    m_pieo=(Communicator::ipe(3)*m_Nt+Communicator::ipe(2)*m_Nz+Communicator::ipe(1)*m_Ny) %2;

  }


  inline size_t ostride(const int mu) const {return m_ostride[mu];}
  inline size_t istride(const int mu) const {return m_istride[mu];}

  size_t scalar_idx(const int& in, const int& Nin,
                    const int& x2, const int& y, const int& z, const int& t,
                    const int& ex) const
  {
    size_t rx2=x2/VLEN;
    size_t ix=x2%VLEN;
    size_t ry=y;
    size_t rz=z;
    size_t rt=t;
    size_t rcoord[4];
    rcoord[0]=x2/VLEN;
    rcoord[1]=y;
    rcoord[2]=z;
    rcoord[3]=t;
    size_t vector_site
      = rcoord[0]*m_ostride[0] + rcoord[1]*m_ostride[1]
      + rcoord[2]*m_ostride[2] + rcoord[3]*m_ostride[3];
    size_t vector_idx = in+Nin*(vector_site + m_Rvol*ex);
    return ix + VLEN*vector_idx;
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


  void vsite_to_rcoord(const size_t &vector_ist, int &rx, int &ry, int &rz, int &rt) const
  {
    // vector_site (=site index to specify vector variable)
    // to reduced coordinate (=coordiante for vector vaiables)
    size_t idx=vector_ist;
    rx = idx % m_rlattice[XDIR];
    idx = idx / m_rlattice[XDIR];

    ry = idx % m_rlattice[YDIR];
    idx = idx / m_rlattice[YDIR];

    rz = idx % m_rlattice[ZDIR];
    idx = idx / m_rlattice[ZDIR];

    rt = idx % m_rlattice[TDIR];
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


  template <typename AFIELD>
  void split(AFIELD& field_e, AFIELD& field_o,
             const AFIELD& field_lex);

  template <typename AFIELD>
  void split_complex(AFIELD& field_e, AFIELD& field_o,
                     const AFIELD& field_lex);

  template <typename AFIELD>
  void merge(AFIELD& field_lex,
             const AFIELD& field_e, const AFIELD& field_o);

  template <typename AFIELD>
  void merge_complex(AFIELD& field_lex,
                     const AFIELD& field_e, const AFIELD& field_o);


};
#endif
