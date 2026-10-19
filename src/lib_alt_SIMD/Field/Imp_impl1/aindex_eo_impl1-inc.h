/*!
        @file    $Id:: aindex_eo_impl1-inc.h #$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  #$

        @version $LastChangedRevision: $
*/


#include <assert.h>

//====================================================================


template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,SIMD>::split(AFIELD& field_e, AFIELD& field_o,
                                     const AFIELD& field_lex)
{
  typedef REALTYPE real_t;
  using vreal_t=ValType_real<real_t>;
  assert(AFIELD::IMPL==SIMD);
  assert(sizeof(typename AFIELD::real_t) == sizeof(real_t));

  Element_type::type etype;
  etype = field_lex.field_element_type();

  {
    int Nin = field_lex.nin();
    int Nex = field_lex.nex();

    size_t Nvol = field_lex.nvol();
    size_t Nvol2 = Nvol/2;

    vout.paranoiac("split:  m_Nvol=%d\n",m_Nvol);
    vout.paranoiac("split:   Nvol2=%d\n",Nvol2);

    assert(Nvol == m_Nvol);
    assert(Nvol2 == m_Nvol2);
    assert(field_e.check_size(Nin, Nvol2, Nex));
    assert(field_o.check_size(Nin, Nvol2, Nex));
    assert(field_e.field_element_type()==etype);
    assert(field_o.field_element_type()==etype);
  }
  if(etype== Element_type::REAL){
    vout.crucial("lex --> evenodd: only complex field is supported\n");
    exit(EXIT_FAILURE);
  } else {
    split_complex(field_e, field_o, field_lex);
  }

}

template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,SIMD>::split_complex(AFIELD& field_e, AFIELD& field_o,
                                             const AFIELD& field_lex)
{
  typedef REALTYPE real_t;
  using complex_t=typename AFIELD::complex_t;
  assert(AFIELD::IMPL==SIMD);
  assert(sizeof(typename AFIELD::real_t) == sizeof(real_t));

  int Nin = field_lex.nin_complex();
  int Nex = field_lex.nex();

  AIndex_lex<REALTYPE,SIMD> index_lex;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  complex_t *out_e=field_e.ptr_complex();
  complex_t *out_o=field_o.ptr_complex();
  const complex_t *in_lex=field_lex.ptr_complex();

  const int Nx2=m_Nx2;
  const int proc_ieo = m_pieo;
  for(int ex = 0; ex < Nex; ++ex){
    for(size_t ist2 = is; ist2 < ns; ++ist2){

      int tmp=ist2;
      const int x2 = tmp % Nx2;
      tmp /= Nx2;
      const int y = tmp % m_Ny;
      tmp/=m_Ny;
      const int z = tmp % m_Nz;
      const int t = tmp / m_Nz;

      { //even
        const int ieo=0;
        const int x=x2*2+(y+z+t+proc_ieo+ieo)%2;
        for(int in = 0; in < Nin; ++in){
          size_t idx_lex=index_lex.scalar_idx(in, Nin, x,y,z,t, ex);
          size_t idx_eo=scalar_idx(in,Nin,x2,y,z,t,ex);
          out_e[idx_eo]=in_lex[idx_lex];
        }
      }
      { //odd
        const int ieo=1;
        int x=x2*2+(y+z+t+proc_ieo+ieo)%2;
        for(int in = 0; in < Nin; ++in){
          size_t idx_lex=index_lex.scalar_idx(in, Nin, x,y,z,t, ex);
          size_t idx_eo=scalar_idx(in,Nin,x2,y,z,t,ex);
          out_o[idx_eo]=in_lex[idx_lex];
        }
      }
    }
  }
} // split


//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,SIMD>::merge(AFIELD& field_lex,
                                     const AFIELD& field_e,
                                     const AFIELD& field_o)
{

  typedef REALTYPE real_t;
  assert(AFIELD::IMPL==SIMD);
  assert(sizeof(typename AFIELD::real_t) == sizeof(real_t));

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();
  Element_type::type etype;
  etype = field_lex.field_element_type();

  size_t Nvol = field_lex.nvol();
  size_t Nvol2 = Nvol/2;

  assert(Nvol2 == m_Nvol2);
  assert(Nvol == m_Nvol);
  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));
  assert(field_e.field_element_type()==etype);
  assert(field_o.field_element_type()==etype);

  if(etype== Element_type::REAL){
    vout.crucial("evenodd --> lex: only complex field is supported\n");
    exit(EXIT_FAILURE);
  } else {
    merge_complex(field_lex, field_e, field_o);
  }
} // merge



//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,SIMD>::merge_complex(AFIELD& field_lex,
                                             const AFIELD& field_e,
                                             const AFIELD& field_o)
{
  typedef REALTYPE real_t;
  typedef typename AFIELD::complex_t complex_t;
  assert(AFIELD::IMPL==SIMD);
  assert(sizeof(typename AFIELD::real_t) == sizeof(real_t));


  int Nin = field_lex.nin_complex();
  int Nex = field_lex.nex();
  AIndex_lex<REALTYPE,SIMD> index_lex;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);
  complex_t *out_lex=field_lex.ptr_complex();
  const complex_t *in_e=field_e.ptr_complex();
  const complex_t *in_o=field_o.ptr_complex();

  const int Nx2=m_Nx2;
  const int proc_ieo = m_pieo;

  for(int ex = 0; ex < Nex; ++ex){
    for(size_t ist2 = is; ist2 < ns; ++ist2){

      int tmp=ist2;
      int x2 = tmp % Nx2;
      tmp /= Nx2;
      int y = tmp % m_Ny;
      tmp/=m_Ny;
      int z = tmp % m_Nz;
      int t=tmp / m_Nz;

      { // even
        const int ieo=0;
        int x=x2*2+(y+z+t+proc_ieo+ieo)%2;
        for(int in = 0; in < Nin; ++in){
          size_t idx_lex=index_lex.scalar_idx(in, Nin, x,y,z,t, ex);
          size_t idx_eo=scalar_idx(in,Nin,x2,y,z,t,ex);
          out_lex[idx_lex]=in_e[idx_eo];
        }
      }

      { // odd
        const int ieo=1;
        int x=x2*2+(y+z+t+proc_ieo+ieo)%2;
        for(int in = 0; in < Nin; ++in){
          size_t idx_lex=index_lex.scalar_idx(in, Nin, x,y,z,t, ex);
          size_t idx_eo=scalar_idx(in,Nin,x2,y,z,t,ex);
          out_lex[idx_lex]=in_o[idx_eo];
        }
      }

    } // ist loop
  } // ex loop
} // merge_complex


//====================================================================
//============================================================END=====
