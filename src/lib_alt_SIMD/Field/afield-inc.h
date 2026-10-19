// This file must not be included in .h file.
/*!
 Time-stamp: <2023-10-20 04:16:01 kanamori>

        @file    $Id: afield-inc.h #$

        @brief   useful functions to access AField

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */


#ifndef AFIELD_SIMD_INC_INCLUDED
#define AFIELD_SIMD_INC_INCLUDED

#include "lib/ResourceManager/threadManager.h"
#include "../inline/vector_arith-inc.h"

//====================================================================


//====================================================================
namespace {
  inline
  void set_threadtask(int& ith, int& nth, int& is, int& ns,
                      const int size)
  {
    nth = ThreadManager_OpenMP::get_num_threads();
    ith = ThreadManager_OpenMP::get_thread_id();

    size_t tmp_is = (size_t)size * ith;
    size_t tmp_ns = (size_t)size * (ith+1);
    is = tmp_is/nth;
    ns = tmp_ns/nth;
  }

  inline void set_threadtask_without_master_afield(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    nth = ThreadManager_OpenMP::get_num_threads();
    ith = ThreadManager_OpenMP::get_thread_id();

    if(nth==1){
      set_threadtask(ith,nth,is,ns,size);
      return;
    }
    if(ith==0){ // assumes that the 0-th thread is the master
      is=0;
      ns=0;
    } else{
      size_t tmp_is = (size_t)size * (ith-1);
      size_t tmp_ns = (size_t)size * (ith);
      is = tmp_is / (nth-1);
      ns = tmp_ns / (nth-1);
    }
  }
}



//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE,SIMD>& v, const AField<REALTYPE,SIMD> &w)
{
  v.copy(w);
}

//====================================================================
template <typename REALTYPE>
void copy(AField<REALTYPE,SIMD>& v, const int ex,
          const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  v.copy(ex, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
REALTYPE norm2(const AField<REALTYPE,SIMD>& v)
{
  return v.norm2();
}


//====================================================================
template <typename REALTYPE>
REALTYPE dot(const AField<REALTYPE,SIMD>& v, const AField<REALTYPE,SIMD>& w)
{
  return v.dot(w);
}


//====================================================================
template <typename REALTYPE>
typename AField<REALTYPE,SIMD>::complex_t dotc(const AField<REALTYPE,SIMD>& v, const AField<REALTYPE,SIMD>& w)
{
  return v.dotc(w);  // <v|w>
}

//====================================================================
template <typename REALTYPE>
typename AField<REALTYPE,SIMD>::complex_t dotc(const AField<REALTYPE,SIMD>& v, const int ex, const AField<REALTYPE,SIMD>& w, const int ex_w)
{
  return v.dotc(ex, w, ex_w);  // <v[ex] | w([ex_w] >
}


//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD>& v, const REALTYPE a, const AField<REALTYPE,SIMD> &w)
{
  // |v> += a |w>
#ifdef DEBUG
  vout.paranoiac("aypx with 3 arguments: v.m_ptr=%p, w.m_ptr=%p\n", v.ptr_vcomplex(), w.ptr_vcomplex());
#endif
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void aypx(const REALTYPE a, AField<REALTYPE,SIMD>& v,  const AField<REALTYPE,SIMD> &w)
{
  // |v> := a |v> + |w>
  v.aypx(a, w);
}


//====================================================================
template <typename REALTYPE>
void aypx(const typename AField<REALTYPE,SIMD>::complex_t  a, AField<REALTYPE,SIMD> & v, const AField<REALTYPE,SIMD> &w)
{
  // |v> := a |v> + |w>
  v.aypx(a, w);
}


//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD>& v, const int ex, const REALTYPE a,
          const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  // |v[ex]> += a|w[ex_w]>
  v.axpy(ex, a, w, ex_w);
}


//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD>& v, const int ex, const typename AField<REALTYPE,SIMD>::complex_t a,
          const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  // |v[ex]> += a|w[ex_w]>
  v.axpy(ex, a, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD>& v, const int ex, const typename AField<REALTYPE,SIMD>::real_t a,
          const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  // |v[ex]> += a|w[ex_w]>
  v.axpy(ex, a, w, ex_w);
}

//====================================================================
template <typename REALTYPE>
void axpy(AField<REALTYPE,SIMD>& v, const typename AField<REALTYPE,SIMD>::complex_t a,
          const AField<REALTYPE,SIMD> &w)
{
  // |v> += a|w>
  v.axpy(a, w);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE,SIMD>& v, const REALTYPE a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
void scal(AField<REALTYPE,SIMD>& v, const typename AField<REALTYPE,SIMD>::complex_t a)
{
  v.scal(a);
}

//====================================================================
template <typename REALTYPE>
void add(AField<REALTYPE,SIMD>& v, const int ex,
         const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  // |v[ex]> += |w[ex_w]>
  v.add(ex, w, ex_w);
}
//====================================================================
template <typename REALTYPE>
void add(AField<REALTYPE,SIMD>& v, const AField<REALTYPE,SIMD> &w)
{
  // |v> += |w>
  assert(v.nex()==w.nex());
  for(int ex=0; ex<v.nex(); ++ex){
    v.add(ex, w, ex);
  }
}
#define AFIELD_HAS_ADD


//====================================================================
template <typename REALTYPE>
void sub(AField<REALTYPE,SIMD>& v, const int ex,
         const AField<REALTYPE,SIMD> &w, const int ex_w)
{
  // |v[ex]> -= |w[ex_w]>
  v.sub(ex, w, ex_w);
}
//====================================================================
template <typename REALTYPE>
void sub(AField<REALTYPE,SIMD>& v, const AField<REALTYPE,SIMD> &w)
{
  // |v> -= |w>
  assert(v.nex()==w.nex());
  for(int ex=0; ex<v.nex(); ++ex){
    v.sub(ex, w, ex);
  }
}
#define AFIELD_HAS_SUB


//====================================================================
template<typename REALTYPE>
void dotc_and_norm2(typename AField<REALTYPE,SIMD>::complex_t &dotc,
                    REALTYPE &v_norm2, REALTYPE &w_norm2,
                    const AField<REALTYPE,SIMD>& v,
                    const AField<REALTYPE,SIMD>& w)
{
  // returns <v|w>, <v|v>, <w|w>
  assert(v.nex()==w.nex());
  assert(v.nvol()==w.nvol());
  assert(v.nin()==w.nin());

  dotc=v.dotc_and_norm2(v_norm2, w_norm2, w);
}
#define AFIELD_HAS_DOTC_AND_NORM2

//====================================================================
template <typename AFIELD>
void gramschmidt(AFIELD& v)
{
  assert(v.field_element_type()==Element_type::COMPLEX);
  //typedef AField<REALTYPE,SIMD> AFIELD;
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vreal_t vreal_t;
  typedef typename AFIELD::complex_t complex_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;

  // size of AFeild
  const int Nin  = v.nin_complex();
  const int Rvol = v.rvol();
  const int Nex  = v.nex();
  const int size = Rvol*Nin;
  int ith, nth, is, ns;

  set_threadtask(ith, nth, is, ns, Rvol);
  is*=Nin;
  ns*=Nin;
  vcomplex_t *pv=v.ptr_vcomplex();

#pragma omp barrier

  for(int i=0; i<Nex; ++i){
    vcomplex_t *pv_i=pv+i*size;
    for(int j=0; j<i; ++j){
      // csum= <j|i>
      vcomplex_t csum;
      set_zero(csum);
      const vcomplex_t *pv_j=pv+j*size;
      for(int s=is; s<ns; ++s){
        vcomplex_t vj=pv_j[s];
        vcomplex_t vi=pv_i[s];
        accum_cmult_conj(csum,vj,vi); // sum+=conj(vj)*vi
      }
      complex_t scalar_csum=reduce(csum);
      real_t re, im;
      re=real(scalar_csum);
      im=imag(scalar_csum);
      real_t tmp[2]={re,im};
      ThreadManager_OpenMP::reduce_sum_global(tmp, 2, ith, nth);
      re=tmp[0];
      im=tmp[1];

        // |i> := |i> - sum_j <j|i> |j>
      vcomplex_t c;
      bcast(c, complex_t(-re, -im));
      for(int s=is; s<ns; ++s){
        vcomplex_t vj=pv_j[s];
        vcomplex_t vi=pv_i[s];
        accum_cmult(vi, c, vj);
        pv_i[s]=vi;
      }
    } // j
    // normalize |i>
    vreal_t sum2;
    set_zero(sum2);
    for(int s=is; s<ns; ++s){
      const vcomplex_t vi=pv_i[s];
      accum_abs2(sum2, vi);
      pv_i[s]=vi;
    }
    real_t scalar_sum2=reduce(sum2);
    ThreadManager_OpenMP::reduce_sum_global(scalar_sum2, ith, nth);
    vreal_t factor;
    bcast(factor, 1.0/sqrt(scalar_sum2));
    for(int s=is; s<ns; ++s){
      scal(pv_i[s], factor);
    }
  } // i
#pragma omp barrier
}





template <class INDEX, class AFIELD>
void convert_strict(INDEX& index, AFIELD& v, const Field& w)
{
  // w = v
  // requires strict matching of the shape
  const size_t Nin2 = w.nin()/2;
  const size_t Nvol = w.nvol();
  const size_t Nex  = w.nex(); // assumes w is COMPLEX
  assert(w.field_element_type()==Element_type::COMPLEX);
  assert(v.field_element_type()==Element_type::COMPLEX);

  // size of AFeild
  const size_t ANin  = v.nin_complex();
  const size_t ANvol = v.nvol();
  const size_t ANex  = v.nex();

  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;


#ifdef DEBUG
  vout.paranoiac("convert:  Nin,  Nvol,  Nex= %d, %d, %d\n",Nin2*2, Nvol, Nex);
  vout.paranoiac("convert: ANin, ANvol, ANex= %d, %d, %d\n",ANin, ANvol, ANex);
#endif

  assert(v.check_size(2*Nin2, Nvol, Nex));


#pragma omp barrier
  int is, ns,  ith, nth;

  set_threadtask(ith, nth, is, ns, v.rvol());
  const int simd_size=vcomplex_t::Nsimd();
  std::vector<size_t> site_array(simd_size);
  aligned_vector<vcomplex_t> buf;
  buf.resize(Nin2);
  for(int ex=0; ex<Nex; ++ex){
    vcomplex_t *pv=v.ptr_vcomplex()+ex*Nin2*v.rvol();

    for(int rs=is; rs<ns; ++rs){
      index.vsite_to_site(rs,site_array);
      real_t *pbuf= (real_t*)(&buf[0]);
      for(int i=0; i<simd_size; ++i){
        size_t site=site_array[i];
        for(int ain = 0; ain < ANin; ++ain){
          size_t iw=ain+Nin2*site;
          int iv= ain*simd_size+i;
          pbuf[2*iv]=w.cmp(2*iw);
          pbuf[2*iv+1]=w.cmp(2*iw+1);
        }
      }
      for(int ain=0; ain<ANin; ++ain){
        sstore(pv+ANin*rs+ain, buf[ain]);
      }
    }
  }
#pragma omp barrier

}

//====================================================================
template <class INDEX, class AFIELD>
void revert_strict(INDEX& index, Field& v, const AFIELD& w)
{
  // v = w
  // requires strict matching of the shape
  size_t Nin2  = v.nin()/2;
  size_t Nvol = v.nvol();
  size_t Nex  = v.nex();

  assert(w.field_element_type()==Element_type::COMPLEX);
  assert(v.field_element_type()==Element_type::COMPLEX);
  w.check_size(2*Nin2, Nvol, Nex);  // assumes v is COMPLEX

  const int Nx = index.Nx();
  const int Ny = index.Ny();
  const int Nz = index.Nz();
  const int Nt = index.Nt();

#pragma omp barrier
#pragma omp master
  {
    for(int ex=0; ex<Nex; ++ex){
      for(int t=0; t< Nt; ++t){
        for(int z=0; z< Nz; ++z){
          for(int y=0; y< Ny; ++y){
            for(int x=0; x< Nx; ++x){
              size_t site=Nx * (Ny * (Nz * t + z) + y) + x;
              for(int in = 0; in < Nin2; ++in){
                size_t iv = in + Nin2 * (site + Nvol * ex);
                size_t iw=index.scalar_idx(in, Nin2, x,y,z,t, ex);
                v.set(2*iv, double( w.e(2*iw)));
                v.set(2*iv+1, double( w.e(2*iw+1)));
              }
            }}}}
    }
  }
#pragma omp barrier
}

template <class INDEX, class AFIELD>
void convert(INDEX& index, AFIELD& v, const Field& w) 
{
  // w = v
  // allows conversion from inner/outer to outer/inner:
  //    in + Nin*ex is kept

  // size of Feild
  const size_t Nin2 = w.nin()/2;
  const size_t Nvol = w.nvol();
  const size_t Nex  = w.nex();

  // size of AFeild
  const size_t ANin  = v.nin_complex();
  const size_t ANvol = v.nvol();
  const size_t ANex  = v.nex();

  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;

  vout.paranoiac("afield_simd, convert: Nin2,  Nvol,  Nex= %d, %d, %d\n",Nin2, Nvol, Nex);
  vout.paranoiac("afield_simd, convert: ANin, ANvol, ANex= %d, %d, %d\n",ANin,ANvol,ANex);

  assert(w.field_element_type()==Element_type::COMPLEX);
  assert(v.field_element_type()==Element_type::COMPLEX);
  assert( Nvol==ANvol && Nin2*Nex==ANin*ANex);

  const int simd_size=vcomplex_t::Nsimd();


#pragma omp barrier
  int is, ns,  ith, nth;

  set_threadtask(ith, nth, is, ns, v.rvol());
  std::vector<size_t> site_array(simd_size);
  aligned_vector<vcomplex_t> buf(ANin);
  for(int aex=0; aex<ANex; ++aex){
    vcomplex_t *pv=v.ptr_vcomplex()+aex*ANin*v.rvol();

    for(int rs=is; rs<ns; ++rs){
      index.vsite_to_site(rs,site_array);  // reduced index to index for Field
      real_t *pbuf= (real_t*)(&buf[0]);
      for(int i=0; i<simd_size; ++i){      // run over inside the simd variable
        size_t site=site_array[i];
        for(int ain = 0; ain < ANin; ++ain){
          size_t id =ain + ANin*aex;
          size_t ex = id/Nin2;
          size_t in = id%Nin2;
          size_t iw = in + Nin2 * (site + Nvol * ex);

          int iv= ain*simd_size+i;
          pbuf[2*iv]=w.cmp(2*iw);
          pbuf[2*iv+1]=w.cmp(2*iw+1);
        }
      }
      for(int ain=0; ain< ANin; ++ain){
        sstore(pv+ANin*rs+ain, buf[ain]);
      }
    }
  }

  vout.paranoiac("convert, done\n");
#pragma omp barrier
}


//====================================================================
template <class INDEX, class AFIELD>
void revert(INDEX& index, Field& v, const AFIELD& w) 
{
  // v = w
  // allows conversion from inner/outer to outer/inner:
  //    in + Nin*ex is kept

  // size of Feild
  const size_t Nin2 = v.nin()/2;
  const size_t Nvol = v.nvol();
  const size_t Nex  = v.nex();

  // size of AFeild
  const size_t ANin  = w.nin_complex();
  const size_t ANvol = w.nvol();
  const size_t ANex  = w.nex();

  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;

  vout.paranoiac("afield_simd, convert: Nin2,  Nvol,  Nex= %d, %d, %d\n",Nin2, Nvol, Nex);
  vout.paranoiac("afield_simd, convert: ANin, ANvol, ANex= %d, %d, %d\n",ANin,ANvol,ANex);

  assert(w.field_element_type()==Element_type::COMPLEX);
  assert(v.field_element_type()==Element_type::COMPLEX);
  assert( Nvol==ANvol && Nin2*Nex==ANin*ANex);

  const int simd_size=vcomplex_t::Nsimd();

#pragma omp barrier
  int is, ns,  ith, nth;

  set_threadtask(ith, nth, is, ns, w.rvol());
  std::vector<size_t> site_array(simd_size);
  aligned_vector<vcomplex_t> buf(ANin);
  for(int aex=0; aex<ANex; ++aex){
    vcomplex_t *pw=w.ptr_vcomplex()+aex*ANin*w.rvol();
    for(int rs=is; rs<ns; ++rs){
      index.vsite_to_site(rs,site_array);  // reduced index to index for Field
      real_t *pbuf= (real_t*)(&buf[0]);
      for(int ain=0; ain< ANin; ++ain){
        buf[ain]=pw[ANin*rs+ain];     // read data in unit of simd variable
      }
      for(int i=0; i<simd_size; i++){ // run over inside the simd variable
        size_t site=site_array[i];
        for(int ain = 0; ain < ANin; ++ain){
          size_t id =ain + ANin*aex;
          size_t ex = id/Nin2;
          size_t in = id%Nin2;
          size_t iv = in + Nin2 * (site + Nvol * ex);
          int iw= ain*simd_size+i;
          v.set(2*iv,pbuf[2*iw]);
          v.set(2*iv+1,pbuf[2*iw+1]);
        }
      }

    }
  }

  vout.paranoiac("revert, done\n");
#pragma omp barrier

}

//====================================================================
// alias
template <class INDEX, class AFIELD>
void reverse(INDEX& index, Field& v, const AFIELD& w)
{
  revert(index, v, w);
}


//====================================================================
template <class INDEX, class FIELD>
void convert_spinor(INDEX& index, FIELD& v, const Field& w) 
{
  int Nin  = w.nin();
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2*Nc*Nd);
  convert(index, v, w);
}


//====================================================================
template <class INDEX2, class FIELD2, class INDEX1, class FIELD1>
void convert(INDEX2& index2, FIELD2& w,
             const INDEX1& index1, const FIELD1& v)
{
  convert_h(index2, w, index1, v);
}



//====================================================================
template <class INDEX1, class FIELD1, class INDEX2, class FIELD2>
void convert_h(INDEX1& index1, FIELD1& v, INDEX2& index2, const FIELD2& w) {

  // size of AField 2
  const size_t ANin2 = w.nin_complex();
  const size_t ANvol2 = w.nvol();
  const size_t ANex2  = w.nex();

  // size of AField 1
  const size_t ANin1  = v.nin_complex();
  const size_t ANvol1 = v.nvol();
  const size_t ANex1  = v.nex();

  typedef typename FIELD1::real_t real1_t;
  typedef typename FIELD1::complex_t complex1_t;
  typedef typename FIELD1::vcomplex_t vcomplex1_t;
  typedef typename FIELD2::real_t real2_t;
  typedef typename FIELD2::complex_t complex2_t;
  typedef typename FIELD2::vcomplex_t vcomplex2_t;

  vout.paranoiac("convert from size %d to %d.\n", sizeof(real2_t), sizeof(real1_t));

  assert(w.field_element_type()==Element_type::COMPLEX);
  assert(v.field_element_type()==Element_type::COMPLEX);
  assert( ANvol1==ANvol2 && ANin1*ANex1==ANin2*ANex2);

  const int simd_size=vcomplex1_t::Nsimd();  // use the simd size of the output
  const complex2_t *pw=w.ptr_complex();
#pragma omp barrier

  int is, ns,  ith, nth;
  set_threadtask(ith, nth, is, ns, v.rvol());
  std::vector<size_t> site_array(simd_size);
  aligned_vector<vcomplex1_t> buf(ANin1);
  for(int aex=0; aex<ANex1; ++aex){
    vcomplex1_t *pv=v.ptr_vcomplex()+aex*ANin1*v.rvol();
    for(int rs=is; rs<ns; ++rs){
      index1.vsite_to_site(rs,site_array);  // reduced index to index for Field
      int iv=0;
      complex1_t *pbuf=(complex1_t*)(&buf[0]);
      for(int ain = 0; ain < ANin1; ++ain){
        size_t id =ain + ANin1*aex;
        size_t ex2 = id/ANin2;
        size_t in2 = id%ANin2;
        for(int i=0; i<simd_size; ++i){      // run over inside the simd variable
          size_t site=site_array[i];
          size_t iw = index2.scalar_idx(in2, ANin2, site, ex2);
          pbuf[iv]=pw[iw];
          iv++;
        }
      }
      for(int ain=0; ain< ANin1; ++ain){
        sstore(pv+ANin1*rs+ain, buf[ain]);
      }
    }
  }

  vout.paranoiac("convert, done\n");
#pragma omp barrier

}

  //====================================================================
template <class INDEX, class FIELD>
void revert_spinor(INDEX& index, Field& v, FIELD& w){  // renamed
  int Nin  = v.nin();
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2*Nc*Nd);
  revert(index, v, w);
}

//====================================================================
template <class INDEX, class FIELD>
void reverse_spinor(INDEX& index, Field& v, FIELD& w){  // renamed
  revert_spinor(index,v,w);
}


template <class INDEX, class AFIELD>
void convert_gauge(INDEX& index, AFIELD& v, const Field& w) {
  convert(index, v, w);
}

template <class INDEX, class FIELD>
void reverse_gauge(INDEX& index, Field& v, FIELD& w){
  revert(index,v,w);
}

#endif
