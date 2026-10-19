/*!

        @file    $Id: afield.cpp #$

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
*/


#include "afield.h"
#include <cassert>

#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_SIMD/inline/vector_arith-inc.h"


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
}



template<typename REALTYPE>
const std::string AField<REALTYPE,SIMD>::class_name = "SIMD::AField";

//template<> Timer AField<double, SIMD>::timer_dotc_norm2("dotc_and_norm2");
//template<> Timer AField<double, SIMD>::timer_dotc_norm2_calc("dotc_and_norm2_calc");
//template<> Timer AField<float, SIMD>::timer_dotc_norm2("dotc_and_norm2");
//template<> Timer AField<float, SIMD>::timer_dotc_norm2_calc("dotc_and_norm2_calc");


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::init(const int nin, const int nvol, const int nex)
{

  vout.paranoiac("AField init, started:   nin=%4d,   nvol=%4d,   nex=%4d\n", nin,nvol,nex);
  vout.paranoiac("AField init,          m_nin=%4d, m_nvol=%4d, m_nex=%4d, m_ptr=%p,\n",m_nin,m_nvol,m_nex, m_ptr);

  if(m_ptr){
    vout.paranoiac("AField init: deleting m_ptr=%p\n",m_ptr);
    delete(m_ptr);
    vout.paranoiac("init: m_ptr has been deleted\n");
  }

  m_nin  = nin;
  m_nvol = nvol;
  m_nex  = nex;
  m_nsize = m_nin * m_nvol * m_nex;

  m_rsize=0;
  //  m_field=0;
  if(m_element_type == Element_type::REAL){
    using Field_t=Field_vreal;
    m_ncomplex=1;
    m_nin_complex=nin;
    m_nsize_complex=m_nsize;
    m_rvol = m_nvol/( vreal_t::Nsimd());
    if(nin*nex>0){
      if(m_rvol * vreal_t::Nsimd() != m_nvol){
        vout.crucial("%s: bad nvol (too small?), must be a multiple of the simd length %d:  (given nvol=%d)\n",
                     class_name.c_str(), vreal_t::Nsimd(), m_nvol);
        exit(EXIT_FAILURE);
      }
      m_rsize=m_nsize/( vreal_t::Nsimd());
      Field_t *ptr_field=new Field_t(m_rsize);
      m_ptr=ptr_field;
    }
  } else if(m_element_type == Element_type::COMPLEX){
    using Field_t=Field_vcomplex;
    m_ncomplex=2;
    m_nin_complex=nin/2;
    m_nsize_complex=m_nsize/2;
    m_rvol = m_nvol/( vcomplex_t::Nsimd() );
    if(nin % 2 ==1){
      vout.crucial("AField init: new as a complex but nin is odd: nin=%d, nex=%d\n", nin, nex);
    }
    if(nin*nex>0){
      if(m_rvol * vcomplex_t::Nsimd() != m_nvol){
        vout.crucial("%s: bad nvol (too small?), must be a multiple of the simd length %d:  (given nvol=%d)\n",
                     class_name.c_str(), vcomplex_t::Nsimd(), m_nvol);
        exit(EXIT_FAILURE);
      }
      m_rsize=m_nsize/( vcomplex_t::Nsimd() );
      Field_t *ptr_field=new Field_t(m_rsize);
      m_ptr=ptr_field;
    }
  } else {
    vout.crucial("%s: cannot happen: improper element_type: %d\n", class_name.c_str(),m_element_type);
    vout.crucial("%s:   should be REAL (=%d)\n", class_name.c_str(), Element_type::REAL);
    vout.crucial("%s:   or     COMPLEX (=%d)\n", class_name.c_str(), Element_type::COMPLEX);
    exit(EXIT_FAILURE);
  }
  if(m_nsize>0){
#ifdef AFIELD_FIRST_TOUCH
#pragma omp parallel
    {
      set(0.0);
    }
#endif
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::tidyup()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);
  if(m_ptr){
    delete m_ptr;  // calls destrctor of Field_body through virtual destructor
    m_ptr=0;
  }
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::set(const REALTYPE a)
{
  if(m_element_type == Element_type::REAL){
    vreal_t *out=ptr_vreal();
    vreal_t va;
    bcast(va,a);
    int ith, nth, is, ns;

    set_threadtask(ith, nth, is, ns, m_rvol);
    is*=m_nin_complex;
    ns*=m_nin_complex;
    size_t size=m_rvol*m_nin_complex;
    for(int ex=0; ex<m_nex; ++ex){
      for(int s=is; s<ns; ++s){
        // sstore(out+s,va);
        out[s]=va;
      }
      out+=size;
    }
  }else if(m_element_type == Element_type::COMPLEX){
    vcomplex_t *out=ptr_vcomplex();
    vcomplex_t vc;
    complex_t c(a,a);  // not very intuitive
    bcast(vc,c);
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_rvol);
    is*=m_nin_complex;
    ns*=m_nin_complex;
    size_t size=m_rvol*m_nin_complex;
    for(int ex=0; ex<m_nex; ++ex){
      for(int s=is; s<ns; ++s){
        //      sstore(out+s,vc);
        out[s]=vc;
      }
      out+=size;
    }
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::setc(const typename AField<REALTYPE,SIMD>::complex_t a)
{
  if(m_element_type == Element_type::REAL){
    vout.crucial("%s: complex parameter for real field %s", class_name.c_str(), __func__);
  }else if(m_element_type == Element_type::COMPLEX){

    vout.paranoiac("AField::set(complex a): a=%f %f, m_nsize_complex=%d\n",a.real(), a.imag(), m_nsize_complex);
    vcomplex_t *out=ptr_vcomplex();
    vcomplex_t va;
    bcast(va,a);
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_rvol);
    is*=m_nin_complex;
    ns*=m_nin_complex;
    size_t size=m_rvol*m_nin_complex;
    for(int ex=0; ex<m_nex; ++ex){
      for(int s=is; s<ns; ++s){
        // sstore(out+s,va);
        out[s]=va;
      }
      out+=size;
    }
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::set(const typename AField<REALTYPE,SIMD>::complex_t a)
{
  setc(a);
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::setc(const REALTYPE a)
{
  complex_t c(a,0);
  setc(a);
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::copy(const AField<REALTYPE,SIMD>& w)
{
  assert(check_size(w));
  assert(m_element_type == w.field_element_type());

  vreal_t *out=ptr_vreal();
  const vreal_t *in=w.ptr_vreal();
  if(out == in){
    return;
  }

  int ith, nth, is, ns;

  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  size_t size=m_rvol*m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      //    sstore(out+s,in[s]);
      out[s]=in[s];
    }
    out+=size;
    in+=size;
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::copy(const int index,
                                 const AField<REALTYPE,SIMD>& w, const int index_w, const int size)
{

  vout.crucial("%s: partial copy is not supported", class_name.c_str());
  exit(EXIT_FAILURE);

  /*
  assert(index + size <= m_nsize);
  assert(index_w + size <= w.size());

  assert(check_size(w));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size);

  for(int i = is; i < ns; ++i){
    m_field[index + i] = w.m_field[index_w + i];
  }
  */

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::copy(const int ex,
                                 const AField<REALTYPE,SIMD>& w, const int ex_w){
  const size_t offset_in=m_nin_complex*m_rvol*ex_w;
  const size_t offset_out=m_nin_complex*m_rvol*ex;

  vreal_t *out=ptr_vreal() + offset_out;
  const vreal_t *in =w.ptr_vreal() + offset_in;
  if(out == in){
    return;
  }

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int s=is; s<ns; ++s){
    //    sstore(out+s,in[s]);
    out[s]=in[s];
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::axpy(const REALTYPE a, const AField<REALTYPE,SIMD>& w)
{
  // |this> += a|w>

#ifdef DEBUG
  vout.paranoiac("AField::axpy (2 arg version) is called\n");
  vout.paranoiac("AField::axpy: m_ptr=%p, w.m_ptr=%p\n", m_ptr, w.m_ptr);
#endif
  assert(check_size(w));

  vreal_t *inout=ptr_vreal();
  const vreal_t *in=w.ptr_vreal();
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  size_t size=m_rvol*m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      const vreal_t x=in[s];
      vreal_t y=inout[s];
      vaxpy(y,va,x);
      // sstore(inout+s, y);
      inout[s]=y;
    }
    in+=size;
    inout+=size;
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::aypx(const REALTYPE a, const AField<REALTYPE,SIMD>& w)
{
  // |this> := a |this> + |w>

#ifdef DEBUG
  vout.paranoiac("AField::aypx (2 arg version) is called\n");
  vout.paranoiac("AField::aypx: m_ptr=%p, w.m_ptr=%p\n", m_ptr, w.m_ptr);
#endif
  assert(check_size(w));

  vreal_t *inout=ptr_vreal();
  const vreal_t *in=w.ptr_vreal();
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_rvol*m_nin_complex;
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      const vreal_t x=in[s];
      vreal_t y=inout[s];
      vaypx(y,va,x);
      // sstore(inout+s,y);
      inout[s]=y;
    }
    in+=size;
    inout+=size;
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::aypx(const complex_t a, const AField<REALTYPE,SIMD>& w)
{
  // |this> := a |this> + |w>

  assert(check_size(w));

  vcomplex_t *inout=ptr_vcomplex();
  const vcomplex_t *in=w.ptr_vcomplex();
  vcomplex_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_rvol*m_nin_complex;
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      const vcomplex_t x=in[s];
      vcomplex_t y=inout[s];
      vaypx(y,va,x);  // y= a*y+x
      // sstore(inout+s,y);
      inout[s]=y;
    }
    in+=size;
    inout+=size;
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::axpy(const int index, const REALTYPE a,
                                 const AField<REALTYPE,SIMD>& w, const int index_w,
                                 const int size)
{

  vout.crucial("%s: partial axpy is not supported", class_name.c_str());
  exit(EXIT_FAILURE);

  /*
  assert(index + size <= m_nsize);
  assert(index_w + size <= w.size());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size);

  for (int i = is; i < ns; ++i) {
    m_field[index + i] += a * w.m_field[index_w + i];
  }
  */

}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::axpy(const int ex, const REALTYPE a,
                                 const AField<REALTYPE,SIMD>& w, const int ex_w)
{
  // |this[ex] >  += a |w[w_ex] >

#ifdef DEBUG
  vout.paranoiac("AField::axpy (4 arg version) is called: a=%f\n",a);
#endif

  // TODO: prefetch

  assert(ex<=m_nex);
  assert(ex_w<=w.nex());
  const size_t offset_in=m_nin_complex*m_rvol*ex_w;
  const size_t offset_out=m_nin_complex*m_rvol*ex;

  vreal_t *out=ptr_vreal() + offset_out;
  const vreal_t *in=w.ptr_vreal()+ offset_in;
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int s=is; s<ns; ++s){
    vreal_t y=out[s];
    vreal_t x=in[s];
    vaxpy(y, va, x); // y := va*x + y
    out[s]=y;
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::axpy(const int ex, const complex_t a,
                                 const AField<REALTYPE,SIMD>& w, const int ex_w)
{
  // |this[ex] >  += a |w[w_ex] >

#ifdef DEBUG
  vout.paranoiac("AField::axpy (4 arg version) is called: a=%f\n",a);
#endif

  // TODO: prefetch

  assert(ex<=m_nex);
  assert(ex_w<=w.nex());
  const size_t offset_in=m_nin_complex*m_rvol*ex_w;
  const size_t offset_out=m_nin_complex*m_rvol*ex;

  vcomplex_t *out=ptr_vcomplex() + offset_out;
  const vcomplex_t *in=w.ptr_vcomplex()+ offset_in;
  vcomplex_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int s=is; s<ns; ++s){
    vcomplex_t y=out[s];
    vcomplex_t x=in[s];
    accum_cmult(y,va,x); // y := va*x + y
    out[s]=y;
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::axpy(const complex_t a,
                                 const AField<REALTYPE,SIMD>& w)
{
  // |this >  += a |w>

#ifdef DEBUG
  vout.paranoiac("AField::axpy (2 arg version) is called: a=%f\n",a);
#endif

  // TODO: prefetch

  assert(m_nex==w.nex());
  assert(m_nin==w.nin());

  vcomplex_t *out_base=ptr_vcomplex();
  const vcomplex_t *in_base=w.ptr_vcomplex();
  vcomplex_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    vcomplex_t *out=out_base+ex*m_rvol*m_nin_complex;
    const vcomplex_t *in=in_base+ex*m_rvol*m_nin_complex;
    for(int s=is; s<ns; ++s){
      vcomplex_t y=out[s];
      vcomplex_t x=in[s];
      accum_cmult(y,va,x); // y := va*x + y
      out[s]=y;
    }
  }
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::scal(const REALTYPE a)
{
  vreal_t *inout=ptr_vreal();
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_rvol*m_nin_complex;
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      vreal_t vv=inout[s];
      //    vstream(inout[s],vv*va);
      inout[s]=vv*va;
    }
    inout+=size;
  }
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::scal(const complex_t a)
{
  vcomplex_t *inout=ptr_vcomplex();
  vcomplex_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_rvol*m_nin_complex;
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    for(int s=is; s<ns; ++s){
      vcomplex_t vv=inout[s];
      //    vstream(inout[s],vv*va);
      inout[s]=vv*va;
    }
    inout+=size;
  }
}

//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,SIMD>::dot(const AField<REALTYPE,SIMD>& w) const
{
  //  dot product as real vectors
  assert(check_size(w));

  // to avid loss of trailing digits (espeically for real_t=float),
  // uses partial summations

  if(true){
    const vreal_t *ptr0=ptr_vreal();
    const vreal_t *w_ptr0=w.ptr_vreal();
    vreal_t sum_all;
    set_zero(sum_all);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_rvol);
    size_t size=m_nin_complex;
    for(int ex=0; ex<m_nex; ++ex){
      const vreal_t *ptr  =ptr0   + is*size + ex*m_rvol*m_nin_complex;
      const vreal_t *w_ptr=w_ptr0 + is*size + ex*m_rvol*m_nin_complex;
      vreal_t sum;
      set_zero(sum);
      for(int s=is; s<ns; ++s){
        vreal_t sum_tmp;
        set_zero(sum_tmp);
        for(int a=0; a<size; ++a){
          vreal_t vv=ptr[a];
          vreal_t ww=w_ptr[a];
          vaxpy(sum_tmp,vv,ww); // sum+=vv*ww+sum
        }
        accum_add(sum, sum_tmp);
        ptr+=size;
        w_ptr+=size;
      }
      accum_add(sum_all, sum);
    }
    double scalar_sum=reduce(sum_all);
    ThreadManager_OpenMP::reduce_sum_global(scalar_sum, ith, nth);
    real_t result=scalar_sum;
    return result;

  } else {
    vout.crucial("%s: dot for complex field is called\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

//====================================================================
template<typename REALTYPE>
typename AField<REALTYPE,SIMD>::complex_t AField<REALTYPE,SIMD>::dotc(const AField<REALTYPE,SIMD>& w) const
{
  // returns <this|w>

  //  dot product as complex vectors
  assert(check_size(w));

  if(m_element_type == Element_type::REAL){
    vout.crucial("%s: dotc for real field is called\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // to avid loss of trailing digits (espeically for real_t=float),
  // uses partial summations

  const vcomplex_t *ptr0=ptr_vcomplex();
  const vcomplex_t *w_ptr0=w.ptr_vcomplex();
  vcomplex_t sum_all;
  set_zero(sum_all);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    vcomplex_t sum;
    set_zero(sum);
    const vcomplex_t *ptr  =ptr0   + is*size + ex*m_rvol*m_nin_complex;
    const vcomplex_t *w_ptr=w_ptr0 + is*size + ex*m_rvol*m_nin_complex;
    for(int s=is; s<ns; ++s){
      vcomplex_t sum_tmp;
      set_zero(sum_tmp);
      for(int a=0; a<size; ++a){
        vcomplex_t vv=ptr[a];
        vcomplex_t ww=w_ptr[a];
        accum_cmult_conj(sum_tmp,vv,ww); // sum_tmp+=conj(vv)*ww
      }
      accum_add(sum, sum_tmp);
      ptr+=size;
      w_ptr+=size;
    }
    accum_add(sum_all, sum);
  }
  dcomplex scalar_sum=reduce(sum_all);
  double re, im;
  re=real(scalar_sum);
  im=imag(scalar_sum);
  double sum_array[2]={re, im};
  ThreadManager_OpenMP::reduce_sum_global(sum_array, 2, ith, nth);
  return complex_t(sum_array[0], sum_array[1]);
}




//====================================================================
template<typename REALTYPE>
typename AField<REALTYPE,SIMD>::complex_t AField<REALTYPE,SIMD>::dotc(const int ex, const AField<REALTYPE,SIMD>& w, const int ex_w) const
{
  // returns <this[ex]|w[ex_w]>

  if(m_element_type == Element_type::REAL){
    vout.crucial("%s: dotc for real field is called\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // to avid loss of trailing digits (espeically for real_t=float),
  // uses partial summations

  const vcomplex_t *ptr0  =ptr_vcomplex()+ex*m_rvol*m_nin_complex;
  const vcomplex_t *w_ptr0=w.ptr_vcomplex()+ex_w*m_rvol*m_nin_complex;;
  vcomplex_t sum;
  set_zero(sum);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_nin_complex;
  const vcomplex_t *ptr  =ptr0  +m_nin_complex*is;
  const vcomplex_t *w_ptr=w_ptr0+m_nin_complex*is;
  for(int s=is; s<ns; ++s){
    vcomplex_t sum_tmp;
    set_zero(sum_tmp);
    for(int a=0; a<size; ++a){
      vcomplex_t vv=ptr[a];
      vcomplex_t ww=w_ptr[a];
      accum_cmult_conj(sum_tmp,vv,ww); // sum+=conj(vv)*ww
    }
    ptr+=size;
    w_ptr+=size;
    accum_add(sum, sum_tmp);
  }
  dcomplex scalar_sum=reduce(sum);
  double re, im;
  re=real(scalar_sum);
  im=imag(scalar_sum);
  double sum_array[2]={re, im};
  ThreadManager_OpenMP::reduce_sum_global(sum_array, 2, ith, nth);
  return complex_t(sum_array[0], sum_array[1]);
}


//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE,SIMD>::norm2() const
{
  // to avid loss of trailing digits (espeically for real_t=float),
  // uses partial summations
  const vreal_t *in0=ptr_vreal();
  vreal_t sum2_all;
  set_zero(sum2_all);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_nin_complex;
  for(int ex=0; ex<m_nex; ++ex){
    vreal_t sum2;
    set_zero(sum2);
    for(int s=is; s<ns; ++s){
      vreal_t sum2_tmp;
      set_zero(sum2_tmp);
      const vreal_t *in=in0+ s*m_nin_complex + ex*m_nin_complex*m_rvol;
      for(int a=0; a<m_nin_complex; ++a){
        vreal_t vv=in[a];
        vaxpy(sum2_tmp,vv,vv);
      }
      accum_add(sum2, sum2_tmp);
    }
    accum_add(sum2_all, sum2);
  }
  double scalar_sum2=reduce(sum2_all);

  ThreadManager_OpenMP::reduce_sum_global(scalar_sum2, ith, nth);
  real_t result=scalar_sum2;
  return result;
}

//====================================================================
template<typename REALTYPE>
typename AField<REALTYPE,SIMD>::complex_t AField<REALTYPE,SIMD>::dotc_and_norm2(
                                           REALTYPE &norm2, REALTYPE &w_norm2,
                                           const AField<REALTYPE,SIMD>& w) const
{
  if(m_element_type == Element_type::REAL){
    vout.crucial("%s: dotc_and_norm2 for real field is called\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  assert(check_size(w));
  //#pragma omp master
  //  {
  //    timer_dotc_norm2.start();
  //  }

  // to avid loss of trailing digits (espeically for real_t=float),
  // uses partial summations

  const vcomplex_t *ptr0=ptr_vcomplex();
  const vcomplex_t *w_ptr0=w.ptr_vcomplex();
  dcomplex dotc_sum=0.0;
  double norm2_sum=0.0;
  double w2_sum=0.0;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  size_t size=m_nin_complex;

  //#pragma omp master
  //  {
  //    timer_dotc_norm2_calc.start();
  //  }

  // loop fission version is slightly faster on KNL
  for(int ex=0; ex<m_nex; ++ex){
    const vcomplex_t *ptr  =ptr0   + is*size + ex*m_rvol*m_nin_complex;
    const vcomplex_t *w_ptr=w_ptr0 + is*size + ex*m_rvol*m_nin_complex;
    vcomplex_t dotc_vsum;
    vreal_t norm2_vsum;
    vreal_t w2_vsum;
    set_zero(dotc_vsum);
    set_zero(norm2_vsum);
    set_zero(w2_vsum);

    for(int s=is; s<ns; ++s){
      vcomplex_t tmp_dotc;
      set_zero(tmp_dotc);
      for(int a=0; a<size; ++a){
        vcomplex_t vv=ptr[a];
        vcomplex_t ww=w_ptr[a];
        accum_cmult_conj(tmp_dotc,vv,ww); // tmp_dotc+=conj(vv)*ww
      }
      accum_add(dotc_vsum, tmp_dotc);

      vreal_t tmp_norm2;
      set_zero(tmp_norm2);
      for(int a=0; a<size; ++a){
        vcomplex_t vv=ptr[a];
        accum_abs2(tmp_norm2, vv);
      }
      accum_add(norm2_vsum, tmp_norm2);

      vreal_t tmp_w2;
      set_zero(tmp_w2);
      for(int a=0; a<size; ++a){
        vcomplex_t ww=w_ptr[a];
        accum_abs2(tmp_w2, ww);
      }
      accum_add(w2_vsum, tmp_w2);

      ptr+=size;
      w_ptr+=size;
    }

    dotc_sum+=reduce(dotc_vsum);
    norm2_sum+=reduce(norm2_vsum);
    w2_sum+=reduce(w2_vsum);
  }
  //#pragma omp master
  //  {
  //    timer_dotc_norm2_calc.stop();
  //  }
  double sum_array[4]={real(dotc_sum), imag(dotc_sum), norm2_sum, w2_sum};
  ThreadManager_OpenMP::reduce_sum_global(sum_array, 4, ith, nth);
  norm2=sum_array[2];
  w_norm2=sum_array[3];
  complex_t dotc(sum_array[0], sum_array[1]);
  //#pragma omp master
  //  {
  //    timer_dotc_norm2_calc.stop();
  //  }
  return dotc;
}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::add(const int ex,
                                const AField<REALTYPE,SIMD>& w, const int ex_w){

  // |this[ex] >  += |w[w_ex] >

#ifdef DEBUG
  vout.paranoiac("AField::add (4 arg version) is called\n");
#endif

  // TODO: prefetch

  assert(ex<=m_nex);
  assert(ex_w<=w.nex());
  const size_t offset_in=m_nin_complex*m_rvol*ex_w;
  const size_t offset_out=m_nin_complex*m_rvol*ex;

  vcomplex_t *out=ptr_vcomplex() + offset_out;
  const vcomplex_t *in=w.ptr_vcomplex()+ offset_in;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int s=is; s<ns; ++s){
    vcomplex_t y=out[s];
    vcomplex_t x=in[s];
    accum_add(y,x); // y := y+x
    out[s]=y;
  }

}

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE,SIMD>::sub(const int ex,
                                const AField<REALTYPE,SIMD>& w, const int ex_w){

  // |this[ex] >  -= |w[w_ex] >

#ifdef DEBUG
  vout.paranoiac("AField::sub (4 arg version) is called\n");
#endif

  // TODO: prefetch

  assert(ex<=m_nex);
  assert(ex_w<=w.nex());
  const size_t offset_in=m_nin_complex*m_rvol*ex_w;
  const size_t offset_out=m_nin_complex*m_rvol*ex;

  vcomplex_t *out=ptr_vcomplex() + offset_out;
  const vcomplex_t *in=w.ptr_vcomplex()+ offset_in;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_rvol);
  is*=m_nin_complex;
  ns*=m_nin_complex;
  for(int s=is; s<ns; ++s){
    vcomplex_t y=out[s];
    vcomplex_t x=in[s];
    accum_sub(y,x); // y := y+x
    out[s]=y;
  }

}




//====================================================================

// explicit instanciation.
template<>
const std::string AField<double,SIMD>::class_name = "AField<double,SIMD>";
template<>
const std::string AField<float,SIMD>::class_name = "AField<float,SIMD>";

template class AField<float,SIMD>;
template class AField<double,SIMD>;


// copy between fields with different precisions.
//====================================================================
void copy(AField<double,SIMD>& v, const AField<float,SIMD> &w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  vout.crucial("AField in lib_alt_SIMD: copy from float to double is not supported");
  exit(EXIT_FAILURE);
  /*
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, v.size());

  for (int i = is; i < ns; ++i) {
    v.set(i, double(w.cmp(i)) );
  }
  */
}

//====================================================================
void copy(AField<float,SIMD>& v, const AField<double,SIMD> &w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  vout.crucial("AField in lib_alt_SIMD: copy from double to float is not supported");
  exit(EXIT_FAILURE);

  /*
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, v.size());

  for (int i = is; i < ns; ++i) {
    v.set(i, float(w.cmp(i)) );
  }
  */
}
//============================================================END=====

