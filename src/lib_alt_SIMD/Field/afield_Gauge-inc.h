// afield_Gauge-inc.h :  function templates [17 Dec 2017 H.Matsufuru]
//    lib_alt_SIMD version: [1 Jun 2020 I.Kanamori]
#ifndef AFIELD_GAUGE_INC_H
#define AFIELD_GAUGE_INC_H



//====================================================================
template <typename AFIELD>
void multadd_Gnn(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw,
                 const typename AFIELD::real_T a)
{
  // u[exu] += v[exv]*w[exw]
  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();
  vcomplex_t *up=u.ptr_vcomplex();
  //vreal_t* up=u.ptr_vreal();  // need to be real to use axpy as real arithmetic
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();

  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gnn(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*i+k];
          vcomplex_t bb=wp[Nc*k+j];
          accum_cmult(tmp, aa, bb);  // out+=v[i,k]*u[k,j]
        }
        axpy(out, va, tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }
}

//====================================================================
template <typename AFIELD>
void mult_Gnn(AFIELD& u, const int exu,
              const AFIELD& v, const int exv,
              const AFIELD& w, const int exw)
{
  // u = v*w
  // i.e.,  u[ij] = v[ik]*w[kj]

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_mult_Gnn(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t out;
        set_zero(out);
        for(int k=0; k<Nc; ++k){
          vcomplex_t a=vp[Nc*i+k];
          vcomplex_t b=wp[Nc*k+j];
          accum_cmult(out, a, b);  // out+=v[i,k]*w[k,j]
        }
        up[Nc*i+j]=out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void multadd_Gnd(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw,
                 const typename AFIELD::real_t a)
{
  // u += a * v*w^dagg
  // i.e.,  u[ij] += a* v[ik]* conj(w[jk])

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc(); 
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  //vreal_t* up=u.ptr_vreal();  // need to be real to use axpy as real arithmetic
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gnd(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*i+k];
          vcomplex_t bb=wp[Nc*j+k];
          accum_cmult_conj(tmp, bb, aa);  // out+=a*conj(b)
        }
        axpy(out,va, tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}



//====================================================================
template <typename AFIELD>
void multadd_Gnd(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw)
{
  // u += v*w^dagg
  // i.e.,  u[ij] += v[ik]* conj(w[jk])

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gnd(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*i+k];
          vcomplex_t bb=wp[Nc*j+k];
          accum_cmult_conj(tmp, bb, aa);  // out+=a*conj(b)
        }
        accum_add(out,tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void multsub_Gnd(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw)
{
  // u -= v*w^dagg
  // i.e.,  u[ij] -= v[ik]* conj(w[jk])

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gnd(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*i+k];
          vcomplex_t bb=wp[Nc*j+k];
          accum_cmult_conj(tmp, bb, aa);  // out+=a*conj(b)
        }
        accum_sub(out,tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void mult_Gnd(AFIELD& u, const int exu,
              const AFIELD& v, const int exv,
              const AFIELD& w, const int exw)
{
  // u = v*w^dagg
  // i.e.,  u[ij] = v[ik]* conj(w[jk])

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_mult_Gnd(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t out;
        set_zero(out);
        for(int k=0; k<Nc; ++k){
          vcomplex_t a=vp[Nc*i+k];
          vcomplex_t b=wp[Nc*j+k];
          accum_cmult_conj(out, b, a);  // out+=a*conj(b)
        }
        up[Nc*i+j]=out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void multadd_Gdn(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw,
                 const typename AFIELD::real_t a)
{
  // u += a* v^dagg*w
  // i.e.,  u[ij] *=  a*conj(v[ki])* w[kj]

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();

  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  //vreal_t* up=u.ptr_vreal();  // need to be real to use axpy as real arithmetic
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();
  vreal_t va;
  bcast(va,a);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gdn(up,vp,wp,Nc)
    //vreal_t *up_vreal=(vreal_t*)up;
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        //vreal_t out=up[Nc*i+j];
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*k+i];
          vcomplex_t bb=wp[Nc*k+j];
          accum_cmult_conj(tmp, aa, bb);  // out+=conj(aa)*bb
        }
        //vreal_t tmp2=(vreal_t)tmp;
        //axpy(out,va, tmp2);
        axpy(out,va, tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}


//====================================================================
template <typename AFIELD>
void multadd_Gdn(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw)
{
  // u +=  v^dagg*w
  // i.e.,  u[ij] +=  conj(v[ki])* w[kj]

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();

  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gdn(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*k+i];
          vcomplex_t bb=wp[Nc*k+j];
          accum_cmult_conj(tmp, aa, bb);  // out+=conj(aa)*bb
        }
        accum_add(out,tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}


//====================================================================
template <typename AFIELD>
void multsub_Gdn(AFIELD& u, const int exu,
                 const AFIELD& v, const int exv,
                 const AFIELD& w, const int exw)
{
  // u -=  v^dagg*w
  // i.e.,  u[ij] -=  conj(v[ki])* w[kj]

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  typedef ValType_real<real_t> vreal_t;
  const int Nc=CommonParameters::Nc();

  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gdn(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t tmp;
        set_zero(tmp);
        vcomplex_t out=up[Nc*i+j];
        for(int k=0; k<Nc; ++k){
          vcomplex_t aa=vp[Nc*k+i];
          vcomplex_t bb=wp[Nc*k+j];
          accum_cmult_conj(tmp, aa, bb);  // out+=conj(aa)*bb
        }
        accum_sub(out,tmp);
        up[Nc*i+j] = out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void mult_Gdn(AFIELD& u, const int exu,
              const AFIELD& v, const int exv,
              const AFIELD& w, const int exw)
{
  // u = v^dagg*w
  // i.e.,  u[ij] = conj(v[ki])*w[kj]

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* up=u.ptr_vcomplex();
  const vcomplex_t* vp=v.ptr_vcomplex();
  const vcomplex_t* wp=w.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  assert(v.nex()==1);
  assert(w.nex()==1);
  up+=Nc*Nc*exu + is*u.nin_complex();
  vp+=Nc*Nc*exv + is*v.nin_complex();
  wp+=Nc*Nc*exw + is*w.nin_complex();


  for(int site = is; site < ns; ++site){

    // matrix_multadd_Gdn(up,vp,wp,Nc)
    for(int i = 0; i < Nc; ++i){
      for(int j = 0; j < Nc; ++j){
        vcomplex_t out;
        set_zero(out);
        for(int k=0; k<Nc; ++k){
          vcomplex_t a=vp[Nc*k+i];
          vcomplex_t b=wp[Nc*k+j];
          accum_cmult_conj(out, a, b);  // out+=conj(a)*b
        }
        up[Nc*i+j]=out;
      }}
    up+=u.nin_complex();
    vp+=v.nin_complex();
    wp+=w.nin_complex();
  }

}

//====================================================================
template <typename AFIELD>
void ah_G(AFIELD& u, const int ex)
{
  // u --> (u-udag)/2
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;
  typedef typename AFIELD::vreal_t vreal_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=u.rvol();

  vcomplex_t* up=u.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  // offset for "external" index and the first coodinate
  // assumes all the "extenral" index is converted into the internal
  assert(u.nex()==1);
  up+=Nc*Nc*ex + is*u.nin_complex();

  real_t half=0.5;
  vreal_t vhalf;
  bcast(vhalf, half);
  vcomplex_t vzero;
  set_zero(vzero);

  for(int site = is; site < ns; ++site){

    // matrix_ah_G(up,Nc)
    for(int i = 0; i < Nc; ++i){
      vcomplex_t diag=up[Nc*i+i];
      complex2im(diag);
      up[Nc*i+i]=diag;
    }
    for(int i = 0; i < Nc; ++i){
      for(int j = i+1; j < Nc; ++j){
        vcomplex_t uij=up[Nc*i+j];
        vcomplex_t uji=up[Nc*j+i];
        accum_sub(uij, conj(uji));
        scal(uij, vhalf);  // uij := 0.5*uij
        uji=sub(vzero,conj(uij));
        up[Nc*i+j]=uij;
        up[Nc*j+i]=uji;
      }}

    up+=u.nin_complex();

  }

}

//====================================================================
template <typename AFIELD>
void add_unit_ex(AFIELD& u, const int ex, typename AFIELD::real_t a)
{         // u = u + a * I (I: unit matrix)

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=u.rvol();
  vcomplex_t* up=u.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  assert(ex<u.nex());
  up+=(ex*Rvol + is)*u.nin_complex();

  vcomplex_t va;
  bcast(va, a);

  for(int site = is; site < ns; ++site){

    // matrix_add_unit
    for(int i=0; i<Nc; ++i){
      vcomplex_t diag=up[Nc*i+i];
      accum_add(diag,va);
      up[Nc*i+i]=diag;
    }

    up+=u.nin_complex();

  }

}

//====================================================================
template <typename AFIELD>
void xI(AFIELD& v, const int ex)
{
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::vcomplex_t vcomplex_t;

  if(v.field_element_type() != Element_type::COMPLEX){
    vout.general("%s: xI is not relevant opearation\n",
                 v.class_name.c_str());
    return;
  }

  size_t Rvol=v.rvol();
  size_t ANin=v.nin_complex();

  int ith, nth, is, ns;
  size_t size=Rvol*ANin;
  set_threadtask(ith, nth, is, ns, Rvol);
  is*=ANin;
  ns*=ANin;
  vcomplex_t *pv=v.ptr_vcomplex() + ex*size;

  for(int i = is; i < ns; ++i){
    vcomplex_t tmp=pv[i];
    pv[i]=TimesI(tmp);
  }

}



//====================================================================
template <typename AFIELD>
void conjg(AFIELD& v)
{ // complex conjugate ( NOT Hermitian conj.)
  if(v.field_element_type() == Element_type::REAL) // do nothing
    return;

  typedef typename AFIELD::real_t real_t;
  typedef ValType_complex<real_t> vcomplex_t;
  const int Nc=CommonParameters::Nc();
  const int Rvol=v.rvol();

  vcomplex_t* vp=v.ptr_vcomplex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Rvol);

  const int ANex=v.nex();
  for(int aex=0; aex<ANex; ++aex){
    vcomplex_t* vp=v.ptr_vcomplex() + v.nin_complex()*(is+aex*v.rvol());
    for(int site = is; site < ns; ++site){
      for(int ij=0; ij<Nc*Nc; ++ij){
        vp[ij]=conj(vp[ij]);
      }
      vp+=v.nin_complex();
    }
  }

}



/*
  setting boundary condition for gauge field
  uses
    set_bondary_?p and apply_boundary
  from vector_shift_*-inc.h
 */
template<typename AFIELD>
void set_boundary_config(AFIELD& U, const int mu, const typename AFIELD::real_t bc, const int *rsize){

 typedef typename AFIELD::real_t real_t;
 typedef typename AFIELD::vreal_t vreal_t;
 typedef typename AFIELD::vcomplex_t vcomplex_t;
 assert(CommonParameters::Ndim()==4);

 vout.detailed("applying boundary condition: mu=%d, bc=%f\n",mu,bc);

 // check if this node is relevant
 int ipe_mu=Communicator::ipe(mu);
 int npe_mu=Communicator::npe(mu);

 if(ipe_mu!=npe_mu-1){
   return;
 }
 static const int DirX=Direction::XDIR;
 static const int DirY=Direction::YDIR;
 static const int DirZ=Direction::ZDIR;
 static const int DirT=Direction::TDIR;
 const int Rx=rsize[DirX];
 const int Ry=rsize[DirY];
 const int Rz=rsize[DirZ];
 const int Rt=rsize[DirT];
 const int Nc=CommonParameters::Nc();
 int ith, nth, is, ns;

 // gauge field
 vcomplex_t *u=U.ptr_vcomplex();

#pragma omp barrier

 // apply the condition in each direction
 if(mu==DirX) {
   set_threadtask(ith, nth, is, ns, Ry*Rz*Rt);
   vreal_t vbc;
   set_boundary_xp(vbc, bc);
   const int offset=Rx-1;
   for(int i=is; i<ns; ++i){
    const int ss=i*Rx+offset;
     vcomplex_t *inout=u+Nc*Nc*(ss*4+mu);
     for(int ab=0; ab<Nc*Nc; ab++){
       apply_boundary(inout[ab], vbc);
     }
   }
 } else if (mu==DirY) {
   set_threadtask(ith, nth, is, ns, Rx*Rz*Rt);
   vreal_t vbc;
   set_boundary_yp(vbc, bc);
   const int offset=Rx*(Ry-1);
   for(int i=is; i<ns; ++i){
     int tmp=i/Rx;
     const int ss= i + Rx*tmp*(Ry-1) + offset;
     vcomplex_t *inout=u+Nc*Nc*(ss*4+mu);
     for(int ab=0; ab<Nc*Nc; ab++){
       apply_boundary(inout[ab], vbc);
     }
   }
 } else if (mu==DirZ) {
   set_threadtask(ith, nth, is, ns, Rx*Ry*Rt);
   vreal_t vbc;
   set_boundary_zp(vbc, bc);
   const int offset=Rx*Ry*(Rz-1);
   for(int i=is; i<ns; ++i){
     const int ss= i%(Rx*Ry) + Rx*Ry*Rz*(i/(Rx*Ry))+offset;
     vcomplex_t *inout=u+Nc*Nc*(ss*4+mu);
     for(int ab=0; ab<Nc*Nc; ab++){
       apply_boundary(inout[ab], vbc);
     }
   }
 } else if (mu==DirT) {
   set_threadtask(ith, nth, is, ns, Rx*Ry*Rz);
   vreal_t vbc;
   set_boundary_tp(vbc, bc);
   const int offset=Rx*Ry*Rz*(Rt-1);
   for(int i=is; i<ns; ++i){
     const int ss= i+offset;
     vcomplex_t *inout=u+Nc*Nc*(ss*4+mu);
     for(int ab=0; ab<Nc*Nc; ab++){
       apply_boundary(inout[ab], vbc);
     }
   }
 } else {
   vout.crucial("cannot happen, bad mu=%d\n",mu);
   abort();
 }
#pragma omp barrier
}


//============================================================END=====
#endif
