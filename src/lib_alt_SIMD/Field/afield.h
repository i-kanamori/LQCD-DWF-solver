/*!

        @file    afield.h

        @brief   alternative field (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================

#ifndef AFIELD_SIMD_INCLUDED
#define AFIELD_SIMD_INCLUDED

#include <cstddef>
#include <cstdio>
#include <vector>
#include <string>

#include "lib_alt/Field/afield_base.h"  // primary template

#include "bridge_defs.h"
#include "bridge_complex.h"
#include "complexTraits.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Field/field.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//#include "lib/Tools/timer.h"

#include "lib_alt_SIMD/alt_simd_def.h"

// hack to treat
//   alinged_vecotr<vcomplex_t> and
//   alinged_vecotr<vreal_t>
// uniformly
class Field_body_base{
  public:
  virtual ~Field_body_base(){};
  virtual void* ptr_void() const = 0;   // a hack to be pure virtual class
};


//====================================================================
template<typename VTYPE>
class Field_body : public Field_body_base{
public:
  aligned_vector<VTYPE> data;
  Field_body(const size_t size): m_ptr(0) {
    if(size>0){
      data.resize(size);
      m_ptr=&(data[0]);
    }
  }
  ~Field_body(){
    m_ptr=0;
  }
  void* ptr_void() const { return reinterpret_cast<void*>(m_ptr); }
  VTYPE *ptr()const{ return m_ptr; }
protected:
  VTYPE *m_ptr;

private:
  Field_body();

};

//====================================================================
template<typename REALTYPE>
class AField<REALTYPE, SIMD>
{

public:
  typedef REALTYPE real_t;
  using vreal_t=ValType_real<REALTYPE>;
  using vcomplex_t=ValType_complex<REALTYPE>;
  using complex_t=typename vcomplex_t::scalar_type;
  using Field_vreal=Field_body<vreal_t>;
  using Field_vcomplex=Field_body<vcomplex_t>;
  using Field_real=Field_body<real_t>;
  using Field_complex=Field_body<complex_t>;

  typedef Element_type::type element_type;

  static constexpr Impl IMPL = SIMD;
  static const std::string class_name;

  //  static Timer timer_dotc_norm2;
  //  static Timer timer_dotc_norm2_calc;

protected:
  Bridge::VerboseLevel m_vl;
  element_type m_element_type=Element_type::COMPLEX;
  Field_body_base *m_ptr;
  //  vcomplex_t *m_ptr_vcomplex=0;
  //  vreal_t *m_ptr_vreal=0;

  int m_nin;
  size_t m_nvol;
  int m_nex;
  int m_ncomplex;    // 1 for real, 2 for complex
  int m_nin_complex; // nin interms of complex number: m_nin/m_ncomplex


  std::size_t m_nsize;  // counted by real
  std::size_t m_nsize_complex; // counted by complex or real (based on the type)
  std::size_t m_rsize; // reduced size (= size in terms of simd vector)
  std::size_t m_rvol;  // reduced lattice volume
  //real_t *m_field;

  int m_size_unit;               //!< unit of reduction
  int m_size_rem;                //!< nsize % m_size_unit
  std::size_t m_size_quot;       //!< nsize/m_size_unit
  std::vector<real_t> m_reduce;  //!< array for reduction


public:
  //! constructor without argument
  AField() : m_vl(CommonParameters::Vlevel()), m_element_type(Element_type::COMPLEX), m_ptr(0){ init(0, 0, 0); }

  //! constructor
  AField(const int nin, const int nvol, const int nex,
         const element_type cmpl = Element_type::COMPLEX)
    : m_vl(CommonParameters::Vlevel()), m_element_type(cmpl), m_ptr(0)
  { init(nin, nvol, nex); }

  //! copy constructor
  AField(const AField<real_t,SIMD>& w)
    : m_vl(w.m_vl), m_element_type(w.field_element_type()), m_ptr(0)
  {
    init(w.nin(), w.nvol(), w.nex());
    copy(w);
  }


  //! destructor
  ~AField(){ tidyup(); };


private:
  void init(const int nin, const int nvol, const int nex);

  void tidyup();


public:

  // resetting object.
  void reset(const int nin, const int nvol, const int nex,
             const element_type cmpl = Element_type::COMPLEX)
  {
    vout.paranoiac("AField: resetting...\n");
    if(check_size(nin,nvol,nex) && m_element_type == cmpl) return;
    m_element_type = cmpl;
    vout.paranoiac("AField: resetting...calling init\n");
    init(nin, nvol, nex);
  }

  //! alias of reset (without element_type): use reset.
  void resize(const int nin, const int nvol, const int nex)
  {
    if(check_size(nin,nvol,nex)) {return;}
    init(nin, nvol, nex);
  }

  //! returning size of inner (on site) d.o.f.
  int nin() const { return m_nin; }

  //! returning size of inner (on site) d.o.f.
  int nin_complex() const { return m_nin_complex; }

  //! returning size of site d.o.f.
  int nvol() const { return m_nvol; }

  //! returning size of extra d.o.f.
  int nex() const { return m_nex; }

  //! returning size of reduced site d.o.f.
  int rvol() const { return m_rvol; }

  //! returning element_type (real or complex).
  element_type field_element_type() const { return m_element_type; }

  //! checking size parameters.
  bool check_size(const int nin, const int nvol, const int nex) const
  {
    bool chk = true;
    if(m_nin != nin || m_nvol != nvol || m_nex != nex) chk = false;
    return chk;
  }

  bool check_size(const AField<REALTYPE,SIMD> &w) const
  {
    bool chk = true;
    if(m_nin != w.nin() || m_nvol != w.nvol() || m_nex != w.nex())
      chk = false;
    return chk;
  }


  //! reference of data element.  to be discarded.
  inline real_t& e(const int index)
  {
    real_t *m_field=(real_t*)m_ptr->ptr_void();
    return m_field[index];
  }

  //! reference of data element.  to be discarded.
  inline const real_t& e(const int index) const
  {
    const real_t *m_field=(real_t*)m_ptr->ptr_void();
    return m_field[index];
  }

  //! return the array size
  inline size_t size(void) const { return m_nsize; }

  //! return the array size
  inline size_t rsize(void) const { return m_rsize; }

  //! return the array size
  inline int size_complex(void) const { return m_nsize_complex; }

  //! return the array size
  inline int ntot(void) const { return m_nsize; }

  real_t *ptr(int i){
    real_t *m_field=(real_t*)m_ptr->ptr_void();
    return &m_field[i];
  }

  const real_t *ptr(int i) const {
    const real_t *m_field=(const real_t*)m_ptr->ptr_void();
    return &m_field[i];
  }

  // pointer as a number: do not check REAL/COMPLEX
  real_t *ptr_real() const {
    //    return (real_t*) m_field;
    return (real_t*) m_ptr->ptr_void();
  }

  complex_t *ptr_complex() const {
    //    return (complex_t*) m_field;
    return (complex_t*) m_ptr->ptr_void();
  }

  vreal_t *ptr_vreal() const {
    //    return m_ptr_vreal;
    return (vreal_t*) m_ptr->ptr_void();
  }

  vcomplex_t *ptr_vcomplex() const {
    //    return m_ptr_vcomplex;
    return (vcomplex_t*) m_ptr->ptr_void();
  }

  void *ptr_void() const {
    return m_ptr->ptr_void();
  }

  real_t cmp(const int index) const {
    const real_t *m_field=(real_t*)m_ptr->ptr_void();
    return m_field[index];
  }

  void set(const int index, const real_t a){
    assert(index>=0);
    assert(index < m_nin*m_nex*m_nvol);
    real_t *m_field=(real_t*)m_ptr->ptr_void();
    m_field[index] = a;
  }

  void set(const real_t a);

  void set(const complex_t a);

  void setc(const real_t a);

  void setc(const complex_t a);

  void copy(const AField<real_t,SIMD>& w);

  void copy(const int index,
            const AField<real_t,SIMD>& w, const int index_w, const int size);

  void copy(const int ex,
            const AField<real_t,SIMD>& w, const int ex_w);

  void axpy(const real_t, const AField<real_t,SIMD>&);

  void axpy(const complex_t, const AField<real_t,SIMD>&);

  void axpy(const int index, const real_t a,
            const AField<real_t,SIMD>& w, const int index_w,
            const int size);

  void axpy(const int ex, const real_t a,
            const AField<real_t,SIMD>& w, const int ex_w);

  void axpy(const int ex, const complex_t a,
            const AField<real_t,SIMD>& w, const int ex_w);

  void aypx(const real_t, const AField<real_t,SIMD>&);

  void aypx(const complex_t, const AField<real_t,SIMD>&);

  void scal(const real_t);

  void scal(const complex_t);

  void add(const int ex,
           const AField<real_t,SIMD>& w, const int ex_w);

  void sub(const int ex,
           const AField<real_t,SIMD>& w, const int ex_w);

  // dot as real vectors
  real_t dot(const AField<real_t,SIMD>&) const;

  // dot as complex vectors
  complex_t dotc(const AField<real_t,SIMD>&) const;

  // dot as complex vectors
  complex_t dotc(const int ex, const AField<real_t,SIMD>&, const int ex_w) const;

  complex_t dotc_and_norm2(real_t &norm2, real_t &w_norm2, const AField<real_t,SIMD> &w) const ;

  real_t norm2(void) const;

  real_t norm(void) const {
    return ::sqrt(this->norm2());
  }

};

// copy between fields with different precisions.
void copy(AField<double,SIMD>&, const AField<float,SIMD>&);
void copy(AField<float,SIMD>&, const AField<double,SIMD>&);

#endif // AFIELD_SIMD_INCLUDED
//============================================================END=====
