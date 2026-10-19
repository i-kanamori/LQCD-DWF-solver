// some arithmetics : using Grid
#ifndef VECTOT_ARITH_GRID_H
#define VECTOT_ARITH_GRID_H

// hack to fix a funny bug in Grid
#ifdef IMCI
#define FLIP_timesI
#endif

#ifdef FLIP_timesI
#define MyTimesI(x) timesMinusI((x))
#define MyTimesMinusI(x) timesI((x))
#else
#define MyTimesI(x) timesI((x))
#define MyTimesMinusI(x) timesMinusI((x))
#endif

/* Immediate value, as Quaternary (4 digits) */
#define IMM_QUART4(A,B,C,D) ((A<<6)|(B<<4)|(C<<2)|(D))
/* Immediate value, as Binary (8 digits) */
#define IMM_BIN8(A,B,C,D,E,F,G,H) ((A<<7)|(B<<6)|(C<<5)|(D<<4)|(E<<3)|(F<<2)|(G<<1)|H)


namespace {
  using namespace Grid;

  ///////////////////////////////////////
  // double
  //////////////////////////////////////

  // broadcast
  static inline void bcast(vrealD_t &va, const realD_t a){
    vsplat(va.vec,a);
  }
  static inline void bcast(vcomplexD_t &va, const complexD_t a){
    vsplat(va.vec,a);
  }

  // streaming store
  static inline void sstore(vcomplexD_t *vbuf, const vcomplexD_t &va){
    vstream((vbuf->vec), va.vec);
  }

  static inline void sstore(vrealD_t *vbuf, const vrealD_t &va){
    vstream((vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealD_t &va){
#ifdef AVX512
    va.vec.v=_mm512_setzero_pd();
#else
    vzero(va.vec);
#endif
  }
  static inline void set_zero(vcomplexD_t &va){
#ifdef AVX512
    va.vec.v=_mm512_setzero_pd();
#else
    vzero(va.vec);
#endif
  }

    // reduction
  static inline realD_t reduce(const vrealD_t &va){
    return Reduce(va.vec);
  }
  static inline complexD_t reduce(const vcomplexD_t &va){
    return Reduce(va.vec);
  }

  // add
  static inline vrealD_t add(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(va.vec+vb.vec);
  }
  static inline vcomplexD_t add(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(va.vec+vb.vec);
  }
  static inline void accum_add(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec+=vb.vec;
  }
  //  inline void operator+=(vcomplexD_t &va,const vcomplexD_t &vb){
  //    va.vec+=vb.vec;
  //  }
  
  // subtraction
  static inline vrealD_t sub(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(va.vec-vb.vec);
  }
  static inline vcomplexD_t sub(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(va.vec-vb.vec);
  }
  static inline void accum_sub(vcomplexD_t &va,const vcomplexD_t &vb){
    va.vec-=vb.vec;
  }
  //  inline void operator-=(vcomplexD_t &va,const vcomplexD_t &vb){
  //    va.vec-=vb.vec;
  //  }

  // mult
  static inline vrealD_t mult(const vrealD_t &va, const vrealD_t &vb){
    return vrealD_t(va.vec*vb.vec);
  }
  static inline vcomplexD_t cmult(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(va.vec*vb.vec);
  }

#if defined(AVX512)
  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    // slight modification of a*b in Grid
    __m512d a_real = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0xFF );
    a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( vb.vec.v, 0x55 ) ); 
    out.vec.v =_mm512_fmaddsub_pd( a_real, vb.vec.v, a_imag );
  }
#else
  static inline void set_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    out=vcomplexD_t(va.vec*vb.vec);
  }
#endif

  // mult:  conj(a)*b
#if defined(AVX512)
  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    // slight modification of a*b in Grid
    __m512d a_real = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0xFF );
    a_imag = _mm512_mul_pd( a_imag, _mm512_permute_pd( vb.vec.v, 0x55 ) ); 
    out.vec.v =_mm512_fmsubadd_pd( a_real, vb.vec.v, a_imag );
  }
  static inline vcomplexD_t cmult_conj(const vcomplexD_t &va, const vcomplexD_t &vb){
    vcomplexD_t out;
    set_cmult_conj(out,va,vb);
    return out;
  }
  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    __m512d a_real = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0x00 );
    __m512d a_imag = _mm512_shuffle_pd( va.vec.v, va.vec.v, 0xFF );
    a_imag = _mm512_fmsubadd_pd( a_imag, _mm512_permute_pd( vb.vec.v, 0x55 ) , out.vec.v);
    out.vec.v =_mm512_fmsubadd_pd( a_real, vb.vec.v, a_imag );
  }

#else
  static inline vcomplexD_t cmult_conj(const vcomplexD_t &va, const vcomplexD_t &vb){
    return vcomplexD_t(conjugate(va.vec)*vb.vec);
  }
  static inline void set_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    out.vec=conjugate(va.vec)*vb.vec;
  }
  static inline void accum_cmult_conj(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    out.vec=out.vec+conjugate(va.vec)*vb.vec;
  }
#endif
  
  // sqrt(-1)
  static inline vcomplexD_t TimesI(const vcomplexD_t &in){
    return vcomplexD_t(MyTimesI(in.vec));
  }
  static inline vcomplexD_t TimesMinusI(const vcomplexD_t &in){
    return vcomplexD_t(MyTimesMinusI(in.vec));
  }
#if defined(AVX512)
  static void inline accum_TimesI(vcomplexD_t &out, const vcomplexD_t &in){
      __m512d tmp=_mm512_shuffle_pd(in.vec.v,in.vec.v,0x55);
      tmp=_mm512_mask_sub_pd(tmp,0x55,_mm512_setzero_pd(),tmp); 
      out.vec.v=_mm512_add_pd(out.vec.v,tmp);
  }
  static void inline accum_TimesMinusI(vcomplexD_t &out, const vcomplexD_t &in){
      __m512d tmp=_mm512_shuffle_pd(in.vec.v,in.vec.v,0x55);
      tmp=_mm512_mask_sub_pd(tmp,0x55,_mm512_setzero_pd(),tmp); 
      out.vec.v=_mm512_sub_pd(out.vec.v,tmp);
  }
#else
  static inline void accum_TimesI(vcomplexD_t &out, const vcomplexD_t &in){
    out.vec += MyTimesI(in.vec);
  }
  static inline void accum_TimesMinusI(vcomplexD_t &out, const vcomplexD_t &in){
    out.vec += MyTimesMinusI(in.vec);
  }
#endif

  // axpy: y := ax + y
  static inline void vaxpy(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    vy.vec+=va.vec*vx.vec;
  }
  static inline void vaxpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec+=va.vec*vx.vec;
  }
  static inline void vaxpy(vcomplexD_t &vy, const vrealD_t &va, const vcomplexD_t &vx){
    Grid::Optimization::Mult mult;
    mult.mac(vy.vec.v, va.vec.v, vx.vec.v);
    //    vy.vec+=va.vec*vx.vec;
  }

  static inline void axpy(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vaxpy(vy,va,vx);
  }
  // out = out+va*vb
  static inline void accum_cmult(vcomplexD_t &out, const vcomplexD_t &va, const vcomplexD_t &vb){
    vaxpy(out,va,vb);
  }



  // aypx: y := a*y+x
  static inline void vaypx(vcomplexD_t &vy, const vcomplexD_t &va, const vcomplexD_t &vx){
    vy.vec=va.vec*vy.vec + vx.vec;
  }
  static inline void vaypx(vrealD_t &vy, const vrealD_t &va, const vrealD_t &vx){
    vy.vec=va.vec*vy.vec + vx.vec;
  }

  // scal
  static inline void scal(vcomplexD_t &vx, const vcomplexD_t &va){
    vx.vec*=va.vec;
  }
  static inline void scal(vrealD_t &vx, const vrealD_t &va){
    vx.vec*=va.vec;
  }
  static inline void scal(vcomplexD_t &vx, const vrealD_t &va){
    vrealD_t *tmp=(vrealD_t*)&vx;
    scal(*tmp,va);
  }

  // (re im) (re im)... --> (re re) (re re)...
#if defined(AVX512)
  static inline void complex2rere(vcomplexD_t &vx){
    vx.vec.v=_mm512_shuffle_pd( vx.vec.v, vx.vec.v, 0x00 );
  }
#else
  static inline void complex2rere(vcomplexD_t &vx){
    vx.vec.v=_mm256_movedup_pd(vx.vec.v); // AVX2
    //    vx.vec.v=_mm256_shuffle_pd(vx.vec.v,vx.vec.v,0x0); // AVX1
    // other
    //    double *tmp=(double*)&vx;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }
  }
#endif

  // (re im) (re im)... --> (im im) (im im)
#if defined(AVX512)
  static inline void complex2imim(vcomplexD_t &vx){
    vx.vec.v=_mm512_shuffle_pd( vx.vec.v, vx.vec.v, 0xFF );
  }
#else
  static inline void complex2imim(vcomplexD_t &vx){
    vx.vec.v=_mm256_shuffle_pd(vx.vec.v,vx.vec.v,0xF); // AVX1/AVX2
    // other
    //    double *tmp=(double*)&vx;
    //    for(int i=0; i<vcomplexD_t::Nsimd(); i++){
    //      tmp[2*i]=tmp[2*i+1];
    //    }
  }
#endif  
  
  ///////////////////////////////////////
  // float
  //////////////////////////////////////

  // broadcast
  static inline void bcast(vrealF_t &va, const realF_t a){
    vsplat(va.vec,a);
  }
  static inline void bcast(vcomplexF_t &va, const complexF_t a){
    vsplat(va.vec,a);
  }

  // streaming store
  static inline void sstore(vcomplexF_t *vbuf, const vcomplexF_t &va){
    vstream((vbuf->vec), va.vec);
  }

  static inline void sstore(vrealF_t *vbuf, const vrealF_t &va){
    vstream((vbuf->vec), va.vec);
  }

  // set zero
  static inline void set_zero(vrealF_t &va){
#ifdef AVX512
    va.vec.v=_mm512_setzero_ps();
#else
    vzero(va.vec);
#endif
  }
  static inline void set_zero(vcomplexF_t &va){
#ifdef AVX512
   va.vec.v=_mm512_setzero_ps();
#else
    vzero(va.vec);
#endif
  }

  // reduction
  static inline realF_t reduce(const vrealF_t &va){
    return Reduce(va.vec);
  }
  static inline complexF_t reduce(const vcomplexF_t &va){
    return Reduce(va.vec);
  }

  // add
  static inline vrealF_t add(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(va.vec+vb.vec);
  }
  static inline vcomplexF_t add(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(va.vec+vb.vec);
  }
  static inline void accum_add(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec+=vb.vec;
  }
  //  inline void operator+=(vcomplexF_t &va,const vcomplexF_t &vb){
  //    va.vec+=vb.vec;
  //  }

  // subtraction
  static inline vrealF_t sub(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(va.vec-vb.vec);
  }
  static inline vcomplexF_t sub(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(va.vec-vb.vec);
  }
  inline void accum_sub(vcomplexF_t &va,const vcomplexF_t &vb){
    va.vec-=vb.vec;
  }
  //  inline void operator-=(vcomplexF_t &va,const vcomplexF_t &vb){
  //    va.vec-=vb.vec;
  //  }

  // mult
  static inline vrealF_t mult(const vrealF_t &va, const vrealF_t &vb){
    return vrealF_t(va.vec*vb.vec);
  }
  static inline vcomplexF_t cmult(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(va.vec*vb.vec);
  }


#if defined(AVX512)
  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    // slight modification of a*b in Grid
    __m512 a_real = _mm512_moveldup_ps( va.vec.v ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( va.vec.v ); // Ai Ai
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( vb.vec.v, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    out.vec.v =_mm512_fmaddsub_ps( a_real, vb.vec.v, a_imag ); // (Ar Br, Ar Bi)-+(Ai Bi, Ai Br) =  (ArBr-AiBi , ArBi+AiBr)

  }
#else
  static inline void set_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    out=vcomplexF_t(va.vec*vb.vec);
  }
#endif

  // mult:  conj(a)*b
  static inline vcomplexF_t cmult_conj(const vcomplexF_t &va, const vcomplexF_t &vb){
    return vcomplexF_t(conjugate(va.vec)*vb.vec);
  }

#ifdef AVX512
  static void inline set_cmult_conj(vcomplexF_t &out, vcomplexF_t a, vcomplexF_t b){
  // slight modification of a*b in Grid
    __m512 a_real = _mm512_moveldup_ps( a.vec.v ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( a.vec.v ); // Ai Ai
    a_imag = _mm512_mul_ps( a_imag, _mm512_permute_ps( b.vec.v, 0xB1 ) );  // (Ai, Ai) * (Bi, Br) = Ai Bi, Ai Br
    out.vec.v =_mm512_fmsubadd_ps( a_real, b.vec.v, a_imag ); // (Ar Br, Ar Bi)+-(Ai Bi, Ai Br) = (ArBr+AiBi, ArBi-AiBr)
  }

  static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    __m512 a_real = _mm512_moveldup_ps( va.vec.v ); // Ar Ar
    __m512 a_imag = _mm512_movehdup_ps( va.vec.v ); // Ai Ai
    a_imag = _mm512_fmsubadd_ps( a_imag, _mm512_permute_ps( vb.vec.v, 0xB1 ), out.vec.v );  // (Ai, Ai) * (Bi, Br) +-(Cr,Ci) = Ai Bi+Cr, Ai Br-Ci
    out.vec.v =_mm512_fmsubadd_ps( a_real, vb.vec.v, a_imag ); // (Ar Br, Ar Bi)+-(Ai Bi +Cr, Ai Br -Ci) = (ArBr+AiBi +Cr, ArBi-AiBr +Ci)
  }

#else
  static inline void set_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    out.vec=conjugate(va.vec)*vb.vec;
  }
  static inline void accum_cmult_conj(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    out.vec=out.vec+conjugate(va.vec)*vb.vec;
  }
#endif
     
  // sqrt(-1)
  static inline vcomplexF_t TimesI(const vcomplexF_t &in){
    return vcomplexF_t(MyTimesI(in.vec));
  }
  static inline vcomplexF_t TimesMinusI(const vcomplexF_t &in){
    return vcomplexF_t(MyTimesMinusI(in.vec));
  }
#ifdef AVX512
  static void inline accum_TimesI(vcomplexF_t &out, const vcomplexF_t &in){
      __m512 tmp=_mm512_shuffle_ps(in.vec.v,in.vec.v,IMM_QUART4(2,3,0,1));
      tmp=_mm512_mask_sub_ps(tmp,0x5555,_mm512_setzero_ps(),tmp); 
      out.vec.v=_mm512_add_ps(out.vec.v, tmp);
  }
  static void inline accum_TimesMinusI(vcomplexF_t &out, const vcomplexF_t &in){
      __m512 tmp=_mm512_shuffle_ps(in.vec.v,in.vec.v,IMM_QUART4(2,3,0,1));
      tmp=_mm512_mask_sub_ps(tmp,0x5555,_mm512_setzero_ps(),tmp); 
      out.vec.v=_mm512_sub_ps(out.vec.v, tmp);
  }
#else
  static inline void accum_TimesI(vcomplexF_t &out, const vcomplexF_t &in){
    out.vec += MyTimesI(in.vec);
  }
  static inline void accum_TimesMinusI(vcomplexF_t &out, const vcomplexF_t &in){
    out.vec += MyTimesMinusI(in.vec);
  }
#endif

  // axpy: y := a*x+y
  static inline void vaxpy(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    vy.vec+=va.vec*vx.vec;
  }
  static inline void vaxpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec+=va.vec*vx.vec;
  }
  //  static inline void axpy(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
  //    vaxpy(vy,va,vx);
  //  }
  static inline void vaxpy(vcomplexF_t &vy, const vrealF_t &va, const vcomplexF_t &vx){
    Grid::Optimization::Mult mult;
    mult.mac(vy.vec.v, va.vec.v, vx.vec.v);
    //    vy.vec+=va.vec*vx.vec;
  }
  // out = out+va*vb
  static inline void accum_cmult(vcomplexF_t &out, const vcomplexF_t &va, const vcomplexF_t &vb){
    vaxpy(out,va,vb);
  }


  // aypx: y := a*y+x
  static inline void vaypx(vcomplexF_t &vy, const vcomplexF_t &va, const vcomplexF_t &vx){
    vy.vec=va.vec*vy.vec + vx.vec;
  }
  static inline void vaypx(vrealF_t &vy, const vrealF_t &va, const vrealF_t &vx){
    vy.vec=va.vec*vy.vec + vx.vec;
  }



  // scal
  static inline void scal(vcomplexF_t &vx, const vcomplexF_t &va){
    vx.vec*=va.vec;
  }
  static inline void scal(vrealF_t &vx, const vrealF_t &va){
    vx.vec*=va.vec;
  }
  static inline void scal(vcomplexF_t &vx, const vrealF_t &va){
    vrealF_t *tmp=(vrealF_t*)&vx;
    scal(*tmp,va);
  }

  // (re im) (re im)... --> (re re) (re re)
#if defined(AVX512)
  static inline void complex2rere(vcomplexF_t &vx){
    vx.vec.v=_mm512_moveldup_ps( vx.vec.v );
  }
#else
  static inline void complex2rere(vcomplexF_t &vx){
    vx.vec.v=_mm256_moveldup_ps(vx.vec.v); // AVX2
    //    vx.vec.v=_mm256_shuffle_ps(vx.vec.v, vx.vec.v, IMM_QUART4(2,2,0,0)); // AVX1
    // other
    //    float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i+1]=tmp[2*i];
    //    }
  }
#endif

  // (re im) (re im)... --> (im im) (im im)
#if defined(AVX512)
  static inline void complex2imim(vcomplexF_t &vx){
    vx.vec.v=_mm512_movehdup_ps( vx.vec.v );
  }
#else
  static inline void complex2imim(vcomplexF_t &vx){
    vx.vec.v=_mm256_movehdup_ps(vx.vec.v);  // AVX2
    //    vx.vec.v=_mm256_shuffle_ps(vx.vec.v, vx.vec.v, IMM_QUART4(3,3,1,1)); // AVX1
    // other
    //    float *tmp=(float*)&vx;
    //    for(int i=0; i<vcomplexF_t::Nsimd(); i++){
    //      tmp[2*i]=tmp[2*i+1];
    //    }
  }
#endif
} // anonymous

#endif
