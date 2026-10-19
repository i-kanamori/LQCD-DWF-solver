/*!

        @file    simd_layout.h

        @brief   4-dim layout of simd usage

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
#ifndef SIMD_LAYOUT_H
#define SIMD_LAYOUT_H

#if defined(IMPLTYPE_ALT_IMPL1)

namespace simd_layout{
#if defined(HHHH_256_HHHH)
  static constexpr int vlen_dcomplex[4]={2,1,1,1};
  static constexpr int vlen_fcomplex[4]={4,1,1,1};
  static constexpr int vlen_dreal[4]={4,1,1,1};
  static constexpr int vlen_freal[4]={8,1,1,1};
#define  VLEND   4
#define  VLEND2  2
#define  VLENS   8
#define  VLENS2  4
#elif defined(HHHH_512_HHHH)
  static constexpr int vlen_dcomplex[4]={4,1,1,1};
  static constexpr int vlen_fcomplex[4]={8,1,1,1};
  static constexpr int vlen_dreal[4]={8,1,1,1};
  static constexpr int vlen_freal[4]={16,1,1,1};
#define  VLEND   8
#define  VLEND2  4
#define  VLENS   16
#define  VLENS2  8
#endif
}

#elif defined(IMPLTYPE_ALT_GRIDSIMD)

namespace simd_layout{
#if defined(HHHH_256_HHHH)
  static constexpr int vlen_dcomplex[4]={1,1,1,2};
  static constexpr int vlen_fcomplex[4]={1,1,2,2};
  static constexpr int vlen_dreal[4]={1,1,2,2};
  static constexpr int vlen_freal[4]={1,2,2,2};
#define  VLEND   4
#define  VLEND2  2
#define  VLENS   8
#define  VLENS2  4
#elif defined(HHHH_512_HHHH)
  static constexpr int vlen_dcomplex[4]={1,1,2,2};
  static constexpr int vlen_fcomplex[4]={1,2,2,2};
  static constexpr int vlen_dreal[4]={1,2,2,2};
  static constexpr int vlen_freal[4]={2,2,2,2};
#define  VLEND   8
#define  VLEND2  4
#define  VLENS   16
#define  VLENS2  8
#endif
}
#endif


template<typename REAL_TYPE> 
struct simd_layout_complex{
  static inline constexpr int vlen_x() {return vlen(0);}
  static inline constexpr int vlen_y() {return vlen(1);}
  static inline constexpr int vlen_z() {return vlen(2);}
  static inline constexpr int vlen_t() {return vlen(3);}
  static inline int vlen(int mu);
  static inline bool has_simd(int mu) {return (vlen(mu)>1);}
  static inline constexpr bool has_simd_x() {return (vlen(0)>1);}
  static inline constexpr bool has_simd_y() {return (vlen(1)>1);}
  static inline constexpr bool has_simd_z() {return (vlen(2)>1);}
  static inline constexpr bool has_simd_t() {return (vlen(3)>1);}
};

template<typename REAL_TYPE> 
struct simd_layout_real{
  static inline constexpr int vlen_x() {return vlen(0);}
  static inline constexpr int vlen_y() {return vlen(1);}
  static inline constexpr int vlen_z() {return vlen(2);}
  static inline constexpr int vlen_t() {return vlen(3);}
  static inline int vlen(int mu);
  static inline bool has_simd(int mu) {return (vlen(mu)>1);}
  static inline constexpr bool has_simd_x() {return (vlen(0)>1);}
  static inline constexpr bool has_simd_y() {return (vlen(1)>1);}
  static inline constexpr bool has_simd_z() {return (vlen(2)>1);}
  static inline constexpr bool has_simd_t() {return (vlen(3)>1);}
};

template<>
inline int simd_layout_complex<double>::vlen(int mu) {return simd_layout::vlen_dcomplex[mu];}

template<>
inline int simd_layout_complex<float>::vlen(int mu) {return simd_layout::vlen_fcomplex[mu];}

template<>
inline int simd_layout_real<double>::vlen(int mu) {return simd_layout::vlen_dreal[mu];}

template<>
inline int simd_layout_real<float>::vlen(int mu) {return simd_layout::vlen_freal[mu];}

#endif
