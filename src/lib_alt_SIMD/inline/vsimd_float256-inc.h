/*!

        @file    vsimd_float256-inc.h

        @brief   simd variables

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================
#ifndef VSIMD_T_H
#define VSIMD_T_H

namespace{
  typedef __m256 Vsimd_t;
  typedef __m256 Vsimd_t_r; // temporal hack
  constexpr int Nsimd=4;
}
#endif
