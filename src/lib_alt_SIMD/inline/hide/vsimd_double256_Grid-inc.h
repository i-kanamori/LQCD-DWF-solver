#ifndef VSIMD_T_H
#define VSIMD_T_H

//#include "GridSimd.h"

namespace{
  // typedef __m256d Vsimd_t;
  using namespace Grid;
  typedef vComplexD Vsimd_t;
  typedef vRealD Vsimd_t_r; // temporal hack
  const int Nsimd=2;
}
#endif
