/*!
        @file    $Id: valtype_512.h  #$

        @brief   variables for AVX-512

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */

#ifndef VALTYPE_512_h
#define VALTYPE_512_h

#include "immintrin.h"

namespace Double{ // force include to inside a different namespace
#ifdef VSIMD_T_H // dirty hack
#undef VSIMD_T_H
#endif
#include "vsimd_double512-inc.h"
} // Double

namespace Float{ // force include to inside a different namespace
#ifdef VSIMD_T_H // dirty hack
#undef VSIMD_T_H
#endif
#include "vsimd_float512-inc.h"
} // Float

#endif
