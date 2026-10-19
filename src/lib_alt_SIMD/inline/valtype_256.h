/*!
        @file    $Id: valtype_256.h  #$

        @brief   variables for AVX/AVX2

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */

#ifndef VALTYPE_256_h
#define VALTYPE_256_h

#include "immintrin.h"
namespace Double{ // force include to inside a different namespace
#ifdef VSIMD_T_H // dirty hack
#undef VSIMD_T_H
#endif
#include "vsimd_double256-inc.h"
} // Double

namespace Float{ // force include to inside a different namespace
#ifdef VSIMD_T_H // dirty hack
#undef VSIMD_T_H
#endif
#include "vsimd_float256-inc.h"
} // Float

#endif
