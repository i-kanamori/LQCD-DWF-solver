/*!
 Time-stamp: <2023-10-19 16:48:56 kanamori>

        @file    $Id: aindex_eo-inc.h  #$

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */


#ifndef AINDEX_EO_SIMD_INC_H
#define AINDEX_EO_SIMD_INC_H

#include "aindex_eo.h"

#ifdef IMPLTYPE_ALT_GRIDSIMD
#include "Imp_gridsimd/aindex_eo_gridsimd-inc.h"
#elif defined(IMPLTYPE_ALT_IMPL1)
#include "Imp_impl1/aindex_eo_impl1-inc.h"
#else
#error no type for index simd layout is specified
#endif

#endif // INDEX_EO_ALT_INC_H
