/*!
 Time-stamp: <2023-10-19 16:50:31 kanamori>

        @file    $Id: aindex_eo.h  #$

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */


#ifndef AINDEX_EO_SIMD_H
#define AINDEX_EO_SIMD_H

#include <vector>

#include "lib_alt/Field/aindex_eo_base.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_SIMD/alt_simd_def.h"
#include "afield.h"

#ifdef IMPLTYPE_ALT_GRIDSIMD
#include "Imp_gridsimd/aindex_eo_gridsimd.h"
#elif defined(IMPLTYPE_ALT_IMPL1)
#include "Imp_impl1/aindex_eo_impl1.h"
#else
#error no type for index simd layout is specified
#endif

#endif // AINDEX_EO_SIMD_H
