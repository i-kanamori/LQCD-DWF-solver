// Time-stamp: <2023-10-20 00:25:51 kanamori>
/*!
        @file    $Id: aindex_lex.h  #$

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */

#ifndef AINDEX_LEX_SIMD_H
#define AINDEX_LEX_SIMD_H

#include "lib_alt/Field/aindex_lex_base.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_SIMD/alt_simd_def.h"

#ifdef IMPLTYPE_ALT_GRIDSIMD
#include "lib_alt_SIMD/Field/Imp_gridsimd/aindex_lex_gridsimd.h"
#elif defined(IMPLTYPE_ALT_IMPL1)
#include "lib_alt_SIMD/Field/Imp_impl1/aindex_lex_impl1.h"
#else
#error no type for index simd layout is specified
#endif

#endif
