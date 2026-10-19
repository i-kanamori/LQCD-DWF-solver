// Time-stamp: <2023-10-19 16:49:10 kanamori>
/*!
        @file    $Id: aindex_lex.cpp  #$

        @brief

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */

#include "lib_alt_SIMD/alt_simd_def.h"
#include "aindex_lex.h"

template class AIndex_lex<double,SIMD>;
template class AIndex_lex<float,SIMD>;


