/*!
        @file    $Id:: shiftAField_lex_float.h #$

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: $

        @date    $LastChangedDate:: 2017-12-21 10:39:00 #$

        @version $LastChangedRevision: 1687 $
*/

#include "shiftAField_lex.h"
#include "afield.h"
#include "lib/ResourceManager/threadManager.h"

#define  VLEN    VLENS
#define  VLEN2   VLENS2

// inline function files
#include "lib_alt_SIMD/inline/vector_arith-inc.h"
#include "lib_alt_SIMD/inline/vector_shift_float-inc.h"

// function template files
#include "afield-inc.h"
#ifdef IMPLTYPE_ALT_GRIDSIMD
#error gridsimd is not ready for shiftAField
//#include "Imp_gridsimd/shiftAField_lex_impl1-tmpl.h"
#elif defined(IMPLTYPE_ALT_IMPL1)
#include "Imp_impl1/shiftAField_lex_impl1-tmpl.h"
#else
#error no type for index simd layout is specified
#endif


template<>
const std::string ShiftAField_lex<AField<float,SIMD> >::class_name
                              = "ShiftAField_lex<AField<float,SIMD> >";

// explicit instanciation.
template class ShiftAField_lex<AField<float,SIMD> >;

//============================================================END=====
