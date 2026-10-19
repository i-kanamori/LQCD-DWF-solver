/*!

        @file    aprecond_Mixedprec.cpp

        @brief   Mixed precision preconditionor

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $

 */
//====================================================================


#include <cassert>

#include "lib_alt/Solver/aprecond_Mixedprec.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/aindex_lex.h"
#include "lib_alt_SIMD/Field/aindex_eo.h"

// function template files
#include "lib_alt_SIMD/Field/afield-inc.h"

// template functions
#include "lib_alt/Solver/aprecond_Mixedprec-tmpl.h"


// explicit instanciation.
template<>
const std::string APrecond_Mixedprec<AField<double,SIMD>, AField<float,SIMD> >::class_name
                              = "APrecond_Mixedprec<AField<double,SIMD>, AField<float,SIMD> >";

template class APrecond_Mixedprec<AField<double,SIMD>, AField<float,SIMD> >;


//============================================================END=====
