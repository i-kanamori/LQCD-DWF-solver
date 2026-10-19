/*!
      @file    fprop_alt_Standard_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2019-12-01 11:57:54 #$
      @version $LastChangedRevision: 2031 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"

#include "lib/ResourceManager/threadManager.h"
//#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/afield-inc.h"

#include "lib_alt_SIMD/Field/aindex_lex.h"
#include "lib_alt_SIMD/Field/aindex_eo.h"
#include "lib_alt_SIMD/Field/aindex_eo-inc.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD>.
template<>
const std::string Fprop_alt_Standard_eo<AField<double,SIMD> >::class_name
                          = "Fprop_alt_Standard_eo<AField<double,SIMD> >";


template class Fprop_alt_Standard_eo<AField<double,SIMD> >;

//====================================================================
// explicit instanciation for AField<float,SIMD>.
template<>
const std::string Fprop_alt_Standard_eo<AField<float,SIMD> >::class_name
                           = "Fprop_alt_Standard_eo<AField<float,SIMD> >";


template class Fprop_alt_Standard_eo<AField<float,SIMD> >;

//============================================================END=====
