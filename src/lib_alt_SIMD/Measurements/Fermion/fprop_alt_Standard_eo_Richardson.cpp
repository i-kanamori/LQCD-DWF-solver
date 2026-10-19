/*!
      @file    fprop_alt_Standard_eo_Richardson.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"

#include "lib/ResourceManager/threadManager.h"
//#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/afield-inc.h"

#include "lib_alt_SIMD/Field/aindex_lex.h"
#include "lib_alt_SIMD/Field/aindex_eo.h"
#include "lib_alt_SIMD/Field/aindex_eo-inc.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_Richardson.h"


#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson-tmpl.h"

template<>
const std::string Fprop_alt_Standard_eo_Richardson<AField<double, SIMD>, AField<float, SIMD> >::class_name
  = "Fprop_alt_Standard_eo_Richardson<AField<double,SIMD>, AField<float,SIMD> >";

// class instanciation
template class Fprop_alt_Standard_eo_Richardson<AField<double, SIMD>, AField<float, SIMD> >;

//============================================================END=====
