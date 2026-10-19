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

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/aindex_eo-inc.h"

#include "lib_alt/Solver/aprecond_Mixedprec.h"
#include "lib_alt/Solver/asolver_Richardson.h"

#include "bridge_real.h"

#ifdef USE_MPREAL
#include "bridge_mpreal.h"
using mpfr::mpreal;
#endif

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson-tmpl.h"

template<>
const std::string Fprop_alt_Standard_eo_Richardson<AField<double, QXS>, AField<float, QXS> >::class_name
  = "Fprop_alt_Standard_eo_Richardson<AField<double,QXS>, AField<float,QXS> >";

// class instanciation
template class Fprop_alt_Standard_eo_Richardson<AField<double, QXS>, AField<float, QXS> >;


#ifdef USE_MPREAL
template<>
const std::string Fprop_alt_Standard_eo_Richardson<AField<double, QXS>, AField<mpreal, QXS> >::class_name
  = "Fprop_alt_Standard_eo_Richardson<AField<double,QXS>, AField<mpreal,QXS> >";

// class instanciation
template class Fprop_alt_Standard_eo_Richardson<AField<double, QXS>, AField<mpreal, QXS> >;
#endif

//============================================================END=====
