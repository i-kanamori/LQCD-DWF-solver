/*!
      @file    asolver_BiCGStab.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 1983 $
*/

#include "lib_alt/Solver/asolver_BiCGStab.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/afield-inc.h"


#include "lib_alt/Solver/asolver_BiCGStab-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD>.
template<>
const std::string ASolver_BiCGStab<AField<double,SIMD> >::class_name
                              = "ASolver_BiCGStab<Afield<double,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab<AField_dev<double,SIMD> >::register_factory();
}
#endif

template class ASolver_BiCGStab<AField<double,SIMD> >;

//====================================================================
// explicit instanciation for AField<float,SIMD>.
template<>
const std::string ASolver_BiCGStab<AField<float,SIMD> >::class_name
                              = "ASolver_BiCGStab<AField<float,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab<AField_dev<float,SIMD> >::register_factory();
}
#endif

template class ASolver_BiCGStab<AField<float,SIMD> >;


//====================================================================
//============================================================END=====
