/*!
      @file    asolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 1983 $
*/

#include "bridge_real.h"
#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD>.
template<>
const std::string ASolver_CG<AField<double,SIMD> >::class_name
                              = "ASolver_CG<Afield<double,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<double,SIMD> >::register_factory();
}
#endif

template class ASolver_CG<AField<double,SIMD> >;

//====================================================================
// explicit instanciation for AField<float,SIMD>.
template<>
const std::string ASolver_CG<AField<float,SIMD> >::class_name
                              = "ASolver_CG<Afield<float,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<float,SIMD> >::register_factory();
}
#endif

template class ASolver_CG<AField<float,SIMD> >;

//============================================================END=====
