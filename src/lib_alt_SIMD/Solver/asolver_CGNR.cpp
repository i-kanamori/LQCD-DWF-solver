/*!
      @file    asolver_CGNR.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 1983 $
*/

#include "bridge_real.h"
#include "lib_alt/Solver/asolver_CGNR.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_SIMD/Field/afield.h"
#include "lib_alt_SIMD/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CGNR-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD>.
template<>
const std::string ASolver_CGNR<AField<double,SIMD> >::class_name
                              = "ASolver_CGNR<Afield<double,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<double,SIMD> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<double,SIMD> >;

//====================================================================
// explicit instanciation for AField<float,SIMD>.
template<>
const std::string ASolver_CGNR<AField<float,SIMD> >::class_name
                              = "ASolver_CGNR<Afield<float,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CGNR<AField<float,SIMD> >::register_factory();
}
#endif

template class ASolver_CGNR<AField<float,SIMD> >;

//============================================================END=====
