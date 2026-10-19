/*!
      @file    asolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "bridge_real.h"
#include "bridge_mpreal.h"

#include "lib_alt/Solver/asolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Solver/asolver_CG-tmpl.h"


//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string ASolver_CG<AField<double, QXS> >::class_name
  = "ASolver_CG<Afield<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<double, QXS> >::register_factory();
}
#endif

template class ASolver_CG<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string ASolver_CG<AField<float, QXS> >::class_name
  = "ASolver_CG<Afield<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<float, QXS> >::register_factory();
}
#endif

template class ASolver_CG<AField<float, QXS> >;


#ifdef USE_MPREAL
//====================================================================
// explicit instanciation for AField<mpreal,QXS>.
using mpfr::mpreal;
template<>
const std::string ASolver_CG<AField<mpreal, QXS> >::class_name
  = "ASolver_CG<Afield<mpreal,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_CG<AField<mpreal, QXS> >::register_factory();
}
#endif

template class ASolver_CG<AField<mpreal, QXS> >;
#endif

//============================================================END=====
