/*!

        @file    asolver_BiCGStab_Cmplx.cpp

        @brief   BiCGStab Solver

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $
*/
//====================================================================

#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_SIMD/Field/afield.h"

// function template files
#include "lib_alt_SIMD/Field/afield-inc.h"

// template functions
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx-tmpl.h"


//====================================================================
// explicit instanciation for AField<double>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<double,SIMD> >::class_name
                              = "ASolver_BiCGStab_Cmplx<Afield<double,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField<double,SIMD> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<double,SIMD> >;

//====================================================================
// explicit instanciation for AField<float>.
template<>
const std::string ASolver_BiCGStab_Cmplx<AField<float,SIMD> >::class_name
                              = "ASolver_BiCGStab_Cmplx<AField<float,SIMD> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = ASolver_BiCGStab_Cmplx<AField_dev<float,SIMD> >::register_factory();
}
#endif

template class ASolver_BiCGStab_Cmplx<AField<float,SIMD> >;

//============================================================END=====
