/*!
        @file    fprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2019-07-30 20:13:19 #$
        @version $LastChangedRevision: 1983 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_SIMD/Field/afield.h"


// explicit instanciation for AField<double>.

template class Fprop_alt<AField<double,SIMD> >;

template class Fprop_alt<AField<float,SIMD> >;

//============================================================END=====
