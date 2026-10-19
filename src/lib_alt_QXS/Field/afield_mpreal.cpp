/*!
        @file    afield_float.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
        @version $LastChangedRevision: 2492 $
*/

#ifdef USE_MPREAL

#include "lib_alt_QXS/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

#include "bridge_mpreal.h"

using mpfr::mpreal;

typedef mpreal real_t;

#include "lib_alt_QXS/inline/vsimd_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_mpreal-inc.h"


// template definition
#include "lib_alt_QXS/Field/afield-tmpl.h"


template<>
const std::string AField<mpreal, QXS>::class_name = "AField<mpreal, QXS>";


// explicit instanciation.
template class AField<mpreal, QXS>;

#endif
//============================================================END=====
