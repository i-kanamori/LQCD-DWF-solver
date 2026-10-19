/*!
      @file    shiftAField_lex_mpreal.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-02-28 16:09:41 #$
      @version $LastChangedRevision: 2492 $
*/

#ifdef USE_MPREAL

#include "bridge_mpreal.h"
#include "lib_alt_QXS/Field/shiftAField_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

using mpfr::mpreal;
typedef mpreal real_t;

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"

// inline function files
#include "lib_alt_QXS/inline/vsimd_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_mpreal-inc.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/shiftAField_lex-tmpl.h"

template<>
const std::string ShiftAField_lex<AField<mpreal, QXS> >::class_name
  = "ShiftAField_lex<AField<mpreal,QXS> >";

// explicit instanciation.
template class ShiftAField_lex<AField<mpreal, QXS> >;

#endif
//============================================================END=====
