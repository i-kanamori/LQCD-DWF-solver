/*!

        @file    $Id: alt_simd_def.h #$

        @brief    header to lib_alt_Simd

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $
 */
////////////////////////////////////////////////////


#ifndef ALT_SIMD_DEF_INCLUDED
#define ALT_SIMD_DEF_INCLUDED


// data layout implementation type
#define IMPLTYPE_ALT_IMPL1
//#define IMPLTYPE_ALT_GRIDSIMD

#include "inline/valtype.h"
#include "inline/simd_layout.h"

#define EXCLUSIVE_THREAD_FOR_COMM
//#define FORCE_COMMUNICATE
//#define AFIELD_FIRST_TOUCH
//#define AFOPR_FIRST_TOUCH_BUF

#endif // include guard
