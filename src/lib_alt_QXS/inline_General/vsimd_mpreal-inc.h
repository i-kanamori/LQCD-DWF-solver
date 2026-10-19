/*!
        @file    vsimd_mpreal-inc.h
        @brief
        @author  Issaku Kanamori (kanamori)
        $LastChangedBy:  $
        @date    $LastChangedDate: $
        @version $LastChangedRevision: 2492 $
*/

#ifdef USE_MPREAL

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

typedef unsigned int uint_t;

typedef struct
{
  mpreal v[VLENS];
} Vsimd_t;

typedef struct
{
  int v[VLENS];
} Isimd_t;

typedef struct
{
  unsigned int v[VLENS];
} Usimd_t;

typedef struct
{
  bool v[VLENS];
} svbool_t;

typedef Vsimd_t svreal_t;

typedef Isimd_t svint_t;

typedef Usimd_t svuint_t;

typedef int     int_t;

#endif

#endif
