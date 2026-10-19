/*!
        @file    init_alt_SIMD.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2020-06-04 19:52:20 #$
        @version $LastChangedRevision: 2072 $
*/

#include "lib/IO/bridgeIO.h"
#include "lib_alt_SIMD/init_alt_SIMD.h"
#include "lib_alt_SIMD/bridge_init_factory_alt_SIMD.h"
using Bridge::vout;

bool init_alt_SIMD()
{
  vout.general("alt_SIMD code is being setup.\n");
#ifdef AVX512
  vout.general("    target instruction set: AVX-512\n");
#elif defined(AVX2)
  vout.general("    target instruction set: AVX2\n");
#elif defined(AVX1)
  vout.general("    target instruction set: AVX\n");
#else
  vout.crucial("    unkonw target for SIMD, aborting.\n");
  abort();
#endif
  bool result = true;

  // Factory initialization
#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
  result &= bridge_init_factory_alt_SIMD();
#endif

#ifdef DEBUG
  bridge_report_factory_alt_SIMD();
#endif

#endif /* USE_FACTORY */

  return result;
}

bool fin_alt_SIMD()
{
  return true;
}
