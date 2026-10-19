/*!
        @file    bridge_init_factory_alt_SIMD.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2019-07-29 21:19:48 #$

        @version $LastChangedRevision: 1982 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_SIMD_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_SIMD_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_SIMD();

#endif /* USE_FACTORY_AUTOREGISTER */

#ifdef DEBUG
void bridge_report_factory_alt_SIMD();
#endif

#endif /* USE_FACTORY */


#endif
