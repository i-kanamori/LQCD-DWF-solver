/*!
        @file    bridge_init_factory_alt_QXS.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2026-05-12 13:50:58 #$

        @version $LastChangedRevision: 2712 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_QXS_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_QXS_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_QXS();

#endif /* USE_FACTORY_AUTOREGISTER */

//#ifdef DEBUG
void bridge_report_factory_alt_QXS();
//#endif

bool bridge_fin_factory_alt_QXS();

#endif /* USE_FACTORY */


#endif
