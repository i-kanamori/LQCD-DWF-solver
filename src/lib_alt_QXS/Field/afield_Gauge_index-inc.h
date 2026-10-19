/*!
        @file    afield_Gauge_index-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-27 12:50:00 #$
        @version $LastChangedRevision: 2695 $
*/

#ifndef QXS_AFIELD_GAUGE_INDEX_INC_INCLUDED
#define QXS_AFIELD_GAUGE_INDEX_INC_INCLUDED

namespace {
  inline int idxGr(int ic1, int ic2){ return     2*(ic1 + NC * ic2); }
  inline int idxGi(int ic1, int ic2){ return 1 + 2*(ic1 + NC * ic2); }
}

#endif
