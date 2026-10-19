#ifndef SITE_LIST_H
#define SITE_LIST_H
// Time-stamp: <2020-06-02 00:50:08 kanamori>
/*!

        @file    $Id: site_list.h #$

        @brief   site index

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */

namespace {
  constexpr int DirX=0;
  constexpr int DirY=1;
  constexpr int DirZ=2;
  constexpr int DirT=3;

  constexpr int DirXp=0;
  constexpr int DirYp=2;
  constexpr int DirZp=4;
  constexpr int DirTp=6;

  constexpr int DirXm=1;
  constexpr int DirYm=3;
  constexpr int DirZm=5;
  constexpr int DirTm=7;
}

// note:
//  ieo_work = (eo_offset+ eo + ry+rz+rt)%2
//           = (eo_offset + eo_x) %2
struct site_list_t{
  long_t here;
  long_t ieo_work, dummy1, dummy2;
  //  long_t coord[4];
  long_t packed[4];
  long_t neighbor[8];
};

struct site_list_pack_t{
  long_t here;
  long_t neighbor;// neighboring site
  long_t s0;      // for prefetch 0
  long_t s1;      // for prefetch 1
};

#endif
