/*!

        @file    $Id: site_list_util.cpp #$

        @brief   site index utility class

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision:  $

 */

#include <assert.h>
#include "lib/ResourceManager/threadManager.h"
#include "site_list_utils.h"

//====================================================================
namespace {
  inline
  void set_threadtask(int& ith, int& nth, int& is, int& ns,
                      const int size)
  {
    nth = ThreadManager_OpenMP::get_num_threads();
    ith = ThreadManager_OpenMP::get_thread_id();

    size_t tmp_is = (size_t)size * ith;
    size_t tmp_ns = (size_t)size * (ith+1);
    is = tmp_is/nth;
    ns = tmp_ns/nth;
  }

}

//====================================================================
void site_list_utils::init(const std::vector<int> &Rlattice, const mpi_bools &bools){

  m_bools=bools;

  m_Rx=Rlattice[DirX];
  m_Ry=Rlattice[DirY];
  m_Rz=Rlattice[DirZ];
  m_Rt=Rlattice[DirT];

  int Rxb=m_Rx;
  int Ryb=m_Ry;
  int Rzb=m_Rz;
  int Rtb=m_Rt;
  if(m_bools.has_mpi_x){
    Rxb=(Rxb>1)?(Rxb-2):0;
  }
  if(m_bools.has_mpi_y){
    Ryb=(Ryb>1)?(Ryb-2):0;
  }
  if(m_bools.has_mpi_z){
    Rzb=(Rzb>1)?(Rzb-2):0;
  }
  if(m_bools.has_mpi_t){
    Rtb=(Rtb>1)?(Rtb-2):0;
  }
  m_Rst=m_Rx*m_Ry*m_Rz*m_Rt;
  m_Rst_bulk=Rxb*Ryb*Rzb*Rtb;
  m_Rst_boundary=m_Rst-m_Rst_bulk;

#ifdef DEBUG
    vout.general(m_vl, "Rst_bulk=%d, Rst_boundary=%d\n", m_Rst_bulk, m_Rst_boundary);
    vout.general(m_vl, "  Rxb,Ryb,Rzb,Rtb=%d,%d,%d,%d\n", Rxb, Ryb, Rzb, Rtb);
    vout.general(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
#endif
}

/*******************************************************************

  lexical case

********************************************************************/

//====================================================================
void site_list_utils::set_site_lex_bulk(vector_site_t &site_list_bulk) const {
  ThreadManager_OpenMP::assert_single_thread("set_site_list_lex");

  // obtain reduced volumes (total in the node, bulk, boundary)
  //  int Rst, Rst_bulk, Rst_boundary;
  //  calc_rvol(Rst, Rst_bulk, Rst_boundary);

  if(m_Rst>0){
    site_list_bulk.resize(m_Rst);
  }
#pragma omp parallel
  {
    int ith, nth, is_bulk, ns_bulk;

    // local copy in the thread
    const mpi_bools bools(m_bools);
    const int Rx=m_Rx;
    const int Ry=m_Ry;
    const int Rz=m_Rz;
    const int Rt=m_Rt;

    set_threadtask(ith, nth, is_bulk, ns_bulk, m_Rst);

    // bulk
    for(int s0=is_bulk; s0<ns_bulk; ++s0){

      // get coordinates
      int idx=s0;
      int rx = idx % Rx;
      idx = idx / Rx;
      int ry = idx % Ry;
      idx = idx / Ry;
      int rz = idx % Rz;
      int rt = idx / Rz;

      const int s=rx+Rx*(ry+Ry*(rz+Rz*rt));
      const int drx_p=(rx+1) % Rx -rx;
      const int drx_m=(rx-1 + Rx) % Rx -rx;
      const int dry_p=(ry+1) % Ry -ry;
      const int dry_m=(ry-1 + Ry) % Ry -ry;
      const int drz_p=(rz+1) % Rz -rz;
      const int drz_m=(rz-1 + Rz) % Rz -rz;
      const int drt_p=(rt+1) % Rt -rt;
      const int drt_m=(rt-1 + Rt) % Rt -rt;

      int sn_xp = s + drx_p;
      int sn_xm = s + drx_m;
      int sn_yp = s + Rx*dry_p;
      int sn_ym = s + Rx*dry_m;
      int sn_zp = s + Rx*Ry*drz_p;
      int sn_zm = s + Rx*Ry*drz_m;
      int sn_tp = s + Rx*Ry*Rz*drt_p;
      int sn_tm = s + Rx*Ry*Rz*drt_m;

#ifdef USE_MPI
      // boundary correctoins
      if(bools.has_mpi_x){
        if(rx==0){    sn_xm=-1; }
        if(rx==Rx-1){ sn_xp=-1; }
      }
      if(bools.has_mpi_y){
        if(ry==0){    sn_ym=-1; }
        if(ry==Ry-1){ sn_yp=-1; }
      }
      if(bools.has_mpi_z){
        if(rz==0){    sn_zm=-1; }
        if(rz==Rz-1){ sn_zp=-1; }
      }
      if(bools.has_mpi_t){
        if(rt==0){    sn_tm=-1; }
        if(rt==Rt-1){ sn_tp=-1; }
      }
#endif // USE_MPI
      site_list_bulk[s0].here=s;
      site_list_bulk[s0].packed[DirX]=ry + Ry*(rz + rt*Rz);
      site_list_bulk[s0].packed[DirY]=rx + Rx*(rz + rt*Rz);
      site_list_bulk[s0].packed[DirZ]=rx + Rx*(ry + rt*Ry);
      site_list_bulk[s0].packed[DirT]=rx + Rx*(ry + rz*Ry);
      site_list_bulk[s0].neighbor[DirXp]=sn_xp;
      site_list_bulk[s0].neighbor[DirXm]=sn_xm;
      site_list_bulk[s0].neighbor[DirYp]=sn_yp;
      site_list_bulk[s0].neighbor[DirYm]=sn_ym;
      site_list_bulk[s0].neighbor[DirZp]=sn_zp;
      site_list_bulk[s0].neighbor[DirZm]=sn_zm;
      site_list_bulk[s0].neighbor[DirTp]=sn_tp;
      site_list_bulk[s0].neighbor[DirTm]=sn_tm;
    }// s0
  } // parallel
}

//====================================================================
void site_list_utils::set_site_lex_boundary(vector_site_t &site_list_boundary,
                                            vector_pack_t &pack_list_x,
                                            vector_pack_t &pack_list_y,
                                            vector_pack_t &pack_list_z,
                                            vector_pack_t &pack_list_t,
                                            const vector_site_t &site_list_bulk) const {

#ifdef USE_MPI

  // obtain reduced volumes (total in the node, bulk, boundary)
  //  int Rst, Rst_bulk, Rst_boundary;
  //  calc_rvol(Rst, Rst_bulk, Rst_boundary);

  // a list of the boundary sites
  aligned_vector<size_t> boundary_idx(m_Rst_boundary);
  size_t s=0;
  for(size_t s0=0; s0<m_Rst; ++s0){
    const long_t *sn=site_list_bulk[s0].neighbor;
    bool flagx=( (sn[DirXp]+1)*(sn[DirXm]+1) == 0);
    bool flagy=( (sn[DirYp]+1)*(sn[DirYm]+1) == 0);
    bool flagz=( (sn[DirZp]+1)*(sn[DirZm]+1) == 0);
    bool flagt=( (sn[DirTp]+1)*(sn[DirTm]+1) == 0);
    if(flagx || flagy || flagz || flagt ){
      assert(s<m_Rst_boundary);
      boundary_idx[s]=site_list_bulk[s0].here;
      s++;
    }
  }
  assert(s==m_Rst_boundary);


  // allocate list vector for bounary sites
  site_list_boundary.resize(m_Rst_boundary);

  // allocate list vector for bounary packing
  size_t Nbx0=0;
  size_t Nby0=0;
  size_t Nbz0=0;
  size_t Nbt0=0;
  if(m_bools.has_mpi_x){
    Nbx0=(size_t)m_Ry*m_Rz*m_Rt;
    pack_list_x.resize(Nbx0);
  }
  if(m_bools.has_mpi_y){
    Nby0=(size_t)m_Rx*m_Rz*m_Rt;
    pack_list_y.resize(Nby0);
  }
  if(m_bools.has_mpi_z){
    Nbz0=(size_t)m_Rx*m_Ry*m_Rt;
    pack_list_z.resize(Nbz0);
  }
  if(m_bools.has_mpi_t){
    Nbt0=(size_t)m_Rx*m_Ry*m_Rz;
    pack_list_t.resize(Nbt0);
  }

  // set list vector for the bounary sites
#pragma omp parallel
  {
    int ith, nth, is_boundary, ns_boundary;
    int isx, nsx, isy, nsy, isz, nsz, ist, nst;

    // local copy in the thread
    const int Rx=m_Rx;
    const int Ry=m_Ry;
    const int Rz=m_Rz;
    const int Rt=m_Rt;

    const bool has_mpi_x=m_bools.has_mpi_x;
    const bool has_mpi_y=m_bools.has_mpi_y;
    const bool has_mpi_z=m_bools.has_mpi_z;
    const bool has_mpi_t=m_bools.has_mpi_t;

    const int dist0=m_dist0;
    const int dist1=m_dist1;

    set_threadtask(ith, nth, is_boundary, ns_boundary, m_Rst_boundary);
    set_threadtask(ith, nth, isx, nsx, Nbx0);
    set_threadtask(ith, nth, isy, nsy, Nby0);
    set_threadtask(ith, nth, isz, nsz, Nbz0);
    set_threadtask(ith, nth, ist, nst, Nbt0);


    for(int s0=is_boundary; s0<ns_boundary; ++s0){
      // get coordinates
      const int s=boundary_idx[s0];
      int idx=s;
      int rx = idx % Rx;
      idx = idx / Rx;
      int ry = idx % Ry;
      idx = idx / Ry;
      int rz = idx % Rz;
      int rt = idx / Rz;

      int sn_xp=-1;
      int sn_xm=-1;
      int sn_yp=-1;
      int sn_ym=-1;
      int sn_zp=-1;
      int sn_zm=-1;
      int sn_tp=-1;
      int sn_tm=-1;

      if(has_mpi_x){
        if(rx==0){    sn_xm=s+(Rx-1); }
        if(rx==Rx-1){ sn_xp=s-(Rx-1); }
      }
      if(has_mpi_y){
        if(ry==0){    sn_ym=s+Rx*(Ry-1); }
        if(ry==Ry-1){ sn_yp=s-Rx*(Ry-1); }
      }
      if(has_mpi_z){
        if(rz==0){    sn_zm=s+Rx*Ry*(Rz-1); }
        if(rz==Rz-1){ sn_zp=s-Rx*Ry*(Rz-1); }
      }
      if(has_mpi_t){
        if(rt==0){    sn_tm=s+Rx*Ry*Rz*(Rt-1); }
        if(rt==Rt-1){ sn_tp=s-Rx*Ry*Rz*(Rt-1); }
      }

      site_list_boundary[s0].here=s;
      site_list_boundary[s0].packed[DirX]=ry + Ry*(rz + Rz*rt);
      site_list_boundary[s0].packed[DirY]=rx + Rx*(rz + Rz*rt);
      site_list_boundary[s0].packed[DirZ]=rx + Rx*(ry + Ry*rt);
      site_list_boundary[s0].packed[DirT]=rx + Rx*(ry + Ry*rz);
      site_list_boundary[s0].neighbor[DirXp]=sn_xp;
      site_list_boundary[s0].neighbor[DirXm]=sn_xm;
      site_list_boundary[s0].neighbor[DirYp]=sn_yp;
      site_list_boundary[s0].neighbor[DirYm]=sn_ym;
      site_list_boundary[s0].neighbor[DirZp]=sn_zp;
      site_list_boundary[s0].neighbor[DirZm]=sn_zm;
      site_list_boundary[s0].neighbor[DirTp]=sn_tp;
      site_list_boundary[s0].neighbor[DirTm]=sn_tm;
    } //s0 for boundary

    { // x-direction  (for packing)
      //  packed = 0        ry +    Ry*rz +    Ry*Rz*rt
      //      ss = 0   + Rx*ry + Rx*Ry*rz + Rx*Ry*Rz*rt
      for(int i=isx; i<nsx; ++i){
        const int s_packed=i;
        const int ss=s_packed*Rx;
        const int s0=((i+dist0)*Rx);
        const int s1=((i+dist1)*Rx);
        pack_list_x[i].here=ss;
        pack_list_x[i].neighbor=-1;
        pack_list_x[i].s0=s0;
        pack_list_x[i].s1=s1;
      }
    } // x-direction

    { // y-direction  (for packing)
      //  packed = rx + 0 +    Rx*(rz + Rz*rt)
      //      ss = rx + 0 + Ry*Rx*(rz + Rz*rt)
      for(int i=isy; i<nsy; ++i){
        const int &s_packed=i;
        const int tmp=s_packed/Rx;  // rz + Rz*rt
        const int tmp0=(i+dist0)/Rx;
        const int tmp1=(i+dist1)/Rx;
        const int ss=s_packed + tmp*Rx*(Ry-1);
        const int s0=i+dist0 + tmp0*Rx*(Ry-1);
        const int s1=i+dist1 + tmp1*Rx*(Ry-1);
        pack_list_y[i].here=ss;
        pack_list_y[i].neighbor=-1;
        pack_list_y[i].s0=s0;
        pack_list_y[i].s1=s1;
      }
    } // y-direction

    { // z-direction  (for packing)
      //  packed = rx + Rx*ry  + 0 + Rx*Ry*rt
      //      ss = rx + Rx*ry  + 0 + Rx*Ry*Rz*rt
      for(int i=isz; i<nsz; ++i){
        const int &s_packed=i;
        int tmpq=s_packed/(Rx*Ry);  // rt
        int tmpr=s_packed%(Rx*Ry);  // rx + Rx*ry
        int tmpq0=(i+dist0)/(Rx*Ry);
        int tmpr0=(i+dist0)%(Rx*Ry);
        int tmpq1=(i+dist1)/(Rx*Ry);
        int tmpr1=(i+dist1)%(Rx*Ry);
        const int ss= tmpr  + Rx*Ry*Rz*tmpq;
        const int s0= tmpr0 + Rx*Ry*Rz*tmpq0;
        const int s1= tmpr1 + Rx*Ry*Rz*tmpq1;
        pack_list_z[i].here=ss;
        pack_list_z[i].neighbor=-1;
        pack_list_z[i].s0=s0;
        pack_list_z[i].s1=s1;
      }
    } // z-direction

    { // t-direction  (for packing)
      //  packed = rx + Rx*ry  + Rx*Ry*rz
      //      ss = rx + Rx*ry  + Rx*Ry*rz
      for(int i=ist; i<nst; ++i){
        const int &s_packed=i;
        const int ss=s_packed;
        const int s0= (i+dist0);
        const int s1= (i+dist1);
        pack_list_t[i].here=ss;
        pack_list_t[i].neighbor=-1;
        pack_list_t[i].s0=s0;
        pack_list_t[i].s1=s1;
      }
    } // t-direction

  } // omp parallel
#endif // USE_MPI

}


/*******************************************************************

  even odd case

********************************************************************/

//====================================================================
void site_list_utils::set_site_eo_bulk(vector_site_t *site_list_bulk, const int eo_offset_) const {
  ThreadManager_OpenMP::assert_single_thread("set_site_eo_bulk");

  // obtain reduced volumes (total in the node, bulk, boundary)
  //  int Rst, Rst_bulk, Rst_boundary;
  //  calc_rvol(Rst, Rst_bulk, Rst_boundary);

  assert(m_Rx%2==0);
  if(m_Rx % 2 != 0){
    vout.crucial(m_vl, "%s: Rx must be even for even-odd site list.\n", class_name.c_str());
    vout.crucial(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
    exit(EXIT_FAILURE);
  }
  if(m_Rx < 2 ){
    vout.crucial("%s: too small Rx for even-odd site list.\n", class_name.c_str());
    vout.crucial(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
    exit(EXIT_FAILURE);
  }
  if(m_Ry < 1 ){
    vout.crucial("%s: too small Ry for even-odd site list.\n", class_name.c_str());
    vout.crucial(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
    exit(EXIT_FAILURE);
  }
  if(m_Rz < 1 ){
    vout.crucial("%s: too small Rz for even-odd site list.\n", class_name.c_str());
    vout.crucial(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
    exit(EXIT_FAILURE);
  }
  if(m_Rt < 1 ){
    vout.crucial("%s: too small Rt for even-odd site list.\n", class_name.c_str());
    vout.crucial(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx, m_Ry, m_Rz, m_Rt);
    exit(EXIT_FAILURE);
  }

  const int Rst2=m_Rst/2;
  assert(Rst2>0);
  site_list_bulk[0].resize(Rst2);
  site_list_bulk[1].resize(Rst2);

#pragma omp parallel
  {
    int ith, nth, is_bulk, ns_bulk;

    const mpi_bools bools(m_bools);
    const int eo_offset=eo_offset_;
    const int Rx=m_Rx;
    const int Rx2=m_Rx/2;
    const int Ry=m_Ry;
    const int Rz=m_Rz;
    const int Rt=m_Rt;

    set_threadtask(ith, nth, is_bulk, ns_bulk, Rst2);

    for(int ieo=0; ieo<2; ieo++){
      int ieo_total=(ieo+eo_offset)%2;

      // bulk
      site_list_t *list_bulk=&site_list_bulk[ieo][0];
      for(int s0=is_bulk; s0<ns_bulk; ++s0){
        // get coordinates
        int idx=s0;
        int rx2 = idx % Rx2;
        idx = idx / Rx2;
        int ry = idx % Ry;
        idx = idx / Ry;
        int rz = idx % Rz;
        int rt = idx / Rz;
        int ieo_work=(ieo_total+ry+rz+rt)%2;
        // int rx=2*rx2+ieo_work;

        const int s=rx2+Rx2*(ry+Ry*(rz+Rz*rt));
        const int drx_p=(rx2+1) % Rx2 -rx2;
        const int drx_m=(rx2-1 + Rx2) % Rx2 -rx2;
        const int dry_p=(ry+1) % Ry -ry;
        const int dry_m=(ry-1 + Ry) % Ry -ry;
        const int drz_p=(rz+1) % Rz -rz;
        const int drz_m=(rz-1 + Rz) % Rz -rz;
        const int drt_p=(rt+1) % Rt -rt;
        const int drt_m=(rt-1 + Rt) % Rt -rt;


        int sn_xp = s + ieo_work*drx_p;
        int sn_xm = s + (1-ieo_work)*drx_m;
        int sn_yp = s + Rx2*dry_p;
        int sn_ym = s + Rx2*dry_m;
        int sn_zp = s + Rx2*Ry*drz_p;
        int sn_zm = s + Rx2*Ry*drz_m;
        int sn_tp = s + Rx2*Ry*Rz*drt_p;
        int sn_tm = s + Rx2*Ry*Rz*drt_m;

#ifdef USE_MPI
        if(bools.has_mpi_x){   // boundary correctoins: need to use rx and Rx
          if(rx2==0){     sn_xm=-1; }
          if(rx2==Rx2-1){ sn_xp=-1; }
        }
        if(bools.has_mpi_y){
          if(ry==0){    sn_ym=-1; }
          if(ry==Ry-1){ sn_yp=-1; }
        }
        if(bools.has_mpi_z){
          if(rz==0){    sn_zm=-1; }
          if(rz==Rz-1){ sn_zp=-1; }
        }
        if(bools.has_mpi_t){
          if(rt==0){    sn_tm=-1; }
          if(rt==Rt-1){ sn_tp=-1; }
        }

#endif // USE_MPI

        list_bulk[s0].here=s;
        list_bulk[s0].ieo_work=ieo_work;
        list_bulk[s0].packed[DirX]=ry  + Ry *(rz + Rz*rt);
        list_bulk[s0].packed[DirY]=rx2 + Rx2*(rz + Rz*rt);
        list_bulk[s0].packed[DirZ]=rx2 + Rx2*(ry + Ry*rt);
        list_bulk[s0].packed[DirT]=rx2 + Rx2*(ry + Ry*rz);
        list_bulk[s0].neighbor[DirXp]=sn_xp;
        list_bulk[s0].neighbor[DirXm]=sn_xm;
        list_bulk[s0].neighbor[DirYp]=sn_yp;
        list_bulk[s0].neighbor[DirYm]=sn_ym;
        list_bulk[s0].neighbor[DirZp]=sn_zp;
        list_bulk[s0].neighbor[DirZm]=sn_zm;
        list_bulk[s0].neighbor[DirTp]=sn_tp;
        list_bulk[s0].neighbor[DirTm]=sn_tm;
      }// s0
    } //ieo
  } // parallel

#ifdef DEBUG
  vout.general(m_vl, "site_list_utils::set_site_bulk, done\n");
#endif
}

//====================================================================
void site_list_utils::set_site_eo_boundary(vector_site_t *site_list_boundary,
                                           vector_pack_t *pack_list_x,
                                           vector_pack_t *pack_list_y,
                                           vector_pack_t *pack_list_z,
                                           vector_pack_t *pack_list_t,
                                           const vector_site_t *site_list_bulk,
                                           const int eo_offset_) const {

#ifdef USE_MPI

  int Rxb=m_Rx;
  int Ryb=m_Ry;
  int Rzb=m_Rz;
  int Rtb=m_Rt;
  if(m_bools.has_mpi_x){
    Rxb=(Rxb>1)?(Rxb-2):0;
  }
  if(m_bools.has_mpi_y){
    Ryb=(Ryb>1)?(Ryb-2):0;
  }
  if(m_bools.has_mpi_z){
    Rzb=(Rzb>1)?(Rzb-2):0;
  }
  if(m_bools.has_mpi_t){
    Rtb=(Rtb>1)?(Rtb-2):0;
  }
  const int Rst=m_Rx*m_Ry*m_Rz*m_Rt;
  const int Rst2=Rst/2;
  const int Rst_bulk=Rxb*Ryb*Rzb*Rtb;
  const int Rst_boundary=Rst-Rst_bulk;
  const int Rst_boundary2=Rst_boundary/2;

#ifdef DEBUG
  vout.general(m_vl, "Rst_bulk=%d, Rst_boundary=%d\n", Rst_bulk, Rst_boundary);
  vout.general(m_vl, "  Rxb,Ryb,Rzb,Rtb=%d,%d,%d,%d: bulk=%d\n", Rxb,Ryb,Rzb,Rtb, Rxb*Ryb*Rzb*Rtb);
  vout.general(m_vl, "  Rx ,Ry ,Rz ,Rt =%d,%d,%d,%d\n", m_Rx,m_Ry,m_Rz,m_Rt);
  fflush(0);
#endif

  if(Rst_boundary==0) { return; }
  assert(Rst2>0);


  // make a temporal list of boundary sites
  aligned_vector<size_t> boundary_idx[2];
  boundary_idx[0].resize(Rst_boundary2);
  boundary_idx[1].resize(Rst_boundary2);

#pragma omp parallel for
  for(int ieo=0; ieo<2; ieo++){
    size_t *this_boundary_idx=&boundary_idx[ieo][0];
    const site_list_t *list_bulk=&site_list_bulk[ieo][0];
    size_t s=0;
    for(size_t s0=0; s0<Rst2; s0++){
      const long_t *sn=list_bulk[s0].neighbor;
      bool flagx=( (sn[DirXp]+1)*(sn[DirXm]+1) == 0);
      bool flagy=( (sn[DirYp]+1)*(sn[DirYm]+1) == 0);
      bool flagz=( (sn[DirZp]+1)*(sn[DirZm]+1) == 0);
      bool flagt=( (sn[DirTp]+1)*(sn[DirTm]+1) == 0);
      if(flagx || flagy || flagz || flagt ){
        assert(s<Rst_boundary2);
        this_boundary_idx[s]=list_bulk[s0].here;
        s++;
      }
    }
  }// ieo, omp parallel for

  // allocate list vector for bounary sites
  if(Rst_boundary2>0){
    site_list_boundary[0].resize(Rst_boundary2);
    site_list_boundary[1].resize(Rst_boundary2);
  }

  // allocate list vector for bounary packing
  size_t Nbx0=0;
  size_t Nby0=0;
  size_t Nbz0=0;
  size_t Nbt0=0;

  if(m_bools.has_mpi_x){ Nbx0=(size_t)m_Ry*m_Rz*m_Rt; }
  if(m_bools.has_mpi_y){ Nby0=(size_t)m_Rx*m_Rz*m_Rt/2; }
  if(m_bools.has_mpi_z){ Nbz0=(size_t)m_Rx*m_Ry*m_Rt/2; }
  if(m_bools.has_mpi_t){ Nbt0=(size_t)m_Rx*m_Ry*m_Rz/2; }
  if(Nbx0>0){
    pack_list_x[0].resize(Nbx0);
    pack_list_x[1].resize(Nbx0);
  }
  if(Nby0>0){
    pack_list_y[0].resize(Nby0);
    pack_list_y[1].resize(Nby0);
  }
  if(Nbz0>0){
    pack_list_z[0].resize(Nbz0);
    pack_list_z[1].resize(Nbz0);
  }
  if(Nbt0>0){
    pack_list_t[0].resize(Nbt0);
    pack_list_t[1].resize(Nbt0);
  }
  // set list vector for the bounary sites

#pragma omp parallel
  {
    int ith, nth, is_boundary, ns_boundary;
    int isx, nsx, isy, nsy, isz, nsz, ist, nst;

    // local copy
    const mpi_bools bools(m_bools);
    const int eo_offset=eo_offset_;
    const int Rx=m_Rx;
    const int Rx2=m_Rx/2;
    const int Ry=m_Ry;
    const int Rz=m_Rz;
    const int Rt=m_Rt;
    const bool has_mpi_x=m_bools.has_mpi_x;
    const bool has_mpi_y=m_bools.has_mpi_y;
    const bool has_mpi_z=m_bools.has_mpi_z;
    const bool has_mpi_t=m_bools.has_mpi_t;
    const int dist0=m_dist0;
    const int dist1=m_dist1;

    set_threadtask(ith, nth, is_boundary, ns_boundary, Rst_boundary2);
    set_threadtask(ith, nth, isx, nsx, Nbx0);
    set_threadtask(ith, nth, isy, nsy, Nby0);
    set_threadtask(ith, nth, isz, nsz, Nbz0);
    set_threadtask(ith, nth, ist, nst, Nbt0);

    // boundary sites
    for(int ieo=0; ieo<2; ieo++){
      site_list_t *list_boundary = &site_list_boundary[ieo][0];
      const size_t *this_boundary_idx = &boundary_idx[ieo][0];
      int ieo_total=(ieo+eo_offset)%2;
      for(int s0=is_boundary; s0<ns_boundary; ++s0){
        // get coordinates
        const int s=this_boundary_idx[s0];
        int idx=s;
        const int rx2 = idx % Rx2;
        idx = idx / Rx2;
        const int ry = idx % Ry;
        idx = idx / Ry;
        const int rz = idx % Rz;
        const int rt = idx / Rz;
        assert(s==rx2+Rx2*(ry+Ry*(rz+Rz*rt)));
        int ieo_work=(ieo_total+ry+rz+rt)%2;

        // boundary correctoins
        int sn_xp=-1;
        int sn_xm=-1;
        int sn_yp=-1;
        int sn_ym=-1;
        int sn_zp=-1;
        int sn_zm=-1;
        int sn_tp=-1;
        int sn_tm=-1;

        if(has_mpi_x){
          if(rx2==0){    sn_xm=s+(1-ieo_work)*(Rx2-1); }
          if(rx2==Rx2-1){ sn_xp=s-ieo_work*(Rx2-1); }
        }
        if(has_mpi_y){
          if(ry==0){    sn_ym=s+Rx2*(Ry-1); }
          if(ry==Ry-1){ sn_yp=s-Rx2*(Ry-1); }
        }
        if(has_mpi_z){
          if(rz==0){    sn_zm=s+Rx2*Ry*(Rz-1); }
          if(rz==Rz-1){ sn_zp=s-Rx2*Ry*(Rz-1); }
        }
        if(has_mpi_t){
          if(rt==0){    sn_tm=s+Rx2*Ry*Rz*(Rt-1); }
          if(rt==Rt-1){ sn_tp=s-Rx2*Ry*Rz*(Rt-1); }
        }

        list_boundary[s0].here=s;
        list_boundary[s0].packed[DirX]=ry  + Ry *(rz + Rz*rt);
        list_boundary[s0].packed[DirY]=rx2 + Rx2*(rz + Rz*rt);
        list_boundary[s0].packed[DirZ]=rx2 + Rx2*(ry + Ry*rt);
        list_boundary[s0].packed[DirT]=rx2 + Rx2*(ry + Ry*rz);
        list_boundary[s0].neighbor[DirXp]=sn_xp;
        list_boundary[s0].neighbor[DirXm]=sn_xm;
        list_boundary[s0].neighbor[DirYp]=sn_yp;
        list_boundary[s0].neighbor[DirYm]=sn_ym;
        list_boundary[s0].neighbor[DirZp]=sn_zp;
        list_boundary[s0].neighbor[DirZm]=sn_zm;
        list_boundary[s0].neighbor[DirTp]=sn_tp;
        list_boundary[s0].neighbor[DirTm]=sn_tm;
      } //s0 for boundary

      if(Nbx0>0){ // x-direction  (for packing)
        //  packed = 0         ry +     Ry*rz +     Ry*Rz*rt
        //      ss = 0   + Rx2*ry + Rx2*Ry*rz + Rx2*Ry*Rz*rt
        site_list_pack_t *list_pack_x = &pack_list_x[ieo][0];
        for(int i=isx; i<nsx; ++i){
          const int s_packed=i;
          const int ss=s_packed*Rx2;
          const int s0=((i+dist0)*Rx2);
          const int s1=((i+dist1)*Rx2);
          list_pack_x[i].here=ss;
          list_pack_x[i].neighbor=-1;
          list_pack_x[i].s0=s0;
          list_pack_x[i].s1=s1;
        }
      } // x-direction

      if(Nby0>0){ // y-direction  (for packing)
        //  packed = rx2 + 0 +    Rx2*(rz + Rz*rt)
        //      ss = rx2 + 0 + Ry*Rx2*(rz + Rz*rt)
        site_list_pack_t *list_pack_y = &pack_list_y[ieo][0];
        for(int i=isy; i<nsy; ++i){
          const int &s_packed=i;
          const int tmp=s_packed/Rx2;  // rz + Rz*rt
          const int tmp0=(i+dist0)/Rx2;
          const int tmp1=(i+dist1)/Rx2;
          const int ss=s_packed + tmp*Rx2*(Ry-1);
          const int s0=i+dist0 + tmp0*Rx2*(Ry-1);
          const int s1=i+dist1 + tmp1*Rx2*(Ry-1);
          list_pack_y[i].here=ss;
          list_pack_y[i].neighbor=-1;
          list_pack_y[i].s0=s0;
          list_pack_y[i].s1=s1;
        }
      } // y-direction

      if(Nbz0>0){ // z-direction  (for packing)
        //  packed = rx2 + Rx2*ry  + 0 + Rx2*Ry*rt
        //      ss = rx2 + Rx2*ry  + 0 + Rx2*Ry*Rz*rt
        site_list_pack_t *list_pack_z = &pack_list_z[ieo][0];
        for(int i=isz; i<nsz; ++i){
          const int &s_packed=i;
          int tmpq=s_packed/(Rx2*Ry);  // rt
          int tmpr=s_packed%(Rx2*Ry);  // rx2 + Rx2*ry
          int tmpq0=(i+dist0)/(Rx2*Ry);
          int tmpr0=(i+dist0)%(Rx2*Ry);
          int tmpq1=(i+dist1)/(Rx2*Ry);
          int tmpr1=(i+dist1)%(Rx2*Ry);
          const int ss= tmpr  + Rx2*Ry*Rz*tmpq;
          const int s0= tmpr0 + Rx2*Ry*Rz*tmpq0;
          const int s1= tmpr1 + Rx2*Ry*Rz*tmpq1;
          list_pack_z[i].here=ss;
          list_pack_z[i].neighbor=-1;
          list_pack_z[i].s0=s0;
          list_pack_z[i].s1=s1;
        }
      } // z-direction

      if(Nbt0>0){ // t-direction  (for packing)
        //  packed = rx2 + Rx2*ry  + Rx2*Ry*rz
        //      ss = rx2 + Rx2*ry  + Rx2*Ry*rz
        site_list_pack_t *list_pack_t = &pack_list_t[ieo][0];
        for(int i=ist; i<nst; ++i){
          const int &s_packed=i;
          const int ss=s_packed;
          const int s0= (i+dist0);
          const int s1= (i+dist1);
          list_pack_t[i].here=ss;
          list_pack_t[i].neighbor=-1;
          list_pack_t[i].s0=s0;
          list_pack_t[i].s1=s1;
        }
      } // t-direction

    } // ieo
  } // omp parallel



#endif // USE_MPI
#ifdef DEBUG
  vout.general(m_vl, "site_list_utils::set_site_eo_boundary, done\n");
#endif
}

