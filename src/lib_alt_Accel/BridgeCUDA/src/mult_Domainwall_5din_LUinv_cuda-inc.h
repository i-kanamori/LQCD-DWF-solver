/*!
      @file    mult_Doainwall_5din_LUinv_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2710 $
*/

#ifndef MULT_DOMAINWALL_5DIN_LUINV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_LUINV_ACC_INCLUDED

//====================================================================
__global__ void mult_domainwall_5din_LUinv_dirac_kernel(
                        real_t* __restrict__ vp, real_t* __restrict__ wp,
                        int Ns, int Nx, int Ny, int Nz, int Nt,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  const int Nst  = Nx * Ny * Nz * Nt;
  const int Nin5 = NVCD * Ns;
  const int idx  = blockIdx.x * blockDim.x + threadIdx.x;

  if ( idx  < Nst * NVC) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;
    //anchor
    __shared__ real_t vt[ND], yt[ND], xt[ND];

    int is = 0;
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = e[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
	      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += e[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    vt[0] += a * (xt[0] + xt[2]);
    vt[1] += a * (xt[1] + xt[3]);
    vt[2] += a * (xt[2] + xt[0]);
    vt[3] += a * (xt[3] + xt[1]);

    vt[0] += -0.5 * (yt[0] - yt[2]);
    vt[1] += -0.5 * (yt[1] - yt[3]);
    vt[2] += -0.5 * (yt[2] - yt[0]);
    vt[3] += -0.5 * (yt[3] - yt[1]);

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = a * vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] + vt[2]);
    yt[1] = 0.5 * (vt[1] + vt[3]);
    yt[2] = 0.5 * (vt[2] + vt[0]);
    yt[3] = 0.5 * (vt[3] + vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += - f[is] * yt[id];
      }

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
      }

    }

  } // idx loop

}

void mult_domainwall_5din_LUinv_dirac(
                        real_t* vp, real_t* wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);
  real_t* e_dev;
  real_t* f_dev;
  real_t* dpinv_dev;
  real_t* dm_dev;

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&e_dev , nsize));
  CHECK(cudaMalloc((void **)&f_dev , nsize));
  CHECK(cudaMalloc((void **)&dpinv_dev , nsize));
  CHECK(cudaMalloc((void **)&dm_dev , nsize));

  CHECK(cudaMemcpy(e_dev, e, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(f_dev, f, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dpinv_dev, dpinv, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dm_dev, dm, nsize, cudaMemcpyHostToDevice));

  CHECK(cudaDeviceSynchronize());
  int blockSize  = 2 * VECTOR_LENGTH;
  int gridSize   = (Nst * NVC + blockSize - 1) / blockSize;

  mult_domainwall_5din_LUinv_dirac_kernel<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, 
                 Ns, Nx, Ny, Nz, Nt,
                 e_dev,f_dev, dpinv_dev, dm_dev);

  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(e_dev));
  CHECK(cudaFree(f_dev));
  CHECK(cudaFree(dpinv_dev));
  CHECK(cudaFree(dm_dev));

}


__global__ void mult_domainwall_5din_LUinv_dirac_alt_kernel(
                        real_t* __restrict__ vp, real_t* __restrict__ wp,
                        int Ns, int Nx, int Ny, int Nz, int Nt,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  const int Nst     = Nx * Ny * Nz * Nt;
  const int Nin5    = NVCD * Ns;
  //const int Nst_pad = CEIL_NWP(Nst);
  const int idx     = blockIdx.x * blockDim.x + threadIdx.x;

  if ( idx < Nst * NVC * ND2 ){

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivcd1   = idx2_wp % (NVC * ND2);
    int idx_out = idx2_wp / (NVC * ND2);
    int site    = idx_in + NWP * idx_out;

    int ivcd2 = ivcd1 + NVC * 2;
    real_t vt1, vt2, yt1, yt2, xt1, xt2;
    int is = 0;

    vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

    yt1 = e[0] * vt1;
    yt2 = e[0] * vt2;

    for (int is = 1; is < Ns-1; ++is) {

      xt1 = vt1;
      xt2 = vt2;

      vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
      vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      vt1 += a * (xt1 + xt2);
      vt2 += a * (xt1 + xt1);

      vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
      vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

      yt1 += e[is] * vt1;
      yt2 += e[is] * vt2;

    }

    is = Ns-1;

    xt1 = vt1;
    xt2 = vt2;

    vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    vt1 += a * (xt1 + xt2);
    vt2 += a * (xt2 + xt1);

    vt1 += -0.5 * (yt1 - yt2);
    vt2 += -0.5 * (yt2 - yt1);

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];

    vt1 = a * vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = a * vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

    yt1 = 0.5 * (vt1 + vt2);
    yt2 = 0.5 * (vt2 + vt1);

    for (int is = Ns-2; is >= 0; --is) {

      xt1 = vt1;
      xt2 = vt2;

      vt1 = vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
      vt2 = vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

      real_t a = real_t(0.5) * dm[is];

      vt1 += a * (xt1 - xt2);
      vt2 += a * (xt2 - xt1);

      vt1 += - f[is] * yt1;
      vt2 += - f[is] * yt2;

      real_t aa = dpinv[is];

      vt1 *= aa;
      vt2 *= aa;

      vp[IDX2(Nin5, (ivcd1 + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ivcd2 + NVCD * is), site)] = vt2;

    }
  } // idx loop
}

void mult_domainwall_5din_LUinv_dirac_alt(
                        real_t* vp, real_t* wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  //int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  real_t* e_dev;
  real_t* f_dev;
  real_t* dpinv_dev;
  real_t* dm_dev;

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&e_dev , nsize));
  CHECK(cudaMalloc((void **)&f_dev , nsize));
  CHECK(cudaMalloc((void **)&dpinv_dev , nsize));
  CHECK(cudaMalloc((void **)&dm_dev , nsize));

  CHECK(cudaMemcpy(e_dev, e, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(f_dev, f, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dpinv_dev, dpinv, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dm_dev, dm, nsize, cudaMemcpyHostToDevice));

  int blockSize = 2* VECTOR_LENGTH;
  int gridSize   = (Nst * NC * Ns + blockSize - 1) / blockSize;

  mult_domainwall_5din_LUinv_dirac_alt_kernel<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, 
                 Ns, Nx, Ny, Nz, Nt,
                 e_dev,f_dev, dpinv_dev, dm_dev);

  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(e_dev));
  CHECK(cudaFree(f_dev));
  CHECK(cudaFree(dpinv_dev));
  CHECK(cudaFree(dm_dev));
}



__global__ void mult_domainwall_5din_LUdaginv_dirac_kernel(
                        real_t * __restrict__ vp, real_t *__restrict__ wp,
                        int Ns, int Nx, int Ny,int Nz, int Nt,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  const int Nst  = Nx * Ny * Nz * Nt;
  const int Nin5 = NVCD * Ns;
  const int idx  = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < Nst * NVC) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc     = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site    = idx_in + NWP * idx_out;

    // anchor
    __shared__ real_t vt[ND], yt[ND], xt[ND];
    int is = 0;

    real_t a = dpinv[0];

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = a * wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = f[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is - 1];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += f[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    a = real_t(0.5) * dm[is - 1];

    vt[0] += a * (xt[0] - xt[2]);
    vt[1] += a * (xt[1] - xt[3]);
    vt[2] += a * (xt[2] - xt[0]);
    vt[3] += a * (xt[3] - xt[1]);

    vt[0] += -0.5 * (yt[0] + yt[2]);
    vt[1] += -0.5 * (yt[1] + yt[3]);
    vt[2] += -0.5 * (yt[2] + yt[0]);
    vt[3] += -0.5 * (yt[3] + yt[1]);

    real_t aa = dpinv[is];

    for(int id = 0; id < ND; ++id){
      vt[id] *= aa;
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    // Udag_inv completed

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] - vt[2]);
    yt[1] = 0.5 * (vt[1] - vt[3]);
    yt[2] = 0.5 * (vt[2] - vt[0]);
    yt[3] = 0.5 * (vt[3] - vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += -e[is] * yt[id];
      }
      
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }
    }

  } // idx loop end

}


void mult_domainwall_5din_LUdaginv_dirac(
                        real_t* vp, real_t* wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  //int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  real_t* e_dev;
  real_t* f_dev;
  real_t* dpinv_dev;
  real_t* dm_dev;

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&e_dev , nsize));
  CHECK(cudaMalloc((void **)&f_dev , nsize));
  CHECK(cudaMalloc((void **)&dpinv_dev , nsize));
  CHECK(cudaMalloc((void **)&dm_dev , nsize));

  CHECK(cudaMemcpy(e_dev, e, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(f_dev, f, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dpinv_dev, dpinv, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(dm_dev, dm, nsize, cudaMemcpyHostToDevice));

  int blockSize = 2 * VECTOR_LENGTH;
  int gridSize   = (Nst * NVC + blockSize - 1) / blockSize;

  mult_domainwall_5din_LUdaginv_dirac_kernel<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, 
                 Ns, Nx, Ny, Nz, Nt,
                 e_dev,f_dev, dpinv_dev, dm_dev);

  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(e_dev));
  CHECK(cudaFree(f_dev));
  CHECK(cudaFree(dpinv_dev));
  CHECK(cudaFree(dm_dev));
}






#endif
//============================================================END=====

