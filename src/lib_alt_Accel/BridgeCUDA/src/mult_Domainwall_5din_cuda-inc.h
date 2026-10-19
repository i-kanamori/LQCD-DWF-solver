/*!
      @file    mult_Domainwall_5din_cuda-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-06-12 16:42:53 #$
      @version $LastChangedRevision: 2710 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC_INCLUDED

//====================================================================
__global__ 
void mult_domainwall_5din_5dir_dirac_kernel(
      real_t *vp, real_t *yp, real_t *wp,
      real_t mq, real_t M0, int Ns, real_t *b_dev, real_t *c_dev, int Nst) {

  int Nin5 = NVCD * Ns;
  extern __shared__ real_t shared_mem[];
  real_t* b = shared_mem;
  real_t* c = &shared_mem[Ns]; 

  int tid = threadIdx.x;
  if(tid < Ns) {
    b[tid] = b_dev[tid];
    c[tid] = c_dev[tid];
  }

  int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site < Nst) {
    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD], wt[NVCD];

      int is_up  = (is + 1) % Ns;
      real_t Fup = 0.5;
      if (is == Ns - 1) Fup = -0.5 * mq;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] = Fup * (wt[ID1 + ivc] - wt[ID3 + ivc]);
        vt[ID2 + ivc] = Fup * (wt[ID2 + ivc] - wt[ID4 + ivc]);
        vt[ID3 + ivc] = Fup * (wt[ID3 + ivc] - wt[ID1 + ivc]);
        vt[ID4 + ivc] = Fup * (wt[ID4 + ivc] - wt[ID2 + ivc]);
      }

      int is_dn = (is - 1 + Ns) % Ns;
      real_t Fdn = 0.5;
      if (is == 0) Fdn = -0.5 * mq;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] += Fdn * (wt[ID1 + ivc] + wt[ID3 + ivc]);
        vt[ID2 + ivc] += Fdn * (wt[ID2 + ivc] + wt[ID4 + ivc]);
        vt[ID3 + ivc] += Fdn * (wt[ID3 + ivc] + wt[ID1 + ivc]);
        vt[ID4 + ivc] += Fdn * (wt[ID4 + ivc] + wt[ID2 + ivc]);
      }

      real_t B_is = b[is] * (4.0 - M0) + 1.0;
      real_t C_is = c[is] * (4.0 - M0) - 1.0;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = B_is * wt[ivcd] + C_is * vt[ivcd];
        yp[IDX2(Nin5, (ivcd + NVCD * is), site)] = -0.5 * (b[is] * wt[ivcd] + c[is] * vt[ivcd]);
      }
    }
  }
}


void mult_domainwall_5din_5dir_dirac(
      real_t *vp, real_t *yp, real_t *wp,
      real_t mq, real_t M0, int Ns, real_t *b, real_t *c,  real_t alpha,
      int *Nsize) {

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  real_t* b_dev;
  real_t* c_dev; 

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&b_dev , nsize));
  CHECK(cudaMalloc((void **)&c_dev , nsize));

  CHECK(cudaMemcpy(b_dev, b, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(c_dev, c, nsize, cudaMemcpyHostToDevice));

  // Define block and grid sizes for CUDA kernel
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst + blockSize - 1)/ blockSize;
  int shared_mem_size = 2 * nsize;

  // Launch the CUDA kernel
  mult_domainwall_5din_5dir_dirac_kernel<<<gridSize, blockSize, shared_mem_size>>>(
      vp_dev, yp_dev, wp_dev, mq, M0, Ns, b_dev, c_dev, Nst);

  // Synchronize to ensure completion of kernel execution
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(b_dev));
  CHECK(cudaFree(c_dev));
}

//====================================================================
__global__ 
void mult_domainwall_5din_5dirdag_dirac_kernel(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns, real_t *b_dev, real_t *c_dev, int Nst) 
{

    int Nin5 = NVCD * Ns;

    extern __shared__ real_t shared_mem[];
    real_t* b = shared_mem;
    real_t* c = &shared_mem[Ns]; 

    int tid = threadIdx.x;
    if(tid < Ns) {
        b[tid] = b_dev[tid];
        c[tid] = c_dev[tid];
    }

    int site = blockIdx.x * blockDim.x + threadIdx.x;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {
            real_t vt[NVCD], xt[NVCD], wt[NVCD], yt[NVCD];

            real_t B1 = b[is] * (4.0 - M0) + 1.0;
            real_t a1 = -0.5 * b[is];

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
                yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] = B1 * wt[ID1 + ivc] + a1 * yt[ID3 + ivc];
                vt[ID2 + ivc] = B1 * wt[ID2 + ivc] + a1 * yt[ID4 + ivc];
                vt[ID3 + ivc] = B1 * wt[ID3 + ivc] + a1 * yt[ID1 + ivc];
                vt[ID4 + ivc] = B1 * wt[ID4 + ivc] + a1 * yt[ID2 + ivc];
            }

            int is_up = (is + 1) % Ns;
            real_t C1 = c[is_up] * (4.0 - M0) - 1.0;
            real_t aup = -0.5 * c[is_up];

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
                yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                xt[ID1 + ivc] = C1 * wt[ID1 + ivc] + aup * yt[ID3 + ivc];
                xt[ID2 + ivc] = C1 * wt[ID2 + ivc] + aup * yt[ID4 + ivc];
                xt[ID3 + ivc] = C1 * wt[ID3 + ivc] + aup * yt[ID1 + ivc];
                xt[ID4 + ivc] = C1 * wt[ID4 + ivc] + aup * yt[ID2 + ivc];
            }

            real_t Fup = 0.5;
            if (is == Ns-1) Fup = -0.5 * mq;

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] += Fup * (xt[ID1 + ivc] + xt[ID3 + ivc]);
                vt[ID2 + ivc] += Fup * (xt[ID2 + ivc] + xt[ID4 + ivc]);
                vt[ID3 + ivc] += Fup * (xt[ID3 + ivc] + xt[ID1 + ivc]);
                vt[ID4 + ivc] += Fup * (xt[ID4 + ivc] + xt[ID2 + ivc]);
            }

            int is_dn = (is - 1 + Ns) % Ns;

            real_t C2  = c[is_dn] * (4.0 - M0) - 1.0;
            real_t adn = -0.5 * c[is_dn];

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
                yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                xt[ID1 + ivc] = C2 * wt[ID1 + ivc] + adn * yt[ID3 + ivc];
                xt[ID2 + ivc] = C2 * wt[ID2 + ivc] + adn * yt[ID4 + ivc];
                xt[ID3 + ivc] = C2 * wt[ID3 + ivc] + adn * yt[ID1 + ivc];
                xt[ID4 + ivc] = C2 * wt[ID4 + ivc] + adn * yt[ID2 + ivc];
            }

            real_t Fdn = 0.5;
            if (is == 0) Fdn = -0.5 * mq;

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] += Fdn * (xt[ID1 + ivc] - xt[ID3 + ivc]);
                vt[ID2 + ivc] += Fdn * (xt[ID2 + ivc] - xt[ID4 + ivc]);
                vt[ID3 + ivc] += Fdn * (xt[ID3 + ivc] - xt[ID1 + ivc]);
                vt[ID4 + ivc] += Fdn * (xt[ID4 + ivc] - xt[ID2 + ivc]);
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_5dirdag_dirac(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, real_t alpha,
       int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst  = Nx * Ny * Nz * Nt;

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* yp_dev = (real_t*)dev_ptr(yp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    real_t* b_dev;
    real_t* c_dev; 

    size_t nsize = Ns * sizeof(real_t);

    CHECK(cudaMalloc((void **)&b_dev , nsize));
    CHECK(cudaMalloc((void **)&c_dev , nsize));

    CHECK(cudaMemcpy(b_dev, b, nsize, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(c_dev, c, nsize, cudaMemcpyHostToDevice));

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_5dirdag_dirac_kernel<<<gridSize, blockSize, 2 * nsize>>>(
        vp_dev, yp_dev, wp_dev, mq, M0, Ns, b_dev, c_dev, Nst);

    CHECK(cudaDeviceSynchronize());
    CHECK(cudaFree(b_dev));
    CHECK(cudaFree(c_dev));
}

//====================================================================
__global__ void mult_domainwall_5din_mult_gm5_dirac_kernel(
         real_t *vp, real_t *wp, int Ns, int Nst) {

    int Nin5 = NVCD * Ns;        
    int site = blockIdx.x * blockDim.x + threadIdx.x; 
    
    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            real_t vt[NVCD], wt[NVCD]; 

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] = wt[ID3 + ivc];
                vt[ID2 + ivc] = wt[ID4 + ivc];
                vt[ID3 + ivc] = wt[ID1 + ivc];
                vt[ID4 + ivc] = wt[ID2 + ivc];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_gm5_dirac(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_gm5_dirac_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst);

    CHECK(cudaDeviceSynchronize());
}


//====================================================================
__global__ void mult_domainwall_5din_mult_R_kernel(
         real_t *vp, real_t *wp, int Ns, int Nst) {

    int Nin5 = NVCD * Ns; 
    int site = blockIdx.x * blockDim.x + threadIdx.x; 

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            real_t vt[NVCD]; 
            int isR = Ns - 1 - is; 

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * isR), site)];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_R(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_R_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst);

    CHECK(cudaDeviceSynchronize());
}
//----------------------------------------------------------------
__global__ void mult_domainwall_5din_mult_gm5R_dirac_kernel(
         real_t *vp, real_t *wp, int Ns, int Nx, int Ny, int Nz, int Nt) {
    
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {

            int isR = Ns - 1 - is;
            real_t vt[NVCD], wt[NVCD];

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * isR), site)];
            }

            for (int ivc = 0; ivc < NVC; ++ivc) {
                vt[ID1 + ivc] = wt[ID3 + ivc];
                vt[ID2 + ivc] = wt[ID4 + ivc];
                vt[ID3 + ivc] = wt[ID1 + ivc];
                vt[ID4 + ivc] = wt[ID2 + ivc];
            }

            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_mult_gm5R_dirac(
         real_t *vp, real_t *wp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt;

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_mult_gm5R_dirac_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nx, Ny, Nz, Nt);

    cudaDeviceSynchronize();
}

//====================================================================
__global__ void mult_domainwall_5din_clear_kernel(real_t *vp, int Ns, int Nst) {

    // Set each element of vp to 0.0
    
    int site = blockIdx.x * blockDim.x + threadIdx.x; 
    int Nin5 = NVCD * Ns;

    if (site < Nst) {
        for (int is = 0; is < Ns; ++is) {
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = 0.0;
            }
        }
    }
}

void mult_domainwall_5din_clear(real_t *vp, int Ns, int *Nsize) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);

    int blockSize = VECTOR_LENGTH;
    int gridSize = (Nst + blockSize - 1)/ blockSize;

    mult_domainwall_5din_clear_kernel<<<gridSize, blockSize>>>(vp_dev, Ns, Nst);

    CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ 
void mult_domainwall_5din_hopb_5D_dirac_kernel(
     real_t *vp, real_t *up, real_t *wp, int Ns, 
     int bc_x, int bc_y, int bc_z, int bc_t,
     int Nx, int Ny, int Nz, int Nt,
     int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
     int flag, int Nst_pad) 
{
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx   * Ny;
    int Nxyz = Nxy  * Nz;
    int Nst  = Nx   * Ny * Nz * Nt;
    int idx  = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t * u_up = up;
    real_t * u_dn = up;

    if (idx < Nst_pad * Ns) {

        int idx2_wp = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx2_wp % Ns;
        int idx_out = idx2_wp / Ns;
        int site    = idx_in + NWP * idx_out;

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt  % Nz;
        int it   = izt  / Nz;
        int ixyz = site % Nxyz;
        int idir;

        __shared__ real_t u_0, u_1, u_2, u_3, u_4, u_5;
        __shared__ real_t u_6, u_7, u_8, u_9, u10, u11;
        __shared__ real_t u12, u13, u14, u15, u16, u17;
        real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
        real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
        real_t wt1r, wt1i, wt2r, wt2i;

        real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
        real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
        real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
        real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

        #include "inc/mult_Domainwall_cuda_flag-inc.h"

        real_t * v1 = wp;
        real_t * v2 = vp;

        #include "inc/mult_Domainwall_cuda_hopb_dirac-inc.h"

        v2[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
        v2[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
        v2[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
        v2[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
        v2[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
        v2[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

        v2[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
        v2[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
        v2[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
        v2[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
        v2[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
        v2[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

        v2[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
        v2[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
        v2[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
        v2[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
        v2[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
        v2[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

        v2[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
        v2[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
        v2[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
        v2[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
        v2[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
        v2[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;
    }
}

void mult_domainwall_5din_hopb_dirac_5d(
     real_t *vp, real_t *up, real_t *wp, int Ns, int *bc, int *Nsize, int *do_comm, int flag) {

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 
    int Nst_pad = ceil_nwp(Nst);  

    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];
    int do_comm_x = do_comm[0];
    int do_comm_y = do_comm[1];
    int do_comm_z = do_comm[2];
    int do_comm_t = do_comm[3];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst_pad * Ns + blockSize - 1) / blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t); 

    mult_domainwall_5din_hopb_5D_dirac_kernel<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns,
    bc_x, bc_y, bc_z, bc_t,
    Nx, Ny, Nz, Nt,
    do_comm_x, do_comm_y, do_comm_z, do_comm_t,
    flag, Nst_pad); 
    
    CHECK(cudaDeviceSynchronize());
}


//====================================================================
__global__ void mult_domainwall_5din_hop1x_5D_dirac_dev(
    real_t *buf_xp, real_t *buf_xm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstx_pad) 
{
    int Nyzt   = Ny * Nz * Nt;
    int Nin5   = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst    = Nx * Ny * Nz * Nt;
    int idx    = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstx_pad * Ns) {

        int idx_wp  = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx_wp % Ns;
        int idx_out = idx_wp / Ns;
        int iyzt    = idx_in + NWP * idx_out;

        {
            int ix = 0;
            int site = ix + Nx * iyzt;
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_xp1(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = bc2 * vt[ivcd];
            }
        }
        {
            int ix = Nx - 1;
            int site = ix + Nx * iyzt;
            load_u(ut, up, site + Nst * 0);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_xm1(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)] = vt[ivcd];
            }
        }
    }
}

__global__ void mult_domainwall_5din_hop1y_5D_dirac_dev(
    real_t *buf_yp, real_t *buf_ym,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nxzt_pad) 
{
    int Nxzt = Nx * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];

    real_t* ut = &sharedMemory[0];
    real_t* wt = &sharedMemory[NVCD];

    if (idx < Nxzt_pad * Ns) {
 
        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixzt     = idx_in + NWP * idx_out;

        {
            int iy   = 0;
            int ix   = ixzt % Nx;
            int izt  = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);

            real_t bc2 = bc;
            real_t vt[NVC * ND2];

            #pragma unroll
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2_5D(NVCD, ivcd, is, Ns,site)];
            }
            mult_wilson_yp1(vt, wt);

            #pragma unroll
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_yp[IDX2_5D(NVC * ND2, ivcd, is, Ns, ixzt)]= bc2 * vt[ivcd];
            }
        }

        {
            int iy = Ny - 1;
            int ix = ixzt % Nx;
            int izt = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);
            load_u(ut, up, site + Nst * 1);

            real_t vt[NVC * ND2];
            #pragma unroll
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2_5D(NVCD, ivcd, is, Ns,site)];
            }
            mult_wilson_ym1(vt, ut, wt);
            
            #pragma unroll
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_ym[IDX2_5D(NVC * ND2, ivcd, is, Ns, ixzt)] = vt[ivcd];
            }
        }
    }
}

__global__ void mult_domainwall_5din_hop1z_5D_dirac_dev(
    real_t *buf_zp, real_t *buf_zm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstz_pad) 
{
    int Nxy = Nx * Ny;
    int Nxyt = Nx * Ny * Nt;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstz_pad * Ns) {

        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixyt     = idx_in + NWP * idx_out;

        {
            int iz = 0;
            int ixy = ixyt % Nxy;
            int it = ixyt / Nxy;
            int site = ixy + Nxy * (iz + Nz * it);
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_zp1(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = bc2 * vt[ivcd];
            }
        }

        {
            int iz = Nz - 1;
            int ixy = ixyt % Nxy;
            int it = ixyt / Nxy;
            int site = ixy + Nxy * (iz + Nz * it);
            load_u(ut, up, site + Nst * 2);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_zm1(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)] = vt[ivcd];
            }
        }
    }
}
__global__ void mult_domainwall_5din_hop1t_5D_dirac_dev(
    real_t *buf_tp, real_t *buf_tm,
    real_t *up, real_t *wp,
    int Ns, int bc, int Nx, int Ny, int Nz, int Nt, int Nstt_pad) 
{
    int Nxyz = Nx * Ny * Nz;
    int Nin5 = NVCD * Ns;
    int Nin5bd = NVC * ND2 * Ns;
    int Nst = Nx * Ny * Nz * Nt;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nstt_pad * Ns) {

        int idx_wp   = idx / NWP;
        int idx_in   = idx % NWP;
        int is       = idx_wp % Ns;
        int idx_out  = idx_wp / Ns;
        int ixyz     = idx_in + NWP * idx_out;

        {
            int it = 0;
            int site = ixyz + Nxyz * it;
            real_t bc2 = bc;

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_tp1_dirac(vt, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = bc2 * vt[ivcd];
            }
        }

        {
            int it = Nt - 1;
            int site = ixyz + Nxyz * it;
            load_u(ut, up, site + Nst * 3);

            real_t wt[NVCD], vt[NVC * ND2];
            for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
                wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
            mult_wilson_tm1_dirac(vt, ut, wt);
            for (int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)] = vt[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_hop1_dirac(
    real_t * buf1_xp, real_t * buf1_xm,
    real_t * buf1_yp, real_t * buf1_ym,
    real_t * buf1_zp, real_t * buf1_zm,
    real_t * buf1_tp, real_t * buf1_tm,
    real_t * up, real_t * wp,
    int Ns, int *bc, int *Nsize, int *do_comm){

    int Nx     = Nsize[0];
    int Ny     = Nsize[1];
    int Nz     = Nsize[2];
    int Nt     = Nsize[3];
    int Nst    = Nx * Ny * Nz * Nt;
    int Nin5bd = NVC * Ns * ND2;

    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];

    int do_comm_x = do_comm[0];
    int do_comm_y = do_comm[1];
    int do_comm_z = do_comm[2];
    int do_comm_t = do_comm[3];

    size_t size_bx = Nin5bd * Ny * Nz * Nt;
    size_t size_by = Nin5bd * Nx * Nz * Nt;
    size_t size_bz = Nin5bd * Ny * Nx * Nt;
    size_t size_bt = Nin5bd * Ny * Nz * Nx;

    real_t * up_dev = (real_t*)dev_ptr(up);
    real_t * wp_dev = (real_t*)dev_ptr(wp);

    //int Nst_pad = ceil_nwp(Nst);  

    if( do_comm_x > 0 ){

        int Nstx =  Ny * Nz * Nt;
        int Nstx_pad = ceil_nwp(Nstx);  

        real_t * buf1_xp_dev = (real_t*)dev_ptr(buf1_xp);
        real_t * buf1_xm_dev = (real_t*)dev_ptr(buf1_xm);

        int blockSize   =  VECTOR_LENGTH;
        int gridSize5D  =  (Nstx * Ns + blockSize - 1)/ blockSize;
        int sharedsize  = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_hop1x_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_xp_dev, buf1_xm_dev,
        up_dev, wp_dev,
        Ns, bc_x, Nx, Ny, Nz, Nt, Nstx_pad);

/*
        int gridSize4D  =  (Nstx + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1x_dirac_dev<<<gridSize4D, blockSize, sharedsize>>>(
        buf1_xp_dev, buf1_xm_dev,
        up_dev, wp_dev,
        Ns, bc_x, Nx, Ny, Nz, Nt);
*/
    }

    if( do_comm_y > 0){

        int Nsty =  Nx * Nz * Nt;
        int Nsty_pad = ceil_nwp(Nsty);  

        real_t * buf1_yp_dev = (real_t*)dev_ptr(buf1_yp);
        real_t * buf1_ym_dev = (real_t*)dev_ptr(buf1_ym);

        int blockSize = VECTOR_LENGTH;
        int gridSize5D  =  (Nsty * Ns + blockSize - 1)/ blockSize;
        int sharedsize = (NDF + NVCD) * blockSize * sizeof(real_t);

/*     
        mult_domainwall_5din_hop1y_5D_try_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_yp_dev, buf1_ym_dev,
        up_dev, wp_dev,
        Ns, bc_y, Nx, Ny, Nz, Nt, Nsty_pad);
*/

        mult_domainwall_5din_hop1y_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_yp_dev, buf1_ym_dev,
        up_dev, wp_dev,
        Ns, bc_y, Nx, Ny, Nz, Nt, Nsty_pad);



    }

    if( do_comm_z > 0 ){

        int Nstz =  Nx * Ny * Nt;
        int Nstz_pad = ceil_nwp(Nstz);  

        real_t * buf1_zp_dev = (real_t*)dev_ptr(buf1_zp);
        real_t * buf1_zm_dev = (real_t*)dev_ptr(buf1_zm);

        int blockSize = VECTOR_LENGTH;
        int gridSize5D  =  (Nstz_pad * Ns + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_hop1z_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_zp_dev, buf1_zm_dev,
        up_dev, wp_dev,
        Ns, bc_z, Nx, Ny, Nz, Nt, Nstz_pad);

/*
        int gridSize4D  = (Nstz + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1z_dirac_dev<<<gridSize4D, blockSize, sharedsize>>>(
        buf1_zp_dev, buf1_zm_dev,
        up_dev, wp_dev,
        Ns, bc_z, Nx, Ny, Nz, Nt);
*/
    }

    if( do_comm_t > 0 ){

        int Nstt =  Nx * Ny * Nz;
        int Nstt_pad = ceil_nwp(Nstt);  

        real_t * buf1_tp_dev = (real_t*)dev_ptr(buf1_tp);
        real_t * buf1_tm_dev = (real_t*)dev_ptr(buf1_tm);

        int blockSize = VECTOR_LENGTH;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        int gridSize5D  =  (Nstt_pad * Ns + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1t_5D_dirac_dev<<<gridSize5D, blockSize, sharedsize>>>(
        buf1_tp_dev, buf1_tm_dev,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt, Nstt_pad);

/*
        int gridSize4D  = (Nstt + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop1t_dirac_dev<<<gridSize4D, blockSize, sharedsize>>>(
        buf1_tp_dev, buf1_tm_dev,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt);
*/
    }

    CHECK(cudaDeviceSynchronize());
#ifdef USE_MPI_GPUAWARE
    // do nothing
#else
    cudaStream_t* cuda_stream = nullptr;
    int nstream = 8;
    CREATE_CUDA_STREAMS(&cuda_stream, nstream);

    if (cuda_stream == nullptr) {
        fprintf(stderr, "Failed to create CUDA streams in hop1.\n");
    }

    if ( do_comm_x > 0) {
        update_host_asy(buf1_xp, 0 , size_bx, cuda_stream[0]);
        update_host_asy(buf1_xm, 0 , size_bx, cuda_stream[1]);
    } 
        
    if ( do_comm_y > 0) {
        update_host_asy(buf1_yp, 0 , size_by, cuda_stream[2]);
        update_host_asy(buf1_ym, 0 , size_by, cuda_stream[3]);
    } 

    if ( do_comm_z > 0) {
        update_host_asy(buf1_zp, 0 , size_bz, cuda_stream[4]);
        update_host_asy(buf1_zm, 0 , size_bz, cuda_stream[5]);
    } 

    if ( do_comm_t > 0) {
        update_host_asy(buf1_tp, 0 , size_bt, cuda_stream[6]);
        update_host_asy(buf1_tm, 0 , size_bt, cuda_stream[7]);
    } 

    DESTROY_CUDA_STREAMS(&cuda_stream, nstream);
#endif

}

//====================================================================
__global__ void mult_domainwall_5din_hop2_4D_dirac_dev(
    real_t *vp, real_t *up, real_t *wp,
    real_t *buf_xp, real_t *buf_xm,
    real_t *buf_yp, real_t *buf_ym,
    real_t *buf_zp, real_t *buf_zm,
    real_t *buf_tp, real_t *buf_tm,
    int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t, 
    int Nx, int Ny, int Nz, int Nt) {

    int site         = blockIdx.x * blockDim.x + threadIdx.x;
    int Nxy          = Nx   * Ny           ;
    int Nxyz         = Nx   * Ny  * Nz     ;
    int Nin5bd       = (NVCD / 2) * Ns     ;
    int Nst          = Nx   * Ny  * Nz * Nt;
    int Nin5         = NVCD * Ns           ;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt  % Nz;
        int it   = izt  / Nz;
        int ixyz = site % Nxyz;
        int idir;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vL[ivcd] = 0.0;
            }
            
            int opr_any = 0;

            idir = 0;
            if (do_comm_x > 0) {

                if (ix == Nx-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xp2(vL, ut, wt);
                    ++opr_any;
                }

                if (ix == 0) {
                    real_t bc2 = bc_x;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xm2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 1;
            if (do_comm_y > 0) {
                int ixzt = ix + Nx * izt;
                if (iy == Ny-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                    }
                    mult_wilson_yp2(vL, ut, wt);
                    ++opr_any;
                }

                if (iy == 0) {
                    real_t bc2 = bc_y;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
                    }
                    mult_wilson_ym2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 2;
            if (do_comm_z > 0) {
                int ixyt = ixy + Nxy * it;
                if (iz == Nz-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zp2(vL, ut, wt);
                    ++opr_any;
                }

                if (iz == 0) {
                    real_t bc2 = bc_z;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zm2(vL, wt);
                    ++opr_any;
                }
            }

            idir = 3;
            if (do_comm_t > 0) {
                if (it == Nt-1) {
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tp2_dirac(vL, ut, wt);
                    ++opr_any;
                }

                if (it == 0) {
                    real_t bc2 = bc_t;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tm2_dirac(vL, wt);
                    ++opr_any;
                }
            }

            if (opr_any > 0) {
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	            }
            }
        }   // is loop
    }
}

__global__ void mult_domainwall_5din_hop2_5D_dirac_dev(
    real_t *vp, real_t *up, real_t *wp,
    real_t *buf_xp, real_t *buf_xm,
    real_t *buf_yp, real_t *buf_ym,
    real_t *buf_zp, real_t *buf_zm,
    real_t *buf_tp, real_t *buf_tm,
    int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t, 
    int Nx, int Ny, int Nz, int Nt, int Nst_pad) 
{
    int Nxy    = Nx * Ny;
    int Nxyz   = Nx * Ny * Nz;
    int Nin5bd = (NVCD / 2) * Ns;
    int Nst    = Nx * Ny * Nz * Nt;
    int Nin5   = NVCD * Ns;
    int idx    = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (idx < Nst * Ns) {

        int idx2_wp = idx / NWP;
        int idx_in  = idx % NWP;
        int is      = idx2_wp % Ns;
        int idx_out = idx2_wp / Ns;
        int site    = idx_in + NWP * idx_out;

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int ixy  = site % Nxy;
        int iy   = iyzt % Ny;
        int izt  = site / Nxy;
        int iz   = izt % Nz;
        int it   = izt / Nz;
        int ixyz = site % Nxyz;

        real_t vL[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
            vL[ivcd] = 0.0;
        }

        int opr_any = 0;
        int idir;

        idir = 0;
        if (do_comm_x > 0) {

            if (ix == Nx-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)];
                }
                mult_wilson_xp2(vL, ut, wt);
                ++opr_any;
            }

            if (ix == 0) {
                real_t bc2 = bc_x;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), iyzt)];
                }
                mult_wilson_xm2(vL, wt);
                ++opr_any;
            }
        }

        idir = 1;
        if (do_comm_y > 0) {
            int ixzt = ix + Nx * izt;
            if (iy == Ny-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                }
                mult_wilson_yp2(vL, ut, wt);
                ++opr_any;
            }

            if (iy == 0) {
                real_t bc2 = bc_y;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                }
                mult_wilson_ym2(vL, wt);
                ++opr_any;
            }
        }

        idir = 2;
        if (do_comm_z > 0) {
            int ixyt = ixy + Nxy * it;
            if (iz == Nz-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
                }
                mult_wilson_zp2(vL, ut, wt);
                ++opr_any;
            }

            if (iz == 0) {
                real_t bc2 = bc_z;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyt)];
                }
                mult_wilson_zm2(vL, wt);
                ++opr_any;
            }
        }

        idir = 3;
        if (do_comm_t > 0) {
            if (it == Nt-1) {
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
                }
                mult_wilson_tp2_dirac(vL, ut, wt);
                ++opr_any;
            }

            if (it == 0) {
                real_t bc2 = bc_t;
                real_t wt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd) {
                    wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixyz)];
                }
                mult_wilson_tm2_dirac(vL, wt);
                ++opr_any;
            }
        }

        if (opr_any > 0) {
            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vp[IDX2(Nin5, (ivcd + NVCD * is), site)] += vL[ivcd];
            }
        }
    }
}

void mult_domainwall_5din_hop2_dirac(
        real_t * vp     , real_t * up, real_t * wp,
        real_t * buf2_xp, real_t * buf2_xm,
        real_t * buf2_yp, real_t * buf2_ym,
        real_t * buf2_zp, real_t * buf2_zm,
        real_t * buf2_tp, real_t * buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int Nin5bd = (NVCD / 2) * Ns;
        int Nst_pad = ceil_nwp(Nst);  

        int bc_x = bc[0];
        int bc_y = bc[1];
        int bc_z = bc[2];
        int bc_t = bc[3];
        int do_comm_x = do_comm[0];
        int do_comm_y = do_comm[1];
        int do_comm_z = do_comm[2];
        int do_comm_t = do_comm[3];
#ifdef USE_MPI_GPUAWARE
        // do nothing
#else
        cudaStream_t* cuda_stream = nullptr;
        int nstream = 4;
        CREATE_CUDA_STREAMS(&cuda_stream, nstream);

        if (cuda_stream == nullptr) {
            fprintf(stderr, "Failed to create CUDA streams.\n");
        }

        if( do_comm[0] > 0){
            size_t size_bx = Nin5bd * Ny * Nz * Nt;
            update_device_asy(buf2_xp, 0 ,size_bx, cuda_stream[0]);
            update_device_asy(buf2_xm, 0 ,size_bx, cuda_stream[0]);
        }

        if( do_comm[1] > 0){
            size_t size_by = Nin5bd * Nx * Nz * Nt;
            update_device_asy(buf2_yp, 0 ,size_by, cuda_stream[1]);
            update_device_asy(buf2_ym, 0 ,size_by, cuda_stream[1]);
        }

        if( do_comm[2] > 0){
            size_t size_bz = Nin5bd * Ny * Nx * Nt;
            update_device_asy(buf2_zp, 0 ,size_bz, cuda_stream[2]);
            update_device_asy(buf2_zm, 0 ,size_bz, cuda_stream[2]);
        }

        if( do_comm[3] > 0){
            size_t size_bt = Nin5bd * Ny * Nx * Nz;
            update_device_asy(buf2_tp, 0 ,size_bt, cuda_stream[3]);
            update_device_asy(buf2_tm, 0 ,size_bt, cuda_stream[3]);
        }

        DESTROY_CUDA_STREAMS(&cuda_stream, nstream);
#endif
        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);

        real_t * buf2_xp_dev = (real_t*)dev_ptr(buf2_xp);
        real_t * buf2_xm_dev = (real_t*)dev_ptr(buf2_xm);
        real_t * buf2_yp_dev = (real_t*)dev_ptr(buf2_yp);
        real_t * buf2_ym_dev = (real_t*)dev_ptr(buf2_ym);
        real_t * buf2_zp_dev = (real_t*)dev_ptr(buf2_zp);
        real_t * buf2_zm_dev = (real_t*)dev_ptr(buf2_zm);
        real_t * buf2_tp_dev = (real_t*)dev_ptr(buf2_tp);
        real_t * buf2_tm_dev = (real_t*)dev_ptr(buf2_tm);

        int blockSize  = VECTOR_LENGTH;

        int sharedsize = NDF * blockSize * sizeof(real_t);
/*
        int gridSize4D   = (Nst + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop2_4D_dirac_dev<<<gridSize4D, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,  wp_dev,
            buf2_xp_dev, buf2_xm_dev,
            buf2_yp_dev, buf2_ym_dev,
            buf2_zp_dev, buf2_zm_dev,
            buf2_tp_dev, buf2_tm_dev,
            Ns,  
            bc_x, bc_y, bc_z, bc_t,
            do_comm_x, do_comm_y, do_comm_z, do_comm_t,
            Nx, Ny, Nz, Nt);
*/
        int gridSize5D   = (Nst_pad * Ns + blockSize - 1)/ blockSize;
        mult_domainwall_5din_hop2_5D_dirac_dev<<<gridSize5D, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,  wp_dev,
            buf2_xp_dev, buf2_xm_dev,
            buf2_yp_dev, buf2_ym_dev,
            buf2_zp_dev, buf2_zm_dev,
            buf2_tp_dev, buf2_tm_dev,
            Ns,  
            bc_x, bc_y, bc_z, bc_t,
            do_comm_x, do_comm_y, do_comm_z, do_comm_t,
            Nx, Ny, Nz, Nt, Nst_pad);


        CHECK(cudaDeviceSynchronize());
}

#endif
//============================================================END=====
