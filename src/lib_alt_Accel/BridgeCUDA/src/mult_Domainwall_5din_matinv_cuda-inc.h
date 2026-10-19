/*!
      @file    mult_Doainwall_5din_matinv_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2710 $
*/

#ifndef MULT_DOMAINWALL_5DIN_MATINV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_MATINV_ACC_INCLUDED

//#define IDX_MAT(idx1, idx2, Ns)  ((idx1) + (ND2 * Ns) * (idx2)))

//====================================================================
__global__
void mult_domainwall_5din_ee_inv_dirac_4d_kernel(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;

  //int site = blockIdx.x * blockDim.x + threadIdx.x;
  //if (site < Nst) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < Nst * NC) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int idx_out = idx2_wp / NC;
    int site    = idx_in + NWP * idx_out;

    for (int is1 = 0; is1 < Ns; ++is1) {
      // for (int ic = 0; ic < NC; ++ic) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2;
            real_t mat1, mat2;

	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
            }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
            }	    

            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
    }
  }

}

void mult_domainwall_5din_ee_inv_dirac_4d(
                        real_t * vp, real_t * wp, int jd,
                        int Ns, real_t * mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);
  //real_t* mat_inv_dev = (real_t*)dev_ptr(mat_inv);

  int mat_size = ND2 * Ns;
  size_t mat_size2 = mat_size * mat_size * sizeof(real_t);
  real_t* mat_inv_dev;
  CHECK(cudaMalloc((void **)&mat_inv_dev, mat_size2));
  CHECK(cudaMemcpy(mat_inv_dev, mat_inv, mat_size2, cudaMemcpyHostToDevice));

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst * NC + blockSize - 1) / blockSize;

  mult_domainwall_5din_ee_inv_dirac_4d_kernel<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, jd, Ns, mat_inv_dev, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

  CHECK(cudaFree(mat_inv_dev));

}

//====================================================================
__global__
void mult_domainwall_5din_ee_inv_dirac_5d_kernel(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
// still not 5d implementation
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;
  int Nst_pad = CEIL_NWP(Nst);
  
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < Nst_pad * NC * Ns) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site = idx_in + NWP * idx_out;

    /*
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is1     = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;
    int site    = idx_in + NWP * idx_out;
    */
    //for (int is1 = 0; is1 < Ns; ++is1) {
      //for (int ic = 0; ic < NC; ++ic) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2;
            real_t mat1, mat2;

	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
            }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
            }	    

            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
      // }
  }

}

__global__
void mult_domainwall_5din_ee_inv_dirac_5d_kernel_opt(
                                   real_t * vp, real_t * wp,
                                   int jd, int Ns, real_t * __restrict__ mat_inv,
                                   int Nx, int Ny, int Nz, int Nt)
{
  const int Nst     = Nx * Ny * Nz * Nt;
  const int Nin5    = NVCD * Ns;
  const int Nst_pad = CEIL_NWP(Nst);
  const int idx     = blockIdx.x * blockDim.x + threadIdx.x;
  //__shared__ real_t shared_wp[8];
  __shared__ real_t mat1, mat2;
  __shared__ real_t wp1r, wp1i, wp2r, wp2i, wp3r, wp3i, wp4r, wp4i;

  if (idx < Nst_pad * NC * Ns) {

    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site    = idx_in + NWP * idx_out;

    for (int id = 0; id < ND2; ++id) {

      real_t vt1r = 0.0;
      real_t vt1i = 0.0;
      real_t vt2r = 0.0;
      real_t vt2i = 0.0;
      int idx1 = id + ND2 * is1;
      for (int is2 = 0; is2 < Ns; ++is2) {
        int idx2, ivcd1;

	      if(jd == 1){
          idx2 = 0 + ND2 * is2;
          mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
          idx2 = 1 + ND2 * is2;
          mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
        }else{
          idx2 = 0 + ND2 * is2;
          mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
          idx2 = 1 + ND2 * is2;
          mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
        }	    

        int ivcd1r = 0 + 2 * (ic + NC * 0);
        int ivcd1i = 0 + 2 * (ic + NC * 2);
        int ivcd2r = 1 + 2 * (ic + NC * 0);
        int ivcd2i = 1 + 2 * (ic + NC * 2);
        int ivcd3r = 0 + 2 * (ic + NC * 1);
        int ivcd3i = 0 + 2 * (ic + NC * 3);
        int ivcd4r = 1 + 2 * (ic + NC * 1);
        int ivcd4i = 1 + 2 * (ic + NC * 3);
        wp1r = wp[IDX2(Nin5, (ivcd1r + NVCD * is2), site)];
        wp1i = wp[IDX2(Nin5, (ivcd1i + NVCD * is2), site)];
        wp2r = wp[IDX2(Nin5, (ivcd2r + NVCD * is2), site)];
        wp2i = wp[IDX2(Nin5, (ivcd2i + NVCD * is2), site)];
        wp3r = wp[IDX2(Nin5, (ivcd3r + NVCD * is2), site)];
        wp3i = wp[IDX2(Nin5, (ivcd3i + NVCD * is2), site)];
        wp4r = wp[IDX2(Nin5, (ivcd4r + NVCD * is2), site)];
        wp4i = wp[IDX2(Nin5, (ivcd4i + NVCD * is2), site)];
        /*
        auto pipe = cuda::make_pipeline();
        pipe.producer_acquire();
        cuda::memcpy_async(&shared_wp[0], &wp[IDX2(Nin5, (ivcd1r + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[1], &wp[IDX2(Nin5, (ivcd1i + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[2], &wp[IDX2(Nin5, (ivcd2r + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[3], &wp[IDX2(Nin5, (ivcd2i + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[4], &wp[IDX2(Nin5, (ivcd3r + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[5], &wp[IDX2(Nin5, (ivcd3i + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[6], &wp[IDX2(Nin5, (ivcd4r + NVCD * is2), site)], sizeof(real_t), pipe);
        cuda::memcpy_async(&shared_wp[7], &wp[IDX2(Nin5, (ivcd4i + NVCD * is2), site)], sizeof(real_t), pipe);
        pipe.producer_commit();
        pipe.consumer_wait();
        __syncthreads();
        real_t wp1r = shared_wp[0];
        real_t wp1i = shared_wp[1];
        real_t wp2r = shared_wp[2];
        real_t wp2i = shared_wp[3];
        real_t wp3r = shared_wp[4];
        real_t wp3i = shared_wp[5];
        real_t wp4r = shared_wp[6];
        real_t wp4i = shared_wp[7];
        */
        vt1r +=  mat1 * wp1r + mat2 * wp1i;
        vt1i +=  mat1 * wp2r + mat2 * wp2i;
        vt2r +=  mat1 * wp3r + mat2 * wp3i;
        vt2i +=  mat1 * wp4r + mat2 * wp4i;
	    }
      int id1 = 2 * id;
      int id2 = 2 * id + 1;
      vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
      vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
      vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
      vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	  }
  }
}

void mult_domainwall_5din_ee_inv_dirac_5d(
                        real_t * vp, real_t * wp, int jd,
                        int Ns, real_t * mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);
  //real_t* mat_inv_dev = (real_t*)dev_ptr(mat_inv);

  int mat_size = ND2 * Ns;
  size_t mat_size2 = mat_size * mat_size * sizeof(real_t);

  real_t* mat_inv_dev;
  CHECK(cudaMalloc((void **)&mat_inv_dev, mat_size2));
  CHECK(cudaMemcpy(mat_inv_dev, mat_inv, mat_size2, cudaMemcpyHostToDevice));

  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst_pad * NC * Ns + blockSize - 1) / blockSize;

  mult_domainwall_5din_ee_inv_dirac_5d_kernel_opt<<<gridSize, blockSize>>>(
                 vp_dev, wp_dev, jd, Ns, mat_inv_dev, Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

  CHECK(cudaFree(mat_inv_dev));
}


#endif
//============================================================END=====

