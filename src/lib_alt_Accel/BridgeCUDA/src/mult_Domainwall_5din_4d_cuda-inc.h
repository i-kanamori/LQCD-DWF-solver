/*!
      @file    mult_Domainwall_5din_cuda-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-06-12 16:42:53 #$
      @version $LastChangedRevision: 2710 $
*/

#ifndef MULT_DOMAINWALL_5DIN_4D_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_4D_ACC_INCLUDED

//====================================================================
__global__ 
void mult_domainwall_5din_hopb_4D_dirac_kernel(
     real_t *vp, real_t *up, real_t *wp, int Ns, 
     int bc_x, int bc_y, int bc_z, int bc_t,
     int Nx, int Ny, int Nz, int Nt,
     int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
     int flag) 
{
    int Nxy   = Nx   * Ny;
    int Nxyz  = Nxy  * Nz;
    int Nst   = Nxyz * Nt;
    int Nin5  = NVCD * Ns; 
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ real_t sharedMemory[];
    real_t* ut = &sharedMemory[0];

    if (site < Nst) {

        int ix    = site % Nx;
        int iyzt  = site / Nx;
        int ixy   = site % Nxy;
        int iy    = iyzt % Ny;
        int izt   = site / Nxy;
        int iz    = izt  % Nz;
        int it    = izt  / Nz;
        int ixyz  = site % Nxyz;
        int idir;

        for(int is = 0; is < Ns; ++is) {
            real_t vL[NVCD];
            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            // X direction computations
            idir = 0;

            // Positive X direction computation
            if ((ix < Nx - 1) || (do_comm_x == 0)) {
                int ix2 = (ix + 1) % Nx;
                int nei = ix2 + Nx * iyzt;

                real_t bc2 = (ix == Nx - 1) ? bc_x : 1.0;
                load_u(ut, up, site + Nst * idir);
                real_t wt[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
                }
                mult_wilson_xpb(vL, ut, wt);
            }

            // Negative X direction computation
            if ((ix > 0) || (do_comm_x == 0)) {
                int ix2 = (ix - 1 + Nx) % Nx;
                int nei = ix2 + Nx * iyzt;

                real_t bc2 = (ix == 0) ? bc_x : 1.0;
                load_u(ut, up, nei + Nst * idir);
                real_t wt[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
                }
                mult_wilson_xmb(vL, ut, wt);
            }

            // Y direction computations
            idir = 1;

            if ((iy < Ny-1) || (do_comm_y == 0)) {
              int iy2 = (iy + 1) % Ny;
              int nei = ix + Nx * (iy2 + Ny * izt);

              real_t bc2 = (iy == Ny-1) ? bc_y : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_ypb(vL, ut, wt);
            }

            if ((iy > 0) || (do_comm_y == 0)) {
              int iy2 = (iy - 1 + Ny) % Ny;
              int nei = ix + Nx * (iy2 + Ny * izt);

              real_t bc2 = (iy == 0) ? bc_y : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_ymb(vL, ut, wt);
            }

            idir = 2;

            if ((iz < Nz-1) || (do_comm_z == 0)) {
              int iz2 = (iz + 1) % Nz;
              int nei = ixy + Nxy * (iz2 + Nz * it);

              real_t bc2 = (iz == Nz-1) ? bc_z : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
              mult_wilson_zpb(vL, ut, wt);
            }

            if ((iz > 0) || (do_comm_z == 0)) {
              int iz2 = (iz - 1 + Nz) % Nz;
              int nei = ixy + Nxy * (iz2 + Nz * it);

              real_t bc2 = (iz == 0) ? bc_z : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
              mult_wilson_zmb(vL, ut, wt);
            }

            idir = 3;

            if ((it < Nt-1) || (do_comm_t == 0)) {
              int it2 = (it + 1) % Nt;
              int nei = ixyz + Nxyz * it2;

              real_t bc2 = (it == Nt-1) ? bc_t : 1.0;
              load_u(ut, up, site + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                  wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
              mult_wilson_tpb_dirac(vL, ut, wt);
            }

            if ((it > 0) || (do_comm_t == 0)) {
              int it2 = (it - 1 + Nt) % Nt;
              int nei = ixyz + Nxyz * it2;

              real_t bc2 = (it == 0) ? bc_t : 1.0;
              load_u(ut, up, nei + Nst * idir);
              real_t wt[NVCD];

              for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	          }
                mult_wilson_tmb_dirac(vL, ut, wt);
              }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_hopb_dirac_4d(
     real_t *vp, real_t *up, real_t *wp, int Ns, int *bc, int *Nsize,
     int *do_comm, int flag) {

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

    mult_domainwall_5din_hopb_4D_dirac_kernel<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns,
    bc_x, bc_y, bc_z, bc_t,
    Nx, Ny, Nz, Nt,
    do_comm_x, do_comm_y, do_comm_z, do_comm_t,
    flag); 
    
  CHECK(cudaDeviceSynchronize());

}


#endif
//============================================================END=====
