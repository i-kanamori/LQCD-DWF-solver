# Additional Information for OpenACC version

see `README.md` for full description of this kernel application.
This note is to provide information of `5din` implementation of the matrix-vector
multiplication kernel.

## data structure

The `5din` is a native implementation of Domainwall operator,
which does not internally call Wilson operator.  The naming comes from
the data layout.  One of the indices of the field is moved from the most outer (as if it is 5th coordinate in addition to 4-dim structured lattice) to the inner.

The indexing of (linearlized) site index and on-site index are defined in 
  src/lib_alt_Accel/BridgeOpenACC/define_index.h
`through ID2_*` macros.  It assume the warp length `NWP`, which is define in
  src/lib_alt_Accel/BridgeOpenACC/define_params.h

The field (i.e., vector) is one-dimensional array but intended the following data structure
```
  site = site_inner + NMP * site_outer
  v[ex][site_outer][is][id][ic][re/im][site_inner]
  ex : not used for the vectors
  is = 0...Ns-1  (typically Ns=8...16; the test case is Ns=8)
  id = 0...3     (Dirac spinor)
  ic = 0...2     (color)
  re/im: real and imaginary part of complex numbers
  (ic and re/im makes ivc = 0...5; id, ic and re/im make ivcd = 0...23)
```

## main stencil kernels

The main kernel for sparse matrix-vector multiplication is
`D_eo` : `src/lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 792`
and
`Ddag_eo` : `src/lib_alt_Acce"#l/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 883`

Among the sub kernels called in these two, `BridgeACC::mult_domainwall_5din_eo_hopb_dirac` is the most important one.

`D_eo` has 4 parts:
  * `BridgeACC::mult_domainwall_5din_eo_5dir_dirac` 
      * given at `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 87`
      * parallel d.o.f: site and ivc (ic and re/im)

  * `BridgeACC::mult_domainwall_5din_eo_hop1_dirac`
      * this is to pack the boundary data for communication
      * given at `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 1420`
      * parallel d.o.f: surface part of the site.
      * is can be also easily parallelized.

  * `BridgeACC::mult_domainwall_5din_eo_hopb_dirac`
      * the bulk part of the stencil computation. _The most important_for the performance on single node.  It contains hopping from 8 directions.
      * in addition the implementation called in the current version, there are several trial implementations.
      * the current one: given at `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 1003`, using site and is for the parallelization
      * trial: `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 880`, using site and 8 hopping directions for parallelization
      * trial: `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 1117`, using site for parallelization
      * trial: `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 1219`, using site for parallelization

  * `BridgeACC::mult_domainwall_5din_eo_hop2_dirac`
      * stencil computation of the boundary.  It can be severe in the load balance.
      * implemented at `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 1658`, using site for parallelization


`Ddag_eo` has 5 parts, and 2 are new:

  * `BridgeACC::mult_domainwall_5din_mult_gm5_dirac`
     * permutation in id space
     * given at `src/lib_alt_Accel/BridgeOpenACC/src/mult_Domainwall_5din_openacc-inc.h: 293`
     * site is used for parallelization, but others d.o,f can be also used.
     * merging to `BridgeACC::mult_domainwall_5din_eo_hop1_dirac` and `BridgeACC::mult_domainwall_5din_eo_hopb_dirac` can be one possibility to improve data locality.

  * `BridgeACC::mult_domainwall_5din_eo_5dirdag_dirac`
     * implemented at `src/lib_alt_Accel/BridgeACC/src/mult_Domainwall_5din_eo_openacc-inc.h: 347`, using site, ic and re/im for parallelization


## other kernels

The followings do not contain hopping. They are also diagonal in ic and re/im space.

  * `L_inv` : `src/lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 1049`, calls `BridgeACC::mult_domainwall_5din_L_inv_dirac` at `src/lib_alt_Accel/BridgeOpenACC/src/mult_Domainwall_5din_openacc-inc.h: 1085`

  * `U_inv` : `src/lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 1067`, calls `BridgeACC::mult_domainwall_5din_U_inv_dirac` at `src/lib_alt_Accel/BridgeOpenACC/src/mult_Domainwall_5din_openacc-inc.h: 1191`

  * `Ldag_inv` : `src/lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 1105`, calls `BridgeACC::mult_domainwall_5din_Ldag_inv_dirac` at `src/lib_alt_Accel/BridgeOpenACC/src/mult_Domainwall_5din_openacc-inc.h: 1278`

  * `Udag_inv` : `src/lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo-tmpl.h: 1086`, calls `BridgeACC::mult_domainwall_5din_Udag_inv_dirac` at `src/lib_alt_Accel/BridgeOpenACC/src/mult_Domainwall_5din_openacc-inc.h: 1356`

They all use site d.o.f for parallelization but ic and re/im can be also parallelized.  Merging `L_inv`and `U_inv`, `Ldag_inv` and `Udag_inv` may also improved the performance.
