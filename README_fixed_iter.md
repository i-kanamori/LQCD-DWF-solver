## Parameters for fixed number of iterations

If the target tolerance is too small (or the maximum iteration count is too small), the solver works as a fixed number iteration solver.  In `test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml` file, in the `Solver:` section, one can set the iteration counts and the target tolerance.

Here is an example:

```
Solver:
  ## solver_type = {CGNE,CGNR,BiCGStab_Cmplx,BiCGStab_L_Cmplx,BiCGStab_DS_L_Cmplx,BiCGStab_IDS_L_Cmplx,GMRES_m_Cmplx}
  solver_type                         : CGNR
  maximum_number_of_iteration_single  : 100
  maximum_number_of_iteration         : 10
  maximum_number_of_restart           : 1
  convergence_criterion_squared       : 0
  convergence_criterion_squared_single: 0
```

* single solver never "converges" and always hits the maximum iterations (`convergence_criterion_squared_single` is 0).  If `convergence_criterion_squared_single` is not given, the default is 1e-8.

* outer double precision iteration never "converges" and always hits the maximum iterations: (`convergence_criterion_squared` is 0)

* iteration count for single solver: fixed to 100  (= `maximum_number_of_iteration_single`).  If `maximum_number_of_iteration_single` is not given, the same value for hte outer iteration is used.

* iteration count for outer iteration: fixed to 10  (= `maximum_number_of_iteration` x `maximum_number_of_restart`)

The log file should be something like as follow:

```
Error at ASolver_CG<Afield<float,QXS> >: not converged
  iter(final):      100   2.302062362957713e-08
APrecond_Mixedprec: nconv = -1  diff = 2.302068e-08
Error at ASolver_CG<Afield<float,QXS> >: not converged
  iter(final):      100   6.037626415491104e-03
...
Error at ASolver_CG<Afield<float,QXS> >: not converged
  iter(final):      100   5.506889522075653e-02
APrecond_Mixedprec: nconv = -1  diff = 5.506889e-02
Error at ASolver_CG<Afield<float,QXS> >: not converged
  iter(final):      100   2.378175407648087e-02
APrecond_Mixedprec: nconv = -1  diff = 2.378175e-02
Richardson NOT converged: iter=10
Richardson  nconv = -1, diff2 =  1.349084086950851e-22
```

* `APrecond_Mixedprec:` information of the single precision solver and `nconv` is the iteration counts.  If not converged (i.e., fixed iteration), `nconv = -1` and  `iter(final): ` refers the iteration counts instaed.

* `Richardson` referes the outer iterations and `nconv` is the iteration counts. Again, `nconv = -1` is for not converged and `iter=10` refers the maximum iterations.
