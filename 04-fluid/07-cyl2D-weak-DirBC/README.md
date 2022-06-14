
# **Problem Description**

Solve the classic 2D flow past circular cylinder problem.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

```
   Add BC: cyl_wall {
      Type: Dir
      Time dependence: Steady
      Value: 0.0
      Weakly applied: t
      Penalty parameter: 1000.0
   }
```

Here, instead of strongly enforcing the no-slip boundary condition on the cylinder wall (i.e., prescribing the zero velocity at the nodal values), we use a penalty function to incorporate this boundary condition into the weak form (i.e., weakly apply the Dirichlet boundary condition.) This practice can help with convergence in some complicated problems. 
