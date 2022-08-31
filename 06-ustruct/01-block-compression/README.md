
# **Problem Description**

Solve block compression problem using velocity-based solid equation.  The problem set-up is as follows [1]:

<p align="center">
   <img src="../../05-struct/01-block-compression/configuration.png" width="600">
</p>


The final displacement is plotted below.

<p align="center">
   <img src="../../05-struct/01-block-compression/displacement.png" width="600">
</p>

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

When using Taylor-Hood element (P2P1) to solve mixed finite-element problem, remember to turn on the following flag `Use Taylor-Hood type basis: t`.

## Reference

1. Liu, Ju, and Alison L. Marsden.  A Unified Continuum and Variational Multiscale Formulation for Fluids, Solids, and Fluid Structure Interaction.  *Computer Methods in Applied Mechanics and Engineering* 337 (August 2018): 549 97. https://doi.org/10.1016/j.cma.2018.03.045.
