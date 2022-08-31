
# **Problem Description**

Solve electrophysiology inside a 3D block [1]. The computational domain is a rectangular block, with dimensions of $3×7×20mm$. A small cluster of cells located at one corner of the block (marked as **S**) will receive an initial stimulus, and the ten-Tusscher-Panfilov model is used to model the activation of the rest of the domain.

<p align="center">
   <img src="https://simvascular.github.io/documentation/svfsi/cep/imgs/cuboid.png" width="600">
</p>

We can observe the signal propagation inside the domain.

<p align="center">
   <img src="https://simvascular.github.io/documentation/svfsi/cep/imgs/ttp_cuboid.gif" width="600">
</p>

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template.

## Reference

1. Niederer, S. A., Kerfoot, E., Benson, A. P., Bernabeu, M. O., Bernus, O., Bradley, C., ... & Smith, N. P. (2011). Verification of cardiac tissue electrophysiology simulators using an N-version benchmark. *Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences*, *369*(1954), 4331-4351.
