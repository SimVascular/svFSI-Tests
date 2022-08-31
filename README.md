
# **Examples for svFSI**

This repository contains examples for using and testing `svFSI`. Examples are organized based on the physics solved. Each example folder has multiple cases with several options that include boundary conditions, linear solver settings, etc. User are encouraged to review [`svFSI_master.inp`](./svFSI_master.inp) for a detailed description of the input file and its settings.

A description of each example is provided within the folder. Examples include:

1. heat: diffusion equation solver
2. lelas: linear elasticity
3. stokes: Stokes flow
4. fluid: incompressible flow (Navier-Stokes equation)
   1. **pipe3D_RCR**: pipe flow with RCR boundary condition
   2. dye_AD: dye transportation with fluid flow
   3. driven-cavity-2D: 2D driven cavity flow
   4. **3D0D-coupling-BC**: generic boundary condition that couples 0D model to 3D flow
   5. nonNewtonian: Carreau-Yasuda and Casson non-Newtonian flow
   6. **channel-flow-2D**: 2D channel flow
   7. **cyl2D-weak_DirBC**: 2D flow past circular cylinder with weakly enforced Dirichlet BC

5. struct: nonlinear displacement-based solid mechanics
   1. **block-compression**: block compression 
   2. **LV-Guccione-passive**: passive inflation of the left ventricular model using Guccion material model

6. ustruct: nonlinear velocity-pressure based solid mechanics
   1. **block-compression**: block compression
   2. **tensile-adventitia_HGO**: simple tension of arterial strip
   3. **LV-Guccione-active**: simulate active contraction of the left ventricular model using active stress model

7. fsi: fluid-structure interaction (ALE, CMM)
   1. ale: examples using the Arbitrary Lagrangian-Eulerian Method
      1. channel-leaflets_2D: 2D heart valve model
      2. **channel-block-flag_2D**: 2D flag behind a block
      3. pipe_3D: propagation of pressure pulse inside aorta
      4. **pipe_prestress**: model FSI in arterial system with prestress

   2. cmm: examples using the Coupled Momentum Method
      1. **pipe_RCR**: flow in elastic pipe with RCR boundary condition
      2. **vein-graft**: FSI simulation in a stenosed vein graft

8. cep: cardiac electrophysiology modeling
   1. 2Dsqr_AP: simulate electrophysiology inside a 2D plane using the Aliev-Panfilov model
   2. 1Dcable_tTP: simulate electrophysiology inside a 1D cable using the ten-Tusscher-Panfilov model
   3. **benchmark_tTP**: simulate electrophysiology inside a 3D block using the ten-Tusscher-Panfilov model
   4. **2Dspiral_BO**: simulate spiral wave in a plate using the Bueno-Orovio-Cherry-Fenton model 
   5. **Purkinje**: simulate signal propagation inside the Purkinje network using the ten-Tusscher-Panfilov model 

9. shells: shell structures. (This part of svFSI is still under active development.)
   1. **plate**: deformation of a 2D plate
   2. **valve**: opening and closing of aortic valve


Note: **Highlighted** cases contain simulated results which users can use to make sure `svFSI` is installed correctly. 
