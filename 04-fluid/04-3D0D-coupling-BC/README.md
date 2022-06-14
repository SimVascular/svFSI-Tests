
# **Problem Description**

Solve the same problem as in [04-fluid/01-pipe3D_RCR](../01-pipe3D_RCR). Instead of using the RCR within the `svFSI` solver, this example demonstrates how to set up RCR boundary condition in a more generalized framework using genBC or cplBC.

## Introduction

Both genBC and cplBC provide a framework to programmatically define custom inflow and outflow boundary conditions for a CFD simulation. The framework allows users to create an arbitrary lumped parameter network (LPN, or 0D model) layout suitable for their application. Some common examples include a lumped parameter heart model that models contraction of the heart chambers to use as an inlet boundary condition, sophisticated models of the downstream circulation for various areas of the body such as the legs and upper body, or a closed-loop formulation where all outflow of the SimVascular model returns back to the inflow after passing through the veins, heart, and pulmonary arteries.

**Essentially, genBC and cplBC are two different implementations of the same functionality, and cplBC is the preferred choice.**  genBC is a legacy code developed for [svSolver](https://github.com/SimVascular/svSolver), and cplBC is developed specifically for `svFSI`. cplBC takes advantage of the data structures in `svFSI`, and is more user-friendly. Still, `svFSI` provides backward compatibility for genBC so that svSolver users can migrate to the new solver easily. 

## Configuration of genBC

There are excellent tutorials online that show the users how to set up genBC step-by-step.
SimVascular website: https://simvascular.github.io/docsGenBC.html
Youtube tutorial: https://www.youtube.com/watch?v=znfV0XLV79s&ab_channel=SimVascular

**We strongly encourage users to go through these tutorials first to become familiar with the concept and the workflow.** The set-up of cplBC is mostly the same with some tweaks.

The input file [svFSI_genBC.inp](./svFSI_genBC.inp) follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

```
   Couple to genBC: SI {
      0D code file path: genBC/genBC.exe
   }
```

This tells the solver that the 0d models will be calculated through genBC. Options to couple 0D codes with svFSI are `N`: none; `I`: implicit; `SI`: semi-implicit; `E`: explicit.

```
   Add BC: lumen_inlet {
      Type: Dir
      Time dependence: Coupled
      Impose flux: t
      Zero out perimeter: t
   }

   Add BC: lumen_outlet {
      Type: Neu
      Time dependence: Coupled
   }
```

Following the online tutorials, we use LPN for both inlet Dirichlet boundary condition and outlet RCR boundary condition. The state variables of the LPN usually include flow rate instead of velocity. Hence, it is important to set `Impose flux` to true.

## Configuration of cplBC

Configuration of cplBC is very similar to genBC, and the following files require user's attention: [svFSI_cplBC.inp](./svFSI_cplBC.inp), [cplBC/foo.ini](./cplBC/foo.ini) and [cplBC/src/USER.f](./cplBC/src/USER.f).

### svFSI Input

The input file for cplBC is more instructive than that of genBC. Many of the same parameters in [genBC/src/USER.f](./genBC/src/USER.f) are moved to `svFSI` input file, so that they can be easily changed without the need to recompile cplBC.

```
   Couple to cplBC: SI {
      Number of unknowns: 2
      0D code file path: cplBC/cplBC.exe
      Unknowns initialization file path: cplBC/foo.ini
      File name for 0D-3D communication: CPLBC_0D_3D.tmp
      File name for saving unknowns: cplBC_AllData
      Number of user-defined outputs: 2
   }
```

`Number of unknowns` represents the total number of unknowns within the 0D system. It can be larger than the number of boundaries. `Unknowns initialization file path` specifies from which file to read the initial condition, and requires modification for your specific case. `Number of user-defined outputs` tells the solver to output certain amount of additional variables, which are defined in [cplBC/src/USER.f](./cplBC/src/USER.f).

### cplBC USER.f

[USER.f](./cplBC/src/USER.f) in cplBC is greatly simplified than that in genBC. Assuming the users are already familiar with the configuration of genBC, we list locations in USER.f that are different.

```fortran
!--------------------------------------------------------------------
!     List of all coupled faces, BC groups and Xptr is specified here
!     face%bGrp : either cplBC_Dir or cplBC_Neu
!     face%name: the name of the face in svFSI input file
!     face%Xptr: the corresponding index of X for this face
!--------------------------------------------------------------------
      face(1)%bGrp  = cplBC_Dir
      face(2)%bGrp  = cplBC_Neu
      face(1)%name  = "lumen_inlet"
      face(2)%name  = "lumen_outlet"
      face(1)%Xptr  = 1
      face(2)%Xptr  = 2
```

The number of faces should equal the number of boundaries that are set to `Time dependence: Coupled` in  svFSI input. `face%bGrp` indicates the boundary condition on this face is either Dirichlet or Neumann. `face%name` should match the boundary name in svFSI input. `face%Xptr` matches the internal variable in cplBC to the face. In this case, `X(1)` is the flow rate on face `lumen_inlet` and `X(2)` is the pressure on face `lumen_outlet`.

```fortran
!     BC
      Rp = 121D0
      C  = 1.5D-4
      Rd = 1212D0

      f(1) = -40D0 * pi * pi * SIN(2D0 * pi * t) 
      f(2)= (1D0/C) * (Q(2) - x(2)/Rd)
      offst(2) = Q(2)*Rp

!     Assign the additional parameters to be printed
      Xw(1)=offst(2)
```

Here, we specify a flow rate with cosine profile as the inlet boundary condition, so `f(1)=dX(1)/dt` is set to a sine wave.  RCR boundary is used at outlet, and the corresponding formulation is defined in `f(2)`. Please note that this implementation is different from genBC. You can find more on this bellow.

## Differences between genBC and cplBC

Though genBC and cplBC have the same functionality, cplBC is considerably easier to use in case of complicated [patient-specific cases](./patient-specific-case). Here, we highlight some major differences between these two.

### Matching faces between 3D and 0D

In genBC, user has to provide face tags in the USER.f file and make sure that these tags match the ones in solver.inp file for svSolver. However, cplBC uses face names to identify these tags automatically. In essence, when using genBC, it is the user's responsibility to match both the bcs provided in the svPresolver and svSolver files, and also in the USER.f file. However, for cplBC, user provides this information in the input file (svFSI.inp) in a simple way.

### Implementation of 0D model

The 0D model or LPN is defined in USER.f in both cases, but the implementations are not compatible. Here are the same boundary conditions implemented in genBC

```fortran
      Rp = 121D0
      C  = 1.5D-4
      Rd = 1212D0

      f(1) = -40D0 * pi * pi * SIN(2D0 * pi * t) 
      f(2)= (1D0/C) * (Q(1) - x(2)/Rd)
      offset(2) = Q(1)*Rp
```

Note that `Q(1)` is used here, while `Q(2)` is used in the cplBC case above. In genBC, `Q(:)` and `P(:)` are defined as flow rates of Neumann faces and pressure of Dirichlet faces, respectively. It is the user's responsibility to carefully match the component `Q(i)` to the corresponding face, which can be error-prone when the number of faces are large. On the other hand, in cplBC, `Q(i)` simply represents the flow rate on the ith face.

### Outputs from 0D model

genBC writes out the results in `AllData` in the current folder, and it includes all the unknowns and the user-specified outputs. In this case, the order will be

```
inlet_flux  outlet_pressure  time  outlet_flux  offset(2)
```

Here, the last three are user-defined outputs.

On the other hand, results from cplBC are generated within `svFSI` `SUBROUTINE TXT` and stored in `cplBC_AllData` in the results folder. By default, it will not only output the unknowns, but also flow rate and pressure on all faces.  Also, `svFSI` will automatically write the first user-defined output in the first column. Therefore, **it is strongly recommended to put time as the first user specified output**. In this case, the results will be

```
time  inlet_flux  outlet_flux  inlet_pressure  outlet_pressure  offst(2)
```

Here, the first and last one are the user defined outputs.

### Initial conditions

In genBC, the initial conditions are specified in USER.f through variable `tZeroX`. Hence, user needs to recompile genBC every time it changes. In cplBC, the initial conditions are provided through an input file and is more convenient.

## Advanced Usage
This folder also contains a [patient-specific case](./patient-specific-closed-loop-case), in which more advanced usage of 3D-0D coupling is demonstrated.

