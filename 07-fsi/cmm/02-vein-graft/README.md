
# **Problem Description**

Perform fluid-structure interaction (FSI) simulation in a stenosed vein graft using coupled momentum method (CMM) and variable wall properties. Unsteady Dirichlet velocities are provided at the inflow and a RCR boundary condition is provided at the outflow.

## Solution workflow

This problem is solved in three major steps:

1. Solve for fluid flow alone assuming rigid walls and using the same flow boundary conditions as the FSI simulation.

2. Inflate the solid domain under a diastolic pressure loading. If the FSI problem involves variable wall properties, the same should be used for the inflation step as well.

3. Perform CMM-type FSI simulation by initializing the flow field from Step 1 and displacements from Step 2.

The input files for the above three steps follow the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Flow BCs (Steps 1 & 3)

In this example, Dirichlet BCs are applied at the inflow by reading `bct.vtp` file generated from SimVascular. It is required that the `Time dependence` is `General` for loading a bct.vtp file as given below:

`Add BC: inflow {
    Type: Dir
    Time dependence: General
    BCT file path: bct.vtp
 }`

RCR BCs are applied at the outflow and details are provided in the example, [`04-fluid/01-pipe3D_RCR`](https://github.com/SimVascular/svFSI-Tests/blob/master/04-fluid/01-pipe3D_RCR/README.md)

## Inflation (Step 2)

Four settings in the input file need closer attention when performing inflation step for CMM:

### (a) Load mesh as shell surface

For CMM-based FSI simulation, vessel wall is modeled as a shell surface and therefore, the inflation step is also performed on the wall set as a shell,

`Add mesh: wall {
    Set mesh as shell: t
    Mesh file path: ../mesh/walls_combined.vtp
 }`

### (b) Set Initialize keyword

Within `Add equation: CMM`, the keyword `Initialize` needs to be set to **inflate** as,

`Add equation: CMM {
    ...
    Initialize: inflate
    ...
 }`

Alternately, one could also prestress the vessel wall as,

`Add equation: CMM {
    ...
    Initialize: prestress
    ...
 }`

### (c) Diastolic loading

CMM is initialized via inflation or prestress approaches by applying a diastolic load on the wall. This is performed either using a constant pressure or traction, or from a spatially-varying pressure or traction computed from the rigid wall flow simulation (**Step 1**). Examples are provided below for constant pressure,

`Add BF: wall {
    Type: Neumann
    Value: 8000.0
 }`

and, spatially-varying traction field as,

`Add BF: wall {
    Type: Traction
    Time dependence: Spatial
    Spatial values file path: rigid_wall_mean_traction.vtp
 }`

Note that because the wall is modeled as a shell surface, the above setting applies pressure/traction as a *body force* on the shell surface using `Add BF` keyword and **not** as a *surface traction* using the `Add BC` keyword.

### (d) Variable wall properties

CMM-based FSI simulation allows loading variable wall properties using the keyword,

`Variable wall properties: wall {
    Wall properties file path: svFSI_var_wall_props.vtp
 }`

`svFSI` expects the provided vtp file to contain two nodally-stored scalar variables, `Elasticity_modulus` and `Thickness`, indicating the Young's modulus and wall thickness of the shell surface, respectively.

## FSI using CMM (Step 3)

For performing FSI simulation using CMM method in `svFSI`, we have to set `Add equation: CMM` and use results from the previous steps (**Steps 1 & 2**). Two places where we load the results from the previous steps are,

### (a) Initialize flow field

For faster convergence, we initialize the flow for the FSI simulation using data from **Step 1** within `Add mesh` keyword as,

`Add mesh: msh {
    ...
    Initial velocities file path: init_flow_rigid_3000.vtu
    Initial pressures file path:  init_flow_rigid_3000.vtu
    ...
 }`

In the above statements, the solution from the rigid wall simulation (**Step 1**) at time step 3000 is used to initialize the velocity and pressure fields.

### (b) CMM BC and initial displacements

Initial displacements from the inflation step (**Step 2**) are loaded on the wall set to be a CMM type boundary condition as,

`Add BC: wall {
    Type: CMM
    Initial displacements file path: svFSI_wall_disps_inflate_50.vtu
 }`

Note that the BC type is set to be CMM for the wall (`Add BC: wall {Type: CMM}`). Further, wall displacements from **Step 2** are loaded at the 50th time step using the keyword `Initial displacements file path`.

Alternately, wall prestress could be provided as input if the initialization in **Step 2** is based on prestressing. An example command is shown below,

`Add BC: wall {
    Type: CMM
    Prestress file path: svFSI_init_wall_prestress.vtu
 }`

Additionally, for the CMM-based FSI simulation, the densities of the fluid and solid domains have to be specified separately as,

`Add equation: CMM {
    ...
    Fluid density: 1.06
    Solid density: 1.0
    ...
}`

It is extremely important that all the boundary conditions and material parameters for the CMM-based FSI (**Step 3**) are same as the ones used in Steps 1 and 2.
