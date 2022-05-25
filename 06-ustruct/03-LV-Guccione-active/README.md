
# **Problem Description**

Solve passive inflation of an idealized left-ventricular model. The problem set-up can be found in Problem 3 of the following publication:

> Land, Sander, Viatcheslav Gurev, Sander Arens, Christoph M. Augustin, Lukas Baron, Robert Blake, Chris Bradley, et al.  Verification of Cardiac Mechanics Software: Benchmark Problems and Solutions for Testing Active and Passive Material Behaviour.  *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences* 471, no. 2184 (December 2015): 20150641. https://doi.org/10.1098/rspa.2015.0641.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Fiber Direction

In this problem, the left ventricle is modeled as the Guccione hyperelastic material. It requires two user-supplied fiber directions. The longitudinal and sheetlet fiber directions are specified using the following directives:

```
   Fiber direction file path: mesh/P1/fibersLong.vtu
   Fiber direction file path: mesh/P1/fibersSheet.vtu
```

The spatial distribution of the fibers are stored in the vector variable `FIB_DIR`.

## Active Cardiac Contraction

The active contraction of the idealized LV is modeled through active stress, which is called fiber reinforcement stress in `svFSI`. 

```
   Fiber reinforcement stress: Unsteady {
      Temporal values file path: fib_stress.dat
      Ramp function: t
   }
```

More on the ramp function, please refer to [01-heat/01-diffusion-line-source](../../01-heat/01-diffusion-line-source/README.md).
