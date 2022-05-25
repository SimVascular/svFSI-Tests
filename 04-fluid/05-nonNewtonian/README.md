
# **Problem Description**

Solve the 3D pipe flow problem with different viscosity models. Currently, `svFSI` supports three viscosity model: Newtonian, Carreau-Yasuda and Casson.

Carreau-Yasuda viscosity model is defined as

![my equation](https://latex.codecogs.com/svg.image?\eta=\eta_\infty&plus;(\eta_0-\eta_\infty)\left[&space;1&plus;\left(&space;\lambda&space;\dot{\gamma}&space;\right)^a&space;\right]^{\frac{n-1}{a}}&space;)

Here

![myequation](https://latex.codecogs.com/svg.image?\eta_\infty:&space;\text{Limiting&space;high&space;shear-rate&space;viscosity}\\&space;\eta_0:&space;\text{Limiting&space;low&space;shear-rate&space;viscosity}&space;\\&space;\lambda:&space;\text{Shear-rate&space;tensor&space;multiplier}&space;\\\dot{\gamma}:&space;\text{shear&space;rate}&space;\\a:&space;\text{Shear-rate&space;tensor&space;exponent}&space;\\n:&space;&space;\text{Power-law&space;index})





The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

For Newtonian fluid:

```
   Viscosity: Constant {
      Vsalue: 0.04
   }
```

For Casson fluid

```
   Viscosity: Cassons {
      Asymptotic viscosity parameter: 0.3953
      Yield stress parameter: 0.22803
      Low shear-rate threshold: 0.5
   }
```

For Carreau-Yasuda fluid

```
   Viscosity: Carreau-Yasuda {
      Limiting high shear-rate viscosity: 0.022
      Limiting low shear-rate viscosity: 0.22
      Shear-rate tensor multiplier (lamda): 0.11
      Shear-rate tensor exponent (a): 0.644
      Power-law index (n): 0.392
   }
```



