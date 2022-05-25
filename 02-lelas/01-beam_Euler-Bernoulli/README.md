
# **Problem Description**

Solve Eulear-Bernoulli beam problem with given traction boundary condition.

This example mainly demonstrates the element types supported by `svFSI`. Several types of tetrahedral and hexahedral elements are included here. Many other types are available for different problems. The complete list of element types supported by `svFSI` is: 
```
Point;
Line: linear, quadratic;
Triangle: linear, quadratic;
Quads: bilinear; serendipity, biquadratic;
Tetrahedron: linear, quadratic;
Hex: trilinear, quadratic/serendipity, triquadratic;
Wedge;
NURBS.
```

Users can generate their own mesh using Gambit or Gmsh, and convert it into `svFSI` compatible format through [mesh_converter](https://github.com/SimVascular/svFSI-Tools).
