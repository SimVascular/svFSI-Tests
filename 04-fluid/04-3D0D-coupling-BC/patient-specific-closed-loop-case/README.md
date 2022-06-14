# **Problem Description**

Solve a patient-specific case of left pulmonary artery stenosis using 3D-0D coupling. The methodology is described in Ref [1]. The diagram of the 0D model is depicted below.

<p align="center">
   <img src="https://ars.els-cdn.com/content/image/1-s2.0-S0022522314021849-gr1bc_lrg.jpg" width="1000">
</p>

Parameters of Patient 2 is used here (see [Table E1](https://www.sciencedirect.com/science/article/pii/S0022522314021849#tblE1)). In terms of implementation, the parameters are stored in [parameters.f](./cplBC/include/parameters.f).

## Reference

1. Schiavazzi, Daniele E., Ethan O. Kung, Alison L. Marsden, Catriona Baker, Giancarlo Pennati, Tain-Yen Hsia, Anthony Hlavacek, and Adam L. Dorfman. “Hemodynamic Effects of Left Pulmonary Artery Stenosis after Superior Cavopulmonary Connection: A Patient-Specific Multiscale Modeling Study.” *The Journal of Thoracic and Cardiovascular Surgery* 149, no. 3 (March 2015): 689-696.e3. https://doi.org/10.1016/j.jtcvs.2014.12.040.
