<hr> 

<h1><p align="center"> Vertex centred code for solids in OpenFOAM
</p></h1>


<br/>

## 1. Introduction    

<p align="justify">
This toolkit is based on a vertex centred Finite Volume Method to predict large deformation in solids. The governing equations comprise of first order hyperbolic conservation laws for linear momentum and deformation gradient tensor of the system. This helps to bridge the gap between Computational Fluid Dynamics and Computational Solid Dynamics. The governing equations are spatially discretised using a low order vetex centred Finite Volume Method along with an upwind Riemann solver. 
</p> 


<br/>
<hr> 

## 2. How to use this toolkit?

The structure of this repository is similar to the [explicitSolidDynamics](https://github.com/jibranhaider/explicitSolidDynamics) toolkit which is based on the cell centred Finite Volume Method. Therefore, the reader is advised to consult the following webpages and adapt accordingly: 
* [Installation](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Installation)
* [Tutorials](https://github.com/jibranhaider/explicitSolidDynamics/wiki/Tutorials)


<br/>
<hr> 

## 3. FAQs

#### Does this toolkit support polyhedral elements?
* At the moment only tetrahedral elements are supported.

#### What additional features can we expect in future releases?
* Advanced constitutive models including plasticity, thermo/visco-elasticity.
* Multiple body contact problems.


<br/>
<hr> 

## 4. Author
This toolkit is developed and maintained by [Jibran Haider](http://jibranhaider.weebly.com/) (Swansea University). 


<br/>
<hr> 

## 5. License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file. 


<br/>
<hr> 

## 6. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
