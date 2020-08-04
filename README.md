<hr> 

<h1><p align="center"> Vertex centred explicit solid dynamics toolkit for OpenFOAM
</p></h1>



<p align="center"> 
  <a href="https://travis-ci.org/jibranhaider/explicitSolidDynamics_vc" target="_blank">
    <img alt="Travis (.org)" src="https://img.shields.io/travis/jibranhaider/explicitsoliddynamics_vc"> &nbsp;
  </a>  
  <img alt="OpenFOAM 6" src="https://img.shields.io/badge/OpenFOAM-v7_| v6_| v5-darkgreen.svg"> &nbsp;
  <a href="https://github.com/jibranhaider/explicitSolidDynamics/blob/master/LICENSE">
    <img alt="GPLv3 license" src="https://img.shields.io/badge/License-GPLv3-orange.svg"> &nbsp;
  </a> 
  <a href="https://openfoam.org/dev/coding-style-guide">
  <img alt="Code style" src="https://img.shields.io/badge/Coding_style-OpenFOAM-yellow.svg"> 
  </a>    
</p>

<p align="center">   
  <img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/jibranhaider/explicitsoliddynamics_vc"> &nbsp;
  <img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/jibranhaider/explicitsoliddynamics_vc?color=lightgrey"> &nbsp;
  <img alt="GitHub repo size in bytes" src="https://img.shields.io/github/repo-size/jibranhaider/explicitsoliddynamics_vc?label=Repository%20size&color=lightgrey"> &nbsp;
  <a href="https://jibranhaider.com/research/explicit-solid-dynamics-in-openfoam/">
    <img alt="Website" src="https://img.shields.io/website-up-down-green-red/https/jibranhaider.weebly.com/research.html.svg?label=Website"> 
  </a>  
</p>


<br/>

## 1. Introduction    

<p align="justify">
This toolkit is based on a vertex centred Finite Volume Method to predict large deformation in solids. The governing equations comprise of a system of first order hyperbolic conservation laws for linear momentum and deformation gradient tensor similar to the ones found in Computational Fluid Dynamics.
</p> 

<br/>

#### Preprocessing
<img src="docs/preprocessing/median_dual.png" width="100%">

<br/>

#### Results
<img src="/docs/results/twistingColumn/0.png" width="12%"> &nbsp; &nbsp; &nbsp; &nbsp;
<img src="/docs/results/twistingColumn/2.png" width="12%"> &nbsp; &nbsp; &nbsp; &nbsp;
<img src="/docs/results/twistingColumn/4.png" width="12%"> &nbsp; &nbsp; &nbsp; &nbsp;
<img src="/docs/results/twistingColumn/6.png" width="12%"> &nbsp; &nbsp; &nbsp; &nbsp;
<img src="/docs/results/twistingColumn/8.png" width="12%"> &nbsp; &nbsp; &nbsp; &nbsp;
<img src="/docs/results/twistingColumn/10.png" width="12%">

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
* Quasi-static formulation.
* Parallel implementation of the vertex centred scheme.
* Multiple body contact problems.


<br/>
<hr> 

## 4. Author
This toolkit is developed and maintained by [Jibran Haider](http://jibranhaider.weebly.com/). The following individuals are acknowledged for their support:
* [Dr. Chun Hean Lee](https://www.gla.ac.uk/schools/engineering/staff/chunheanlee/)
* [Prof. Antonio J. Gil](https://www.swansea.ac.uk/staff/engineering/a.j.gil/)
* [Prof. Javier Bonet](https://www.researchgate.net/profile/Javier_Bonet)


<br/>
<hr> 

## 5. License
This toolkit is released under the GNU General Public License (version 3). More details can be found in the [LICENSE](LICENSE) file. 


<br/>
<hr> 

## 6. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.

