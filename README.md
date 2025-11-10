# SHinversion-relief

Code accompanying the manuscript:

> Li, J. et al., *A spherical harmonic method for modeling interface topography with 3D variable density contrasts from potential coefficients*, submitted to **Computers & Geosciences**.

This repository provides a MATLAB implementation of the spherical harmonic method
for modeling subsurface density-interface topography from gravitational potential
coefficients, allowing for fully 3D variable density contrasts at regional to
global scales.

## Requirements
- MATLAB R2023a or later (earlier versions may also work).
## Repository structure
```text
/README.md                 – This file
/LICENSE                   – License information
/m_code/                   – MATLAB codes
/m_code/shell.m              – Example: spherical shell test 
/m_code/seafloor.m           – Example: global bathymetry test
/m_code/Mohoinv.m            – Example: regional Moho test
/m_code/H2Qnm.m              – Core forward modeling routine
/m_code/Qnm2H.m              – Core inversion routine
/data/                     – The data used in the codes
/data/Earth2014_truncate_360degree.mat
- Spherical harmonic coefficients extracted from Earth2014_TBI2014
/data/ECM1_crust_60E_150E_5S_65N.mat 
- The crustal structure data are extracted from ECM1
/data/GOCO06s.gfc
- GOCO06s gravity model download from https://icgem.gfz-potsdam.de/tom_longtime
