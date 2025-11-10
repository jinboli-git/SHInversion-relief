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
- Spherical harmonic coefficients extracted from Earth2014_TBI2014 (Hirt C, Rexer M (2015) Earth2014: 1 arc-min shape, topography, bedrock and ice-sheet models – Available as gridded data and degree-10,800 spherical harmonics. Int J Appl Earth Obs Geoinf 39:103–112. https://doi.org/10.1016/j.jag.2015.03.001)

/data/ECM1_crust_60E_150E_5S_65N.mat 
- The crustal structure data are extracted from ECM1 (Mooney WD, Barrera-Lopez C, Suárez MG, Castelblanco MA (2023) Earth Crustal Model 1 (ECM1): A 1° x 1° Global Seismic and Density Model. Earth Sci Rev 243:104493. https://doi.org/10.1016/j.earscirev.2023.104493)

/data/GOCO06s.gfc
- GOCO06s gravity model download from https://icgem.gfz-potsdam.de/tom_longtime
