[![DOI](https://zenodo.org/badge/753852599.svg)](https://zenodo.org/doi/10.5281/zenodo.10974974)

# Standardizing eccentricity between astrophysical simulations and gravitational-wave signal modeling 
**Aditya Vijaykumar <sup>1,2</sup>, Alexandra G. Hanselman<sup>3</sup>, Michael Zevin<sup>4,5,6,7</sup>**

<sub>1. Canadian Institute for Theoretical Astrophysics, University of Toronto, 60 St George St,  Toronto, ON M5S 3H8, Canada</sub>  
<sub>2. International Centre for Theoretical Sciences, Tata Institute of Fundamental Research, Bangalore  560089, India</sub>  
<sub>3. Department of Physics, The University of Chicago, 5640 South Ellis Avenue, Chicago, Illinois 60637, USA </sub>  
<sub>4. Kavli Institute for Cosmological Physics, The University of Chicago, 5640 South Ellis Avenue, Chicago, Illinois 60637, USA </sub>  
<sub>5. Enrico Fermi Institute, The University of Chicago, 933 East 56th Street, Chicago, Illinois 60637, USA </sub>  
<sub>6. Adler Planetarium, 1300 South DuSable Lake Shore Drive, Chicago, IL, 60605, USA </sub>  
<sub>7. Center for Interdisciplinary Exploration and Research in Astrophysics (CIERA), Northwestern University, Evanston, IL, 60201, USA </sub>  

## Summary of scripts

- `system_of_eqns.py` contains the 2PN evolution equations for the eccentricity, mean anomaly (not used in the calculation), and velocity.
- `utils.py` contains some basic utility functions
- `example_script.py` takes an example end state of a binary from astrophysical simulations and illustrates how to extract eccentricity at a specified $f_{22}$ or $M f_{22}$. To reiterate, the general algorithm is:

  1. Use the source frame quantities (initial separation, initial eccentricity, masses) to get the $e$ vs $f_{22}$ (source-frame) curve. Interpolate this curve.
  2. If you want to extract $e$ at fixed $f_{22}$ (det-frame), evaluate the interpolant at $f_{22} \times (1 + z)$, where $z$ is the redshift of the source.
  3. If instead you want $e$ at fixed $M f_{22}$, then evaluate the interpolant at $M f_{22} / M^{\rm tot}_{\rm src}$.


## Citation

We encourage use of the scripts here in other works. If you use the material provided here, please cite:

1. The paper https://ui.adsabs.harvard.edu/abs/2024arXiv240207892V/abstract
2. The zenodo entry (with DOI) for this repository https://zenodo.org/doi/10.5281/zenodo.10974974
3. The zenodo entry for the data release https://zenodo.org/doi/10.5281/zenodo.10633004

Please get in touch on email/github if you have any questions or comments!
