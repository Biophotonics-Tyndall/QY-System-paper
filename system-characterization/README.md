# Quantum Yield Analysis

**Author**: Jean Matias \
**Email**: jean.matias@tyndall.ie 

This code is Matlab (v. R2019b) based used to analyse data acquired with the QY setup built in Tyndall National Institute. 

## How to run

All the equations and data treatment are made with the classes in the *./loaders/* folder. The scripts in the root path are simplified scripts that uses the main classes.  

To generate the final figures (Fig. 3 and Fig. 4) the following scripts must be run in following order:

1. **A-beamprofile.m**
2. **B_quantumyield.m**
3. **C_spectra.m**
4. **D_emissionspectra.m**
5. **E_autocorrelation.m**
6. **F_makefigures.m**

## Folder structure

+ **datalake** contains all the data used by the scripts,
+ **plots**, figures presented in the paper,
+ **utils** contains a class with useful equations,
+ **loaders**, classes which load, transform and plot data.
+ **aux-files**, json files with initial parameters consumed by the main scripts.