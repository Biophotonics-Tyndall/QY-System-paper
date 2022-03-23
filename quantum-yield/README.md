# QY evaluation of UCNPS using a QY system built at Tyndall National institute

Author: Jean Matias \
Email: jean.matias@tyndall.ie \
Publication: Optics Express

## Installing instructions

The code is fully built in Python v3.9 and it is recommended to start a virtual environment to
avoid compatibility issues. The following steps will guide you for that.

1. First install Python 3.9 and the package manager pip
2. Install the virtual environment package in a cmd prompt (Windows) or a terminal (Unix): `pip install virtualenv`
3. Start a virtual environment: `virtualenv -p <path\to\>python3.9 qyvenv`
4. Activate it: `.\qyvenev\Scripts\activate`
5. Install requirements: `pip -r install .\requirements.txt`
6. Add the qyvenv to the jupyter kernels list: `python -m ipykernel install --user --name=qyvenv`
7. If you need to remove it: `jupyter-kernelspec uninstall qyvenv`
8. Start jupyter notebook: `jupyter-notebook`
9. Navigate to the notebooks and select the active kernel: qyvenv
10. Once you finish, deactivate the qyvenv: `deactivate`

## Folder structure

+ **data** contains all the data used to calculate the QY of the UCNPs, including data to calibrate the system,
+ **docs**, relevant documentation and a pdf of the main jupyter notebook with the QY analysis,
+ **notebooks**, jupyter notebooks with with analysis and calibration of the QY system,
+ **plots**, figures presented in the paper and supplementary material,
+ **scripts**, Python scripts with all the functions, equations, and calculations needed for the analysis. 

