# SFennica-Lignum-SSM
Repository for the supplementary code of the paper 'Data-based stochastic modeling of tree growth and structure formation' in _Silva Fennica_, 2015 (anonymous).

The paper offers a novel approach to model structure of trees by iteratively optimizing empirical distributions describing structure of the data tree (as obtained from the laser scanning) and those of a Stochastic Structural tree Model (SSM, here stochastic LIGNUM). The resulting SSM has 'optimal' parameter values to produce trees statistically similar to the data tree and to each other, that is the SSM sample trees are not exact copies of each other. Thus, unlimited number of sample trees can be generated.

####Structure of the repository

**LPFG-LIGNUM**: the L+C code of the stochastic LIGNUM (one needs LPFG simulator to run it, see below)

**LPFG-LIGNUM-MATLAB**: the Matlab interface to run the LPFG-Lignum code

**OPTIM-TOOL-MATLAB**: optimization routine code including the distance functions

**QSM-LIBRARY**: library of .mat (Matlab data) files containing structural data for the real trees

**LIB**: technical tools for geometry, scatter manipulations, tree class, structural distance functions etc.

####Requirements

1. Matlab, The MathWorks. The paper simluations were run in Matlab ver. R2014b.
2. LPFG, simulator for plant growth modeling. The simulator is available as a part of L-studio (Windows) or VLab (Unix/Mac) simulation environments. Both can be downloaded at <http://algorithmicbotany.org/virtual_laboratory/> along with the documentation and detailed installation guides.

####How to run

1. Navigate to the downloaded directory.
2. _test_sim.m_, _test_sim16.m_, and _real_sim.m_ files are the script files for the **test case** (Fig. 3 and Table 2), the **test case with 16 parameters** to estimate, and the **real case** (Fig. 4), respectively. You may want to run one of these just by typing their names at the Matlab prompt.
3. The resulting figures will appear (what they are can be found in the script files' comments).

####Common pitfalls

I. 'NOT FOUND LIBRARY/SYMBOL' error

MATLAB environment variables are different from those of the system's shell (command line). This produces a 'not found library
or symbol' error in MATLAB, when it calls LPFG simulator (LPFG is foreign to MATLAB). To fix:

1.Run the LPFG simulator within LPFG-LIGNUM folder by typing at the system's command line:
```
lpfg lignum.l view.v material.mat
```
If this runs successfully, the LPFG works. If not, perhaps system's PATH variable is not correctly set or other installation problems occurred.

2.Adjust MATLAB environment variables to those of the system's environment. Use 'setenv' and 'getenv' commands in MATLAB. You may want to look at the PATH, LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, DYLD_FRAMEWORK_PATH variables (these are common in Mac OS systems, use corresponding variables on other systems).

For example, in Mac OS the following was proven to work when typing at the MATLAB command line or in the startup.m script:
```
setenv('PATH',[getenv('PATH') ':/Applications/browser.app/Contents/MacOS/dbin' ':/Applications/browser.app/Contents/MacOS/bin']);% sets the path to the LPFG executables
setenv('DYLD_FRAMEWORK_PATH','<put here the system variable value>');
```
i.e. if in the system environment 
```
echo $DYLD_FRAMEWORK_PATH 
```
gives '/Users/X/lib', then the last command would be:
```
setenv('DYLD_FRAMEWORK_PATH','/Users/X/lib');
```

