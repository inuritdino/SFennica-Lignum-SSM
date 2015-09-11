# SFennica-Lignum-SSM
Repository for the supplementary code of the paper 'Data-based stochastic modeling of tree growth and structure formation' in Silva Fennica, 2015.

The paper offers a novel approach to model structure of trees by iteratively optimizing empirical distributions describing structure of the data tree (as obtained from the laser scanning) and those of a Stochastic Structural tree Model (SSM, here stochastic LIGNUM). The resulting SSM has 'optimal' parameter values to produce trees statistically similar to the data tree and to each other, that is the SSM sample trees are not exact copies of each other. Thus, unlimited number of sample trees can be generated.

Structure of the repository:

LPFG-LIGNUM: the L+C code of the stochastic LIGNUM (one needs LPFG simulator to run it, see below)

LPFG-LIGNUM-MATLAB: the Matlab interface to run the LPFG-Lignum code

OPTIM-TOOL-MATLAB: optimization routine code including the distance functions

Requirements:

1. Matlab, The MathWorks. The paper simluations were run in Matlab ver. R2014b.
2. LPFG, simulator for plant growth modeling. The simulator is available as a part of L-studio (Windows) or VLab (Unix/Mac) simulations environments. Both can be downloaded at <http://algorithmicbotany.org/virtual_laboratory/> along with documentation and detailed installation guides.
