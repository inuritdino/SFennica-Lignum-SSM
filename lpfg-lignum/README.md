LPFG-LIGNUM is a model to simulate pine tree (and others) development over time.
It takes the physiology and radiation into account to produce the 3D shape of a tree.
This model based on that reported in:
1) Perttunen et al, LIGNUM: a model combining the structure and the functioning of trees,
*Ecol. Modelling*, **108**, 189-198, 1998.
2) Sievanen et al, Towards extension of a single tree model to stand level,
*Func. Plant Biol.*, **35**, 964-975 2008.

One main difference from the model present in the aforementioned papers is that LPFG-LIGNUM
uses so called shadow propagation model for radiation conditions (see Palubicki et al, Self-organizing
tree models for image synthesis, *ACM Transactions on Graphics* **28**(3), 58:1-10, 2009).
For this, the voxel space is created and then the radition conditions for each segment are calculated
by means of the shadow, propagating down the voxel layers from the shading segment.

This folder contains the actual LPFG-LIGNUM code (C++ and L+C) to run directly in the LPFG simulator.
The LPFG simulator is available under L-studio (Win) or VLab (Mac/Unix) simulation environment.
See <http://algorithmicbotany.org/virtual_laboratory/> for the software documentation, download, and updates.
