# Single molecule tools (Alex Payne-Dwyer)

Forked from [single-molecule-tools](https://awollman.github.io/single-molecule-tools/) (authors Adam Wollman and Isabel Llorente-Garcia, c. 2018).

This repository comprises a series of Matlab functions and scripts called **ADEMScode** for tracking and quantification of spots in images, deconvolution, segmentation and simulation. It is also capable of estimating single molecule concentrations from volumetric models and 2D images. It is broadly useful for localisation microscopy, step-wise photobleaching, STORM, PALM, FRAP and regular fluorescence microscopy data in 16-bit .tif, .czi and other BioFormat image stack types.  

## SlimVar: repository for ADEMScode v2.2

This section introduces the version of ADEMScode software used to support SlimVar microscopy analysis in the following publications:

1. [SlimVar: rapid in vivo single-molecule tracking of chromatin regulators in plants](https://www.biorxiv.org/content/10.1101/2024.05.17.594710.abstract)  
AL Payne-Dwyer, GJ Jang, C Dean and MC Leake  
bioRxiv, 2024.05.17.594710

2. [Multifunctional polymerization domains determine the onset of epigenetic silencing in Arabidopsis](https://www.biorxiv.org/content/10.1101/2024.02.15.580496.abstract)  
A Schulten, GJ Jang, A Payne-Dwyer, M Fiedler, ML Nielsen, M Bienz, MC Leake and C Dean  
bioRxiv, 2024.02. 15.580496

3. [Modular in vivo assembly of Arabidopsis FCA oligomers into condensates competent for RNA 3’ processing](https://doi.org/10.1038/s44318-025-00394-4)  
GJ Jang, AL Payne-Dwyer, R Maple, Z Wu, F Liu, SG Lopez, Y Wang, X Fang, MC Leake and C Dean  
**Now published in EMBO J** (2025) 44: 2056 - 2074

4. [Single-molecular quantification of flowering control proteins within nuclear condensates in live whole Arabidopsis root](https://link.springer.com/protocol/10.1007/978-1-0716-2221-6_21)
AL Payne-Dwyer and MC Leake
Methods in Molecular Biology (Volume 2476): Chromosome Architecture pp 311-328 (2022)
Also available via [open access preprint](https://arxiv.org/pdf/2108.13743)


Novel functionality here includes distribution-based estimates of periodic clustering in molecular assemblies, as well as scripts for visualising correlated properties of trajectories colocalised across two channels.

### Guide to setting up github source control in matlab


Here is a useful [MathWorks protocol](https://blogs.mathworks.com/community/2014/10/20/matlab-and-git/) to contribute to or update this repository directly from Matlab.
