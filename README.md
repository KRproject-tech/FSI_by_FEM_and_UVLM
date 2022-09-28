# FSI_by_FEM_and_UVLM
Fluid-Structure Interaction Analysis Using FEM and UVLM:

source code for Matlab (Windows): FSI analysis for the flapping sheet (This code is for MATLAB R2007b – R2018a)

## Overview
Fluid flow is modeled using the unsteady vortex lattice method (UVLM). A flexible sheet is modeled using the finite element method (FEM) with absolute nodal coordinate formulation (ANCF) for the shell element. This is done to reproduce the deformation of a plate while considering the spanwise deformation and geometrical nonlinearity.

## Preparation before analysis
__[Step 1]__ Install the ToolBoxes

The following ToolBoxes in “./ ToolBoxes/” are required,
*	“Meshing a plate using four noded elements” by KSSV
https://jp.mathworks.com/matlabcentral/fileexchange/33731-meshing-a-plate-using-four-noded-elements
*	“mmwrite” by Micah Richert:
https://jp.mathworks.com/matlabcentral/fileexchange/15881-mmwrite
*	“Quiver 5” by Bertrand Dano:
https://jp.mathworks.com/matlabcentral/fileexchange/22351-quiver-5?s_tid=FX_rc3_behav
*	“Sparse sub access” by Bruno Luong: 
https://jp.mathworks.com/matlabcentral/fileexchange/23488-sparse-sub-access
*	“TriStream” by Matthew Wolinsky:
https://jp.mathworks.com/matlabcentral/fileexchange/11278-tristream
*	“Vectorized Multi-Dimensional Matrix Multiplication” by Darin Koblick:
https://jp.mathworks.com/matlabcentral/fileexchange/47092-vectorized-multi-dimensional-matrix-multiplication?s_tid=prof_contriblnk

__[Step 1.2]__ Modify the source code



__[Step 2]__ Start GUI form

Open the “GUI.fig” from MATLAB.

__[Step 3]__ Start analysis

Push the “exe” button and wait until the finish of the analysis.

__[Step 4]__ Plot results

Push the “plot” button.


## Image
![Velocity_field](https://user-images.githubusercontent.com/114337358/192750314-cb1e90ff-6000-4cc9-8b85-8bcad371dddc.png)
Streamlines around flapping sheets (not pathlines)
