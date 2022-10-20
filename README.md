# FSI_by_FEM_and_UVLM
Fluid-Structure Interaction Analysis Using FEM and UVLM:

source code for Matlab (Windows): FSI analysis for the flapping sheet (This code is for MATLAB R2007b – )

## Overview
Fluid flow is modeled using the unsteady vortex lattice method (UVLM). A flexible sheet is modeled using the finite element method (FEM) with absolute nodal coordinate formulation (ANCF) for the shell element. This is done to reproduce the deformation of a plate while considering the spanwise deformation and geometrical nonlinearity.

## Preparation before analysis
__[Step 1] Install the ToolBoxes__

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

__[Step 1.2] Add path to installed ToolBoxes__

Modify "add_pathes.m".

__[Step 1.3] Modify the source code in the “TriStream” ToolBox__

"TriStream.m" must be modified to plot the stram line Line.

````
Line 49: x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=u(:)'; v=v(:)';
````
→
````
 x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=double(u(:)'); v=double(v(:)');
````
and 
````
line 63: TRI = tsearch(x,y,tri',Xbeg,Ybeg);
````
→ 
````
TRI = tsearchn([x.' y.'], tri',[Xbeg.' Ybeg.']);
````

__[Step 2] Start GUI form__

Open the “GUI.fig” from MATLAB.

![タイトルなし](https://user-images.githubusercontent.com/114337358/192756887-25b36670-8faa-423f-b535-63a536ced8c8.png)

__[Step 2.1] Pre-setting__

Push the "Parameters" buttun and edit parameters.

__[Step 3] Start analysis__

Push the “exe” button and wait until the finish of the analysis.

__[Step 4] Plot results__

Push the “plot” button.


## Parameters

Analytical condisions are in "./save/param_setting.m"

````
%% Analytical conditions
End_Time = 20;                                          %% Nondimensional analysis time [-]
d_t = 1.0e-3;                                           %% Nondimensional step time [-]
core_num = 6;                                           %% Core number [-]
speed_check = 0;                                        %% 1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:implicit solver，0:explicit solver [-]

Ma = 1.0;                                              %% Mass ratio [-]
Ua = 15.0;                                            	%% Nondimensional flow velocity [-]
theta_a_vec = 0e-1*[ 0 10];                           	%% Nondimensional material damping [-]
````
where

* __Mass ratio $M^*$__ is the density ratio of the fluid and sheet
* __Nondimensional flow velocity $U^*$__ represents the free-stream velocity nondimensionalized by the flag rigidity and inertia

 Initial disturbances on two sheets to break the trivial equilibrium are written as, 

````
q_in_norm = @( time)( 0.5*sin( pi*time/0.2).*( time < 0.2 ) );              %% Initial disturbance (upper sheet)
q_in_norm_1 = @( time)( 0.0*sin( pi*time/0.2).*( time < 0.2 ) );        	%% Initial disturbance (lower sheet)
q_in_vec = [ 0 0 1].';                                                      %% Force direction [-]  
````


## Image

![Velocity_field](https://user-images.githubusercontent.com/114337358/192750314-cb1e90ff-6000-4cc9-8b85-8bcad371dddc.png)
Streamlines around flapping sheets (not pathlines)

![displacement](https://user-images.githubusercontent.com/114337358/195290706-db1c0575-f07e-4d27-9696-a1ede965fedd.png)
Wake behind sheets

![snapshot](https://user-images.githubusercontent.com/114337358/195290303-25102659-399d-45a3-b99a-e861bdb5a68e.png)
Snapshot of two flapping sheets

## Demonstration movie

https://youtu.be/9FOBByPBSeA

## References

[1] Influence of the aspect ratio of the sheet for an electric generator utilizing the rotation of a flapping sheet, Mechanical Engineering Journal, Vol. 8, No. 1 (2021).

https://doi.org/10.1299/mej.20-00459

[2] Flow-induced vibration and energy-harvesting performance analysis for parallelized two flutter-mills considering span-wise plate deformation with geometrical nonlinearity and three-dimensional flow, International Journal of Structural Stability and Dynamics, Vol. 22, No. 14, (2022).

https://doi.org/10.1142/S0219455422501632

[3] Influence of boundary conditions on a flutter-mill, Journal of Sound and Vibration, Vol. 478, No. 21 (2020).

https://doi.org/10.1016/j.jsv.2020.115359



