![図1](https://user-images.githubusercontent.com/114337358/219939559-f0153a0e-2202-4760-9e00-2342ade1eb8a.png)
# <p align=center>FSI_by_FEM_and_UVLM</p>

**Communication**

<a style="text-decoration: none" href="https://twitter.com/hogelungfish_" target="_blank">
    <img src="https://img.shields.io/badge/twitter-%40hogelungfish_-1da1f2.svg" alt="Twitter">
</a>
<p>

**Language**
<p>
<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/matlab/matlab-original.svg" width="60"/>
<p>

__Strong-coupled Fluid-Structure Interaction Analysis (FSI) Using FEM and UVLM to analyze the flapping sheet under the post-flutter [^1][^2][^3]:  
source code for Matlab (Windows): FSI analysis for the flapping sheet under the post-flutter (This code is validated with MATLAB R2007b or later versions)__

## Overview
This is the numerical simulation code for a limit cycle oscillation of a rectangular sheet flapping in a three-dimensional flow. 
Fluid flow is modeled using the unsteady vortex lattice method (UVLM) [^5]. A flexible sheet is modeled using the finite element method (FEM) with absolute nodal coordinate formulation (ANCF) for the shell element [^4]. This is done to reproduce the deformation of a plate while considering the spanwise deformation and geometrical nonlinearity. Robust FSI analysis under a large mass ratio is achieved by the strong coupling between the fluid solver and the structure solver by introducing the explicit added mass calculation. 

Strong coupled FSI can achieve more robust numerical analysis under large fluid density than a loose coupling scheme.

[![](https://img.youtube.com/vi/heaMrV6I3RQ/0.jpg)](https://www.youtube.com/watch?v=heaMrV6I3RQ)

## Publications

If you use this work in an academic context, please cite the following publication(s):

* Influence of the aspect ratio of the sheet for an electric generator utilizing the rotation of a flapping sheet, Mechanical Engineering Journal, Vol. 8, No. 1 (2021).  
https://doi.org/10.1299/mej.20-00459

````
@article{Akio YAMANO202120-00459,
    title={Influence of the aspect ratio of the sheet for an electric generator utilizing the rotation of a flapping sheet},
    author={Akio YAMANO and Hiroshi IJIMA and Atsuhiko SHINTANI and Chihiro NAKAGAWA and Tomohiro ITO},
    journal={Mechanical Engineering Journal},
    volume={8},
    number={1},
    pages={20-00459-20-00459},
    year={2021},
    doi={10.1299/mej.20-00459}
}
````

* Flow-induced vibration and energy-harvesting performance analysis for parallelized two flutter-mills considering span-wise plate deformation with geometrical nonlinearity and three-dimensional flow, International Journal of Structural Stability and Dynamics, Vol. 22, No. 14, (2022).  
https://doi.org/10.1142/S0219455422501632

````
@article{doi:10.1142/S0219455422501632,
    author = {Yamano, Akio and Chiba, Masakatsu},
    title = {Flow-Induced Vibration and Energy-Harvesting Performance Analysis for Parallelized Two Flutter-Mills Considering Span-Wise Plate Deformation with Geometrical Nonlinearity and Three-Dimensional Flow},
    journal = {International Journal of Structural Stability and Dynamics},
    volume = {22},
    number = {14},
    pages = {2250163},
    year = {2022},
    doi = {10.1142/S0219455422501632}
}
````

* Influence of boundary conditions on a flutter-mill, Journal of Sound and Vibration, Vol. 478, No. 21 (2020).  
https://doi.org/10.1016/j.jsv.2020.115359

````
@article{YAMANO2020115359,
    title = {Influence of boundary conditions on a flutter-mill},
    journal = {Journal of Sound and Vibration},
    volume = {478},
    pages = {115359},
    year = {2020},
    doi = {https://doi.org/10.1016/j.jsv.2020.115359},
    author = {A. Yamano and A. Shintani and T. Ito and C. Nakagawa and H. Ijima}
}
````

## Directory    
<pre>
├─double_sheets
│  ├─cores
│  │  ├─functions
│  │  │  ├─fluid
│  │  │  └─structure
│  │  └─solver
│  │      ├─fluid
│  │      └─structure
│  ├─save
│  │  └─fig
│  │      └─modes
│  └─ToolBoxes
└─single_sheet
    ├─functions
    │  ├─fluid
    │  └─structure
    ├─save
    │  └─fig
    │      └─modes
    ├─solver
    │  ├─fluid
    │  └─structure
    └─ToolBoxes
</pre>
    
## Preparation before analysis
__[Step 1] Install the ToolBoxes__

The following ToolBoxes in “./XXXX/ToolBoxes/” ("XXXX" is "double_sheets" and "single_sheet") are required,

__For numerical analysis:__
*	__“Meshing a plate using four noded elements”__ by KSSV:  
https://jp.mathworks.com/matlabcentral/fileexchange/33731-meshing-a-plate-using-four-noded-elements

*	__“Sparse sub access”__ by Bruno Luong:  
https://jp.mathworks.com/matlabcentral/fileexchange/23488-sparse-sub-access

*	__“Vectorized Multi-Dimensional Matrix Multiplication”__ by Darin Koblick:  
https://jp.mathworks.com/matlabcentral/fileexchange/47092-vectorized-multi-dimensional-matrix-multiplication?s_tid=prof_contriblnk

__For plotting results:__
*	__“TriStream”__ by Matthew Wolinsky:  
https://jp.mathworks.com/matlabcentral/fileexchange/11278-tristream

*	__“mmwrite”__ by Micah Richert:  
https://jp.mathworks.com/matlabcentral/fileexchange/15881-mmwrite

*	__“Quiver 5”__ by Bertrand Dano:  
https://jp.mathworks.com/matlabcentral/fileexchange/22351-quiver-5?s_tid=FX_rc3_behav


__[Step 1.2] Add path to installed ToolBoxes__

Modify "add_pathes.m" to add path to abovementined installed ToolBoxes as follows,
````
addpath ./ToolBoxes/XX;
````
where `XX` is the name of folder of the installed ToolBox.

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
    
__[Step 5] View plotted results__

Results (figures and movie) plotted by [Step 4] are in "./XXXX/save" directory.
    

## Parameters

Analytical condisions are in "./XXXX/save/param_setting.m"

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


* __Mass ratio $M^*$__ is the density ratio of the fluid and sheet, which is defined by,

$$M^* := \dfrac{1}{\mu} = \dfrac{\rho_f L}{\rho_m h},$$

where $\rho_f$ and $\rho_m$ are the density of fluid and sheet, respectively. $L$ is the total length of sheet, and $h$ is the thickness of sheet.


* __Nondimensional flow velocity $U^*$__ represents the free-stream velocity nondimensionalized by the flag rigidity and inertia [^1][^2][^3],
    
$$U^* := \sqrt{ \dfrac{\mu}{\eta}} = \sqrt{ \dfrac{ \rho_m h L^2 H U^2_\infty }{ EI }}, $$

where $H$ is the width, $E$ is the Young's modulus, and $I := Hh^3/12$ is second moment of area of the sheet, respectively. 
    

 Initial disturbances on two sheets to break the trivial equilibrium are written as, 
````
q_in_norm = @( time)( 0.5*sin( pi*time/0.2).*( time < 0.2 ) );              %% Initial disturbance (upper sheet)
q_in_norm_1 = @( time)( 0.0*sin( pi*time/0.2).*( time < 0.2 ) );        	%% Initial disturbance (lower sheet)
q_in_vec = [ 0 0 1].';                                                      %% Force direction [-]  
````


Dimensions of sheets are defined by,

````
Length = 1.0;                                   %% (Nondimensional) length [-]
Width = 1.0;                                    %% aspect ratio [-]
Height = 2.0;                              	    %% distance between two sheets [-]
thick = 1e-3;                                   %% thickness [-]

````
where the aspect ratio `Width` is $H^* := H/L$, the nondimensional thickness `thick` is $h^* := h/L$, and the nondimensional distance between two sheets `Height` is $D^* := D/L$.
 

Boundary conditions for two sheets are written as,

* __Clamped at the leading-edge__
````
%% Boundary conditions for two sheets

%%[0] Upper sheet
node_r_0 = [ 1:Ny+1 ];                                                      %% Node number giving the displacement constraint [-]
node_dxr_0 = [ 1:Ny+1 ];                                                    %% Node number giving x-directional gradient constraint [-]
node_dyr_0 = [ 1:Ny+1 ];                                                    %% Node number giving y-directional gradient constraint [-]

%%[1] Lower sheet
node_r_0_1 = [ 1:Ny+1 ];                                                    %% Node number giving the displacement constraint [-]
node_dxr_0_1 = [ 1:Ny+1 ];                                                  %% Node number giving x-directional gradient constraint [-]
node_dyr_0_1 = [ 1:Ny+1 ];                                                  %% Node number giving y-directional gradient constraint [-]

````

* __Pinned at the leading-edge__
````
%% Boundary conditions for two sheets

%%[0] Upper sheet
node_r_0 = [ 1:Ny+1 ];                                                      %% Node number giving the displacement constraint [-]
node_dxr_0 = [ ];                                                           %% Node number giving x-directional gradient constraint [-]
node_dyr_0 = [ 1:Ny+1 ];                                                    %% Node number giving y-directional gradient constraint [-]

%%[1] Lower sheet
node_r_0_1 = [ 1:Ny+1 ];                                                    %% Node number giving the displacement constraint [-]
node_dxr_0_1 = [ ];                                                         %% Node number giving x-directional gradient constraint [-]
node_dyr_0_1 = [ 1:Ny+1 ];                                                  %% Node number giving y-directional gradient constraint [-]

````
where index in vector shows the node index around a plate element to apply boundary conditions.

![タイトルなし](https://user-images.githubusercontent.com/114337358/196866330-b2dec9e7-cacc-441a-9c69-da409bc73a81.png)
(b) Shell elements on a plate. A flexible sheet is partitioned by shell elements when the number of elements is $N_{elem} := N_x N_y$ [^2].


## Gallery

![Velocity_field](https://user-images.githubusercontent.com/114337358/192750314-cb1e90ff-6000-4cc9-8b85-8bcad371dddc.png)
__Streamlines around flapping sheets (not streaklines)__

![displacement](https://user-images.githubusercontent.com/114337358/195290706-db1c0575-f07e-4d27-9696-a1ede965fedd.png)
__Wake behind sheets__

![snapshot](https://user-images.githubusercontent.com/114337358/195290303-25102659-399d-45a3-b99a-e861bdb5a68e.png)
__Snapshot of two flapping sheets__

![comparison](https://user-images.githubusercontent.com/114337358/233918767-5ddeca5e-bd84-4662-b5f2-995c8b75cd7e.png)
__Comparisons of snapshot of a flapping sheet under various $U^*$ between numerical results (above) [^3] and experimental results (below) [^6]__


## Demonstration movie

[![](https://img.youtube.com/vi/9FOBByPBSeA/0.jpg)](https://youtu.be/9FOBByPBSeA)

### References  
[^1]: Influence of the aspect ratio of the sheet for an electric generator utilizing the rotation of a flapping sheet, Mechanical Engineering Journal, Vol. 8, No. 1 (2021).  
https://doi.org/10.1299/mej.20-00459

[^2]: Flow-induced vibration and energy-harvesting performance analysis for parallelized two flutter-mills considering span-wise plate deformation with geometrical nonlinearity and three-dimensional flow, International Journal of Structural Stability and Dynamics, Vol. 22, No. 14, (2022).  
https://doi.org/10.1142/S0219455422501632

[^3]: Influence of boundary conditions on a flutter-mill, Journal of Sound and Vibration, Vol. 478, No. 21 (2020).  
https://doi.org/10.1016/j.jsv.2020.115359

[^4]: A. Shabana, Computational Continuum Mechanics, Chap. 6 (Cambridge University Press, 2008), pp. 231–285.

[^5]: J. Katz and A. Plotkin, Low-Speed Aerodynamics, (Cambridge University Press, New York, 2001).

[^6]: M. Chen, L. Jia, Y. Wu, X. Yin, Y. Ma, Bifurcation and chaos of a flag in an inviscid flow, J. Fluid Struct. 45 (2014b) 124-137.



