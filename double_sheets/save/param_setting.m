%% parameters
param_ver = 1.0;                                        %% parameter file version

%% 解析条件
End_Time = 20;                                          %% Nondimensional analysis time [-]
d_t = 1.0e-3;                                           %% Nondimensional step time [-]
core_num = 6;                                           %% Core number [-]
speed_check = 0;                                        %% 1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:implicit solver，0:explicit solver [-]

Ma = 1.0;                                               %% Mass ratio [-]
Ua = 15.0;                                            	%% Nondimensional flow velocity [-]
theta_a_vec = 0e-1*[ 0 10];                           	%% Nondimensional material damping [-]



q_in_norm = @( time)( 0.5*sin( pi*time/0.2).*( time < 0.2 ) );              %% Initial disturbance (upper sheet)
q_in_norm_1 = @( time)( 0.0*sin( pi*time/0.2).*( time < 0.2 ) );        	%% Initial disturbance (lower sheet)
q_in_vec = [ 0 0 1].';                                                      %% Force direction [-]  

mode_num = 5;

%% plot conditions

i_snapshot = 50;                            %% Snapshot plot ticking
Snapshot_tmin = 34;                         %% Snapshot plot start time [-].
Snapshot_tmax = End_Time;                   %% End time [-]
panel_node_plot = 0;                        %% 1:ON, 0:OFF [-]
pressure_interp_plot = 0;                   %%

% movie_format = 'avi';
movie_format = 'wmv';                      





%% 平板パラメータ

mu_m = 1/Ma;                                    %% 1/"mass ratio" [-]
nu = 0.3;                                       %% Poisson's ratio [-]

Length = 1.0;                                   %% (Nondimensional) length [-]
Width = 1.0;                                    %% aspect ratio [-]
Height = 2.0;                              	%% distance between two sheets [-]
thick = 1e-3;                                   %% thickness [-]
Aa = Width*thick;                               %%  [-]
Ia = Width*thick^3/12;                          %%  [-]
eta_m = mu_m/Ua^2;                              %%  [-]
zeta_m = Aa/Ia*eta_m*Length^2;                	%%  [-]
Nx = 14;                                        %% Number of elements in x direction [-]
Ny = 10;                                     	%% Number of elements in y direction [-]
n_LW = 1.0;                                     %%  [-]
N_gauss = 5;                                    %% Order of Gauss-Legendre quadrature [-]
N_gauss_RI = 5;                               	%% Order of Gauss-Legendre quadrature (for membrane stiffness calculation) [-]



F_in = 0*[ 0 1 0].';                            %% nondimensional volume force [-]



x_vec = ( (0:Nx)/Nx ).^n_LW*Length;             %%  [-]
y_vec = (0:Ny)/Ny*Width;                        %%  [-]
N_element = Nx*Ny;                             	%% Number of all elements [-]

Dp_mat = 1/(1 - nu^2)*[	1   nu  0;
                        nu  1   0;
                        0   0   (1 - nu)/2];





%% 流体パラメータ

U_in = 1.0;                                     %% Nondimensional velocity [-]

V_in = ones(N_element,1)*U_in*[ 1 0 0];         %% Velocity vector [-]
V_in_all = repmat( V_in, [ 2 1]);
%%% W. Hoydonck et al, Validity of Viscous Core Correction Models for
%%% Self-Induced Velocity Calculations, the Journal of the American Helicopter Society, pp. 1-8, 2011.
r_eps.fine = 1e-6;                              %% Vortex core size [-]: (Collocation points)
r_eps.rough = 10e-2;                          	%% Vortex core size [-]: (For Wake)
Ncore = 2;                                      %% Dimension of vortex core model
%%% R. Leuthold, Multiple-Wake Vortex Lattice Method for Membrane-Wing
%%% Kites, Master of Science Thesis, p. 90, 2015.
eps_v = 1e-9;         


dL_vec = diff( x_vec);
dt_wake = dL_vec(end)/U_in;                                         %% Wake step time [-] (dt_wake = dL/Uin)
dt_wake_per_dt = ceil( dt_wake/d_t);                                %%  [-]

%% Boundary conditions for two sheets

%%[0] Upper sheet
node_r_0 = [ 1:Ny+1 ];                                                      %% Node number giving the displacement constraint [-]
node_dxr_0 = [ 1:Ny+1 ];                                                    %% Node number giving x-directional gradient constraint [-]
node_dyr_0 = [ 1:Ny+1 ];                                                    %% Node number giving y-directional gradient constraint [-]

%%[1] Lower sheet
node_r_0_1 = [ 1:Ny+1 ];                                                    %% Node number giving the displacement constraint [-]
node_dxr_0_1 = [ 1:Ny+1 ];                                                  %% Node number giving x-directional gradient constraint [-]
node_dyr_0_1 = [ 1:Ny+1 ];                                                  %% Node number giving y-directional gradient constraint [-]

%%[2] 
node_r_marge = [ ];                                                       	%% Node number giving the displacement constraint [-]
node_dxr_marge = [ ];                                                       %% Node number giving x-directional gradient constraint [-]
node_dyr_marge = [ ];                                                       %% Node number giving y-directional gradient constraint [-]


%% Wake parameters
R_wake_x_threshold = 5.5*Length;                                            %% Wake end position threshold [-]
R_wake_x_threshold_no_change = R_wake_x_threshold - 1.5*Length;             %% Threshold for Wake end position that allows deformation [-]


%% global 
global var_param

var_param.Length = Length;
var_param.Nx = Nx;
var_param.d_t = d_t;
var_param.alpha_v = alpha_v;
var_param.theta_a_vec = theta_a_vec;


var_param.node_r_0 = node_r_0;
var_param.node_dxr_0 = node_dxr_0;
var_param.node_dyr_0 = node_dyr_0;

var_param.node_r_0_1 = node_r_0_1;
var_param.node_dxr_0_1 = node_dxr_0_1;
var_param.node_dyr_0_1 = node_dyr_0_1;

var_param.node_r_marge = node_r_marge;
var_param.node_dxr_marge = node_dxr_marge;
var_param.node_dyr_marge = node_dyr_marge;


var_param.r_eps = r_eps;
var_param.Ncore = Ncore;
var_param.eps_v = eps_v;


