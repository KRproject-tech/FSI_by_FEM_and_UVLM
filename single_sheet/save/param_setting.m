%% parameters
param_ver = 9.0;                                        %% parameter file version


%% ��͏���
End_Time = 30;                                          %% Nondimensional analysis time [-]
d_t = 1.5e-3;                                           %% Nondimensional step time [-]
core_num = 6;                                           %% Core number [-]
speed_check = 0;                                        %% 1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:implicit solver�C0:explicit solver [-]

Ma = 1.0;                                              %% Mass ratio [-]
Ua = 8;                                               %% Nondimensional flow velocity [-]
theta_a = 0*1e-1;                                      %% Nondimensional material damping [-]: �� = D(�� + ��^* dt��), ��^* := (�� [1/s])*(U_in/L [s])
C_theta_a = 0*1e-2;                                  	%% C��^* := C��/(��f*L^3*W*Uin) [-]
J_a = 0;                                                %% J^* := J/(��f*L^4*W) [-]

dt_rz_end = 0*3e-1;                                     %% ���[(X=1 [-])�ł̖��������������x [-] (X�����ɔ��I�ɕ��z) ---------------- �������x��^����ƃC���p���X�I�ȉ����ɂȂ�D
q_in_norm = @( time)( 0.1*sin( pi*time/0.2).*( time < 0.2 ) );              %% Initial disturbance
q_in_vec = [ 0 0 1].';                                                      %% Force direction [-]  

mode_num = 10;                                          %% �v�Z���[�h��

%% plot����

i_snapshot = 50;                           %% �X�i�b�v�V���b�gplot����
Snapshot_tmin = 4;                         %% �X�i�b�v�V���b�gplot�J�n���� [-]
Snapshot_tmax = End_Time;                   %% �X�i�b�v�V���b�gplot�I������ [-]
panel_node_plot = 0;                        %% ���̃p�l���m�[�h�E�R���P�[�V�����_�E���̗̓x�N�g�� ��plot: 1:ON, 0:OFF [-]
pressure_interp_plot = 0;                   %% ��ԗ��̗̓x�N�g�� [Pa]
% movie_format = 'mpeg';                      %% ����ۑ��`�� [-]
% movie_format = 'avi';
movie_format = 'wmv';                      



%% ���p�����[�^

mu_m = 1/Ma;                                    %% �������������x [-]
nu = 0.3;                                       %% �|�A�\���� [-]

Length = 1.0;                                   %% (Nondimensional) length [-]
Width = 2.00*Length;                           	%% aspect ratio [-]
thick = 1.0e-3*Length;                          %% thickness [-]
Aa = Width*thick;                               %% �f�ʐ� [-]
Ia = Width*thick^3/12;                          %% �f�ʓ񎟃��[�����g [-]
eta_m = mu_m/Ua^2;                              %% ���������Ȃ����� [-]
zeta_m = Aa/Ia*eta_m*Length^2;                	%% ���������L�э��� [-]
Nx = 8;                                        %% x�����v�f�� [-]
Ny = 14;                                         %% y�����v�f�� [-]
n_LW = 1.0;                                     %% �s���Ԋu���b�V���p�����[�^ [-]
N_gauss = 5;                                    %% Gauss-Legendre���ς̎��� [-]



F_in = -0*[ 0 0 1].';                           %% ���������̐ϗ� [-]



x_vec = ( (0:Nx)/Nx ).^n_LW*Length;             %% �m�[�h��x���W(�v�f���W�n)�̍��ݔ͈� [-]
y_vec = (0:Ny)/Ny*Width;                        %% �m�[�h��y���W(�v�f���W�n)�̍��ݔ͈� [-]
N_element = Nx*Ny;                             	%% �S�v�f�� [-]

Dp_mat = 1/(1 - nu^2)*[	1   nu  0;
                        nu  1   0;
                        0   0   (1 - nu)/2];





%% ���̃p�����[�^

U_in = 1.0;                                     %% X�����������x [-]

V_in = ones(N_element,1)*U_in*[ 1 0 0];         %% ���������x�N�g�� [-]
%%[*] �Q�_�����̋ߐڂɗR��������ْl������
%%% W. Hoydonck et al, Validity of Viscous Core Correction Models for
%%% Self-Induced Velocity Calculations, the Journal of the American Helicopter Society, pp. 1-8, 2011.
r_eps.fine = 1e-6;                             %% Wake�i�q�̍ő咷����ɂ��� Vortex core size [-]: (����R���P�[�V�����_��̏z�Z�o��)
r_eps.rough = 10e-2;                          	%% Wake�i�q�̍ő咷����ɂ��� Vortex core size [-]: (Wake�̎��Ԕ��W��)
Ncore = 2;                                      %% Dimension of vortex core model
%%[*] ����̊O�ς̓��ٓ_������
%%% R. Leuthold, Multiple-Wake Vortex Lattice Method for Membrane-Wing
%%% Kites, Master of Science Thesis, p. 90, 2015.
%%% (�Q�R�A���a����10^-10���傫���l�ɂ���D10^-10���傫���l�łȂ��Ɣ�Ώ̗̂����������N��������D)
eps_v = 1e-9;         


dL_vec = diff( x_vec);
dt_wake = dL_vec(end)/U_in;                            	%% Wake���o���ݎ��� [-] (dt_wake = dL/Uin)
dt_wake_per_dt = ceil( dt_wake/d_t);                    %% ���̃\���o�[���ݎ���/�\���\���o�[���ݎ��� [-]

%% �����E����

%%[0] 0�l���
node_r_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                       %% Node number giving the displacement constraint [-]
node_dxr_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                     %% Node number giving x-directional gradient constraint [-]
node_dyr_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                     %% Node number giving y-directional gradient constraint [-]

%%[1] ����l�ɍS�� (��]���R�̍��̖_�ɏ_��V�[�g���Œ肳��Ă���ꍇ)
node_dxr_theta_c = [ ];                             %% x�������z�𓯈�ɍS������m�[�h�ԍ� [-]
node_dyr_theta_c = [ ];                             %% y�������z�𓯈�ɍS������m�[�h�ԍ� [-]

%%[2] �O����]�_���p�̎x�����
element_C_theta = 1:Ny;                             	%% �O����]�_���p�Ɏx�����ꂽ�v�f�ԍ� [-]

%% Wake��̓p�����[�^
R_wake_x_threshold = 5.5*Length;                                    %% Wake���[�ʒu��臒l [m]
R_wake_x_threshold_no_change = R_wake_x_threshold - 1.5*Length;     %% �ό`�����e����Wake���[�ʒu��臒l [m]


%% global �ϐ�
global var_param

var_param.Length = Length;
var_param.Nx = Nx;
var_param.d_t = d_t;
var_param.alpha_v = alpha_v;
var_param.theta_a = theta_a;
var_param.C_theta_a = C_theta_a;
var_param.J_a = J_a;

var_param.node_r_0 = node_r_0;
var_param.node_dxr_0 = node_dxr_0;
var_param.node_dyr_0 = node_dyr_0;
var_param.node_dxr_theta_c = node_dxr_theta_c;
var_param.node_dyr_theta_c = node_dyr_theta_c;

var_param.r_eps = r_eps;
var_param.Ncore = Ncore;
var_param.eps_v = eps_v;


