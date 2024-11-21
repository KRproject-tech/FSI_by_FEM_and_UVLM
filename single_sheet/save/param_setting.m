%% parameters
param_ver = 12.0;                                        %% parameter file version

%% ��͏���
End_Time = 30;                                          %% ����������͎��� [-]
% End_Time = 10;                                          %% ����������͎��� [-]
d_t = 1.5e-3;                                           %% �����������ݎ��� [-]
core_num = 8;                                           %% CPU�� [-]
speed_check = 0;                                        %% �v�Z���x�m�F�F1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:�A��@�C0:�z��@ [-]
coupling_flag = 1;                                      %% 1: ���A���i�t�����ʗz�I�v�Z�j�C0: ��A���i�݂��Ⴂ�v�Z�j


Ma = 1.0;                                              %% ���ʔ� [-]
Ua = 25;                                                %% ������������ [-]
% theta_a = 1e-2;                                       %% �\�������p�����[�^ [-]: �� = D(�� + ��^* dt��), ��^* := (�� [1/s])*(U_in/L [s])
theta_v = 7.8e-4/25;                                   %% �L�����\�������W�� [s]�~�W�� [1/s] (��^*�~��_n(7)^*/2 = 0.001, ��_n(7)^* = 2.574 [-] @ U^* = 25 �� ��^* = 0.00078 @ U^* = 25)
theta_a = 0*(Ua/Ma^2)*theta_v;                           %% �\�������p�����[�^ [-]: �� = D(�� + ��^* dt��), ��^* := (�� [1/s])*(U_in/L [s]) = (U^*/M*^2)*��(��_f^4 E/12��_m^5 H^2)�� �� (U^*/M*^2)*��
C_theta_a = 0*1e-2;                                     	%% �O����]�_���p�W�� [-]: C��^* := C��/(��f*L^3*W*Uin) [-]
J_a = 0;                                                %% ���������������[�����g [-]: J^* := J/(��f*L^4*W) [-]

dt_rz_end = 0*3e-1;                                     %% ���[(X=1 [-])�ł̖��������������x [-] (X�����ɔ��I�ɕ��z) ---------------- �������x��^����ƃC���p���X�I�ȉ����ɂȂ�D
q_in_norm = @( time)( 0.5*sin( pi*time/0.2).*( time < 0.2 ) );              %% �O���̐ϗ� [-] ------------------------------------------- �Ȃ߂炩�ȉ����ɂȂ�D
q_in_vec = [ 0 0 1].';                                                      %% �O���̐ϗ͂̍�p���� [-]  

mode_num = 5;

%% plot����

i_snapshot = 50;                           %% �X�i�b�v�V���b�gplot����
Snapshot_tmin = 25;                         %% �X�i�b�v�V���b�gplot�J�n���� [-]
Snapshot_tmax = End_Time;                   %% �X�i�b�v�V���b�gplot�I������ [-]
panel_node_plot = 0;                        %% ���̃p�l���m�[�h�E�R���P�[�V�����_�E���̗̓x�N�g�� ��plot: 1:ON, 0:OFF [-]
pressure_interp_plot = 0;                   %% ��ԗ��̗̓x�N�g�� [Pa]
% movie_format = 'mpeg';                      %% ����ۑ��`�� [-]
% movie_format = 'avi';
movie_format = 'wmv';                      



%% ���p�����[�^

mu_m = 1/Ma;                                    %% �������������x [-]
nu = 0.3;                                       %% �|�A�\���� [-]

Length = 1.0;                                   %% ���������S�� [-]
Width = 1.00*Length;                           	%% ����������(�A�X�y�N�g��) [-]
% Width = 50.0*Length;                           	%% ����������(�A�X�y�N�g��) [-]

thick = 1e-3;%(1.21/1.39e+3)*1/Ma;                    %% ������������ [-] (�����ȂقǍ\����͂̍��ݎ��Ԃ�傫���ł���D)
Aa = Width*thick;                               %% �f�ʐ� [-]
Ia = Width*thick^3/12;                          %% �f�ʓ񎟃��[�����g [-]
eta_m = mu_m/Ua^2;                              %% ���������Ȃ����� [-]
zeta_m = Aa/Ia*eta_m*Length^2;                	%% ���������L�э��� [-]
Nx = 15;                                        %% x�����v�f�� [-]
Ny = 10;                                         %% y�����v�f�� [-]
% Nx = 14;                                        %% x�����v�f�� [-]
% Ny = 40;                                         %% y�����v�f�� [-]

n_LW = 1.0;                                     %% �s���Ԋu���b�V���p�����[�^ [-]
N_gauss = 5;                                    %% Gauss-Legendre���ς̎��� [-]



k_gravity = 0.136^4/(1.21^3*4.989e-4);
% F_in = -9.807*k_gravity*(Ma/Ua)^2*[ 0 1 0].';                           %% ���������̐ϗ� [-]
F_in = 0*[ 0 1 0].';                           %% ���������̐ϗ� [-]


x_vec = ( (0:Nx)/Nx ).^n_LW*Length;             %% �m�[�h��x���W(�v�f���W�n)�̍��ݔ͈� [-]
y_vec = (0:Ny)/Ny*Width;                        %% �m�[�h��y���W(�v�f���W�n)�̍��ݔ͈� [-]
N_element = Nx*Ny;                             	%% �S�v�f�� [-]

Dp_mat = 1/(1 - nu^2)*[	1   nu  0;
                        nu  1   0;
                        0   0   (1 - nu)/2];



%% ���̗͌v�Z�̂�
flag_fluid_bench = 0;                                                       %% 1: ���̗͌v�Z�̂݁i�\���͍��̕��ŉ�́j�C0:FSI���

k_omega = 1/2;
Theta_pitch = pi/6.0;                                                        %% [rad]

omega_pitch = 2*k_omega;

theta_pitch_time = @( time)( Theta_pitch*sin( omega_pitch*time) );                          %% Pitch�p�x [rad]
dt_theta_pitch_time = @( time)( Theta_pitch*omega_pitch.*cos( omega_pitch*time));          	%% Pitch�p���x [rad/-]
dtt_theta_pitch_time = @( time)( -Theta_pitch*omega_pitch.^2.*sin( omega_pitch*time));     	%% Pitch�p�����x [rad/-]

L_pitch_center =  Length/2;                                                 %% Pitch��]���S�ʒu [-]


R_pitch = @( theta)[ ...
	cos( theta)     0   sin( theta) ;
    0               1   0           ;
    -sin( theta)    0   cos( theta) ];                  %% 3D��]�s��

dt_R_pitch = @( theta, omega) ...
    omega*[ ...
	-sin( theta)   	0   cos( theta) ;
    0               0   0           ;
    -cos( theta)    0   -sin( theta) ];                  %% 3D��]�s��


H_1_2nd = @( k)besselh( 1, 2, k); 
H_0_2nd = @( k)besselh( 0, 2, k); 

C_theodorsen = @( k)( H_1_2nd( k)/( H_1_2nd( k) + 1i*H_0_2nd( k) ) ); 


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
node_r_0 = [ 1:Ny+1];                                     %% �ψʍS����^����m�[�h�ԍ� [-]
node_dxr_0 = [ 1:Ny+1];                                  	%% x�������z�S����^����m�[�h�ԍ� [-]
node_dyr_0 = [ 1:Ny+1];                                   %% y�������z�S����^����m�[�h�ԍ� [-]

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


