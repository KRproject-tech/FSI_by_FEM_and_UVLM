%% parameters
param_ver = 9.0;                                        %% parameter file version


%% 解析条件
End_Time = 30;                                          %% Nondimensional analysis time [-]
d_t = 1.5e-3;                                           %% Nondimensional step time [-]
core_num = 6;                                           %% Core number [-]
speed_check = 0;                                        %% 1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:implicit solver，0:explicit solver [-]

Ma = 1.0;                                              %% Mass ratio [-]
Ua = 8;                                               %% Nondimensional flow velocity [-]
theta_a = 0*1e-1;                                      %% Nondimensional material damping [-]: σ = D(ε + θ^* dtε), θ^* := (θ [1/s])*(U_in/L [s])
C_theta_a = 0*1e-2;                                  	%% Cθ^* := Cθ/(ρf*L^3*W*Uin) [-]
J_a = 0;                                                %% J^* := J/(ρf*L^4*W) [-]

dt_rz_end = 0*3e-1;                                     %% 末端(X=1 [-])での無次元化初期速度 [-] (X方向に比例的に分布) ---------------- 初期速度を与えるとインパルス的な応答になる．
q_in_norm = @( time)( 0.1*sin( pi*time/0.2).*( time < 0.2 ) );              %% Initial disturbance
q_in_vec = [ 0 0 1].';                                                      %% Force direction [-]  

mode_num = 10;                                          %% 計算モード数

%% plot条件

i_snapshot = 50;                           %% スナップショットplot刻み
Snapshot_tmin = 4;                         %% スナップショットplot開始時間 [-]
Snapshot_tmax = End_Time;                   %% スナップショットplot終了時間 [-]
panel_node_plot = 0;                        %% 平板のパネルノード・コロケーション点・流体力ベクトル のplot: 1:ON, 0:OFF [-]
pressure_interp_plot = 0;                   %% 補間流体力ベクトル [Pa]
% movie_format = 'mpeg';                      %% 動画保存形式 [-]
% movie_format = 'avi';
movie_format = 'wmv';                      



%% 平板パラメータ

mu_m = 1/Ma;                                    %% 無次元化平板密度 [-]
nu = 0.3;                                       %% ポアソン比 [-]

Length = 1.0;                                   %% (Nondimensional) length [-]
Width = 2.00*Length;                           	%% aspect ratio [-]
thick = 1.0e-3*Length;                          %% thickness [-]
Aa = Width*thick;                               %% 断面積 [-]
Ia = Width*thick^3/12;                          %% 断面二次モーメント [-]
eta_m = mu_m/Ua^2;                              %% 無次元化曲げ剛性 [-]
zeta_m = Aa/Ia*eta_m*Length^2;                	%% 無次元化伸び剛性 [-]
Nx = 8;                                        %% x方向要素数 [-]
Ny = 14;                                         %% y方向要素数 [-]
n_LW = 1.0;                                     %% 不等間隔メッシュパラメータ [-]
N_gauss = 5;                                    %% Gauss-Legendre求積の次数 [-]



F_in = -0*[ 0 0 1].';                           %% 無次元化体積力 [-]



x_vec = ( (0:Nx)/Nx ).^n_LW*Length;             %% ノードのx座標(要素座標系)の刻み範囲 [-]
y_vec = (0:Ny)/Ny*Width;                        %% ノードのy座標(要素座標系)の刻み範囲 [-]
N_element = Nx*Ny;                             	%% 全要素数 [-]

Dp_mat = 1/(1 - nu^2)*[	1   nu  0;
                        nu  1   0;
                        0   0   (1 - nu)/2];





%% 流体パラメータ

U_in = 1.0;                                     %% X方向流入速度 [-]

V_in = ones(N_element,1)*U_in*[ 1 0 0];         %% 流入流速ベクトル [-]
%%[*] 渦点同氏の近接に由来する特異値を除去
%%% W. Hoydonck et al, Validity of Viscous Core Correction Models for
%%% Self-Induced Velocity Calculations, the Journal of the American Helicopter Society, pp. 1-8, 2011.
r_eps.fine = 1e-6;                             %% Wake格子の最大長を基準にした Vortex core size [-]: (平板上コロケーション点上の循環算出時)
r_eps.rough = 10e-2;                          	%% Wake格子の最大長を基準にした Vortex core size [-]: (Wakeの時間発展時)
Ncore = 2;                                      %% Dimension of vortex core model
%%[*] 分母の外積の特異点を除去
%%% R. Leuthold, Multiple-Wake Vortex Lattice Method for Membrane-Wing
%%% Kites, Master of Science Thesis, p. 90, 2015.
%%% (渦コア半径未満10^-10より大きい値にする．10^-10より大きい値でないと非対称の流れ場を引き起こしうる．)
eps_v = 1e-9;         


dL_vec = diff( x_vec);
dt_wake = dL_vec(end)/U_in;                            	%% Wake放出刻み時間 [-] (dt_wake = dL/Uin)
dt_wake_per_dt = ceil( dt_wake/d_t);                    %% 流体ソルバー刻み時間/構造ソルバー刻み時間 [-]

%% 平板境界条件

%%[0] 0値代入
node_r_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                       %% Node number giving the displacement constraint [-]
node_dxr_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                     %% Node number giving x-directional gradient constraint [-]
node_dyr_0 = [ 1:Ny+1:(Nx+1)*(Ny+1)];                                     %% Node number giving y-directional gradient constraint [-]

%%[1] 同一値に拘束 (回転自由の剛体棒に柔軟シートが固定されている場合)
node_dxr_theta_c = [ ];                             %% x方向勾配を同一に拘束するノード番号 [-]
node_dyr_theta_c = [ ];                             %% y方向勾配を同一に拘束するノード番号 [-]

%%[2] 前縁回転ダンパの支持区間
element_C_theta = 1:Ny;                             	%% 前縁回転ダンパに支持された要素番号 [-]

%% Wake解析パラメータ
R_wake_x_threshold = 5.5*Length;                                    %% Wake末端位置の閾値 [m]
R_wake_x_threshold_no_change = R_wake_x_threshold - 1.5*Length;     %% 変形を許容するWake末端位置の閾値 [m]


%% global 変数
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


