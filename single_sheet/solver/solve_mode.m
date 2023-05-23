%% solve mode

%% 平衡点　q0 = q(0), dt_q0 = 0
q_vec = h_X_vec(1:N_q_all,1);
dt_q_vec = 0*q_vec;


%% 粘弾性ベクトルの組立 Qe^(n)
flag_output = 1;                                                            %% 剛性行列算出の有効化
theta_a = 1;                                                                %% ∂q_ε^TD∂q_ε, ∂q_κ^TD∂q_κ算出の有効化
generate_stiff_matrices;
theta_a = var_param.theta_a;

%% 線形化剛性行列

%%[0-0] 膜剛性
dq_e_Dp_dq_e_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);  

    int_dq_e_Dp_dq_e = Qd_eps_mat_i(:,:,ii);    
    dq_e_Dp_dq_e_global(i_vec,i_vec) = dq_e_Dp_dq_e_global(i_vec,i_vec) + squeeze( zeta_m*int_dq_e_Dp_dq_e);
    
end   
    
K_e_m_mat =  dq_e_Dp_dq_e_global;

%%[0-1] 曲げ剛性
dq_k_Dp_dq_k_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);  

    int_dq_k_Dp_dq_k = Qd_k_mat_i(:,:,ii);    
    dq_k_Dp_dq_k_global(i_vec,i_vec) = dq_k_Dp_dq_k_global(i_vec,i_vec) + squeeze( eta_m*int_dq_k_Dp_dq_k);
    
end   
    
K_e_k_mat =  dq_k_Dp_dq_k_global;


%% 境界条件 (前縁角度の同一値拘束は考慮しない)

%%[*] 0値固定
%%[1-0] 変位境界条件
if ~isempty( node_r_0)
    i_r = repmat( ( N_qi*(node_r_0 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0) 1]);                  %% 変位拘束をかけるノードに対応するx,y変位成分番号(z=0 [m])
    i_r = reshape(i_r.',1,[]);
else
    i_r = [];
end

%%[1-1] 勾配境界条件 (x方向)
if ~isempty( node_dxr_0)
    i_dx_r = repmat( ( N_qi*(node_dxr_0 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0) 1]);           %% dx_r = [1 0 0]^T．
    i_dx_r = reshape(i_dx_r.',1,[]);
else
    i_dx_r = [];
end

%%[1-2] 勾配境界条件 (y方向)
if ~isempty( node_dyr_0)
    i_dy_r = repmat( ( N_qi*(node_dyr_0 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0) 1]);           %% dy_r = [0 1 0]^T．
    i_dy_r = reshape(i_dy_r.',1,[]);
else
    i_dy_r = [];
end


i_vec = [ i_r i_dx_r i_dy_r];

not_i_vec = (1:N_q_all);
not_i_vec(i_vec) = [];

M_global_BC = M_global;
M_global_BC(i_vec,:) = [];
M_global_BC(:,i_vec) = [];

K_e_m_mat_BC = K_e_m_mat;
K_e_m_mat_BC(i_vec,:) = [];
K_e_m_mat_BC(:,i_vec) = [];

K_e_k_mat_BC = K_e_k_mat;
K_e_k_mat_BC(i_vec,:) = [];
K_e_k_mat_BC(:,i_vec) = [];


%% modal analysis

%%[2-0] 小さい順にmode_num個の固有値を算出
[ Phi_dq_mat, omega_a2] = eigs( (mu_m*M_global_BC)\(K_e_m_mat_BC + K_e_k_mat_BC), mode_num, 'SM');   

%%[2-1] 無次元化固有振動数：ω^* := L/U_in*ω [-]
omega_a = sqrt( diag( omega_a2));

%%[2-1] 無次元化変位固有モード: Δq(t)に対応
Phi_dq_mat_BC = zeros(N_q_all,mode_num);
Phi_dq_mat_BC(not_i_vec,:) = Phi_dq_mat;

%%[2-2] 無次元化位置固有モード: q(t) = q0 + Δq(t)に対応
Phi_q_mat_BC = repmat( q_vec, [ 1 mode_num]) + Phi_dq_mat_BC;



%% Modal damping ratio

%%[3-0] モード減衰比 [-]
zeta_n = theta_a*omega_a/2;                      



disp( 'Nondimensional natural frequency: ')
disp( omega_a)
disp( 'Modal damping ratio: ')
disp( zeta_n)


