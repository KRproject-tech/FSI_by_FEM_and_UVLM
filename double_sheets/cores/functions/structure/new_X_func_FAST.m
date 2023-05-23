function [ out, out1] = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, stage, out1)


%% [0] 変数抽出
node_r_0 = var_param.node_r_0;
node_dxr_0 = var_param.node_dxr_0;
node_dyr_0 = var_param.node_dyr_0;

node_r_0_1 = var_param.node_r_0_1;
node_dxr_0_1 = var_param.node_dxr_0_1;
node_dyr_0_1 = var_param.node_dyr_0_1;

node_r_marge = var_param.node_r_marge;
node_dxr_marge = var_param.node_dxr_marge;
node_dyr_marge = var_param.node_dyr_marge;


coordinates = var_param.coordinates;

d_t = var_param.d_t;
alpha_v = var_param.alpha_v;
theta_a_vec = var_param.theta_a_vec;

N_qi = var_param.N_qi; 
N_q_all = var_param.N_q_all;


flag_theta_a = prod( abs( theta_a_vec));                                    %% いずれかのシートの減衰が0の時，0値をとる．


%%[*] 1枚目
M_global = m_global_struct.M_global;                                        %% 質量行列     : M (質量＋自身のシートからの付加質量効果 [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_1*d_t^2 q_1)
M_global1 = m_global_struct.M_global1;                                      %% 質量行列     : M ( 2枚目のシートからの付加質量効果: [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_2*d_t^2 q_2 )
Qf_global = qf_global_struct.Qf_global;                                     %% 体積力       : Q_f
dq_Qe_global = dq_qe_global_struct.dq_Qe_global;                            %% dq_Qe(q(n))
Qe_global = qe_global_struct.Qe_global;                                     %% 剛性力       : Q_e
Qd_global = qd_global_struct.Qd_global;                                     %% 減衰力       : Q_d

%%[*] 2枚目
M_global_1 = m_global_struct.M_global_1;                                 	%% 質量行列     : M (質量＋自身のシートからの付加質量効果 [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_2*d_t^2 q_2)
M_global1_1 = m_global_struct.M_global1_1;                                 	%% 質量行列     : M (1枚目のシートからの付加質量効果 [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_1*d_t^2 q_1)
Qf_global_1 = qf_global_struct.Qf_global_1;                              	%% 体積力       : Q_f
dq_Qe_global_1 = dq_qe_global_struct.dq_Qe_global_1;                       	%% dq_Qe(q(n))
Qe_global_1 = qe_global_struct.Qe_global_1;                              	%% 剛性力       : Q_e
Qd_global_1 = qd_global_struct.Qd_global_1;                              	%% 減衰力       : Q_d


%% [1] 境界条件

%%[*] 0値固定 (シート1)
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





%%[*] 0値固定 (シート2)
%%[1-0] 変位境界条件
if ~isempty( node_r_0_1)
    i_r_1 = repmat( ( N_qi*(node_r_0_1 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0_1) 1]);                  %% 変位拘束をかけるノードに対応するx,y変位成分番号(z=0 [m])
    i_r_1 = reshape(i_r_1.',1,[]);
else
    i_r_1 = [];
end

%%[1-1] 勾配境界条件 (x方向)
if ~isempty( node_dxr_0_1)
    i_dx_r_1 = repmat( ( N_qi*(node_dxr_0_1 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0_1) 1]);           %% dx_r = [1 0 0]^T．
    i_dx_r_1 = reshape(i_dx_r_1.',1,[]);
else
    i_dx_r_1 = [];
end

%%[1-2] 勾配境界条件 (y方向)
if ~isempty( node_dyr_0_1)
    i_dy_r_1 = repmat( ( N_qi*(node_dyr_0_1 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0_1) 1]);           %% dy_r = [0 1 0]^T．
    i_dy_r_1 = reshape(i_dy_r_1.',1,[]);
else
    i_dy_r_1 = [];
end





%%[*] 節点値共有
%%[1-0] 変位境界条件
if ~isempty( node_r_marge)
    i_r_marge = repmat( ( N_qi*(node_r_marge - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_marge) 1]);                  
    i_r_marge = reshape(i_r_marge.',1,[]);
else
    i_r_marge = [];
end

%%[1-1] 勾配境界条件 (x方向)
if ~isempty( node_dxr_marge)
    i_dx_r_marge = repmat( ( N_qi*(node_dxr_marge - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_marge) 1]);           
    i_dx_r_marge = reshape(i_dx_r_marge.',1,[]);
else
    i_dx_r_marge = [];
end

%%[1-2] 勾配境界条件 (y方向)
if ~isempty( node_dyr_marge)
    i_dy_r_marge = repmat( ( N_qi*(node_dyr_marge - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_marge) 1]);           
    i_dy_r_marge = reshape(i_dy_r_marge.',1,[]);
else
    i_dy_r_marge = [];
end


%%[*] 節点値の固定
%%[*-0] シート1枚目
i_vec = [ i_r i_dx_r i_dy_r];
%%[*-1] シート2枚目
i_vec_1 = [ i_r_1 i_dx_r_1 i_dy_r_1];


%%[*] 節点共有    
i_vec_marge = [ i_r_marge i_dx_r_marge i_dy_r_marge];



i_vec_all = union( i_vec_1, i_vec_marge).';
if size( i_vec_all, 2) == 1
    i_vec_all = i_vec_all.';
end



%% グローバル行列の組立と境界条件の反映


 

%%[*] 質量行列
M_mat_BC = M_global;
M_mat_BC_1 = M_global_1;
%%[*] 節点共有の境界条件
%%[*] 1枚目
U_M_mat = M_global1;
U_M_mat(i_vec_marge,:) = M_mat_BC_1(i_vec_marge,:);
M_mat_BC(i_vec_marge,i_vec_marge) = M_mat_BC(i_vec_marge,i_vec_marge) + M_mat_BC_1(i_vec_marge,i_vec_marge);

%%[*] 2枚目
L_M_mat = M_global1_1;
L_M_mat(:,i_vec_marge) = M_mat_BC_1(:,i_vec_marge);

%%[*] 節点値固定の境界条件
M_GLOBAL = [    M_mat_BC  	U_M_mat;
                L_M_mat    	M_mat_BC_1];
            
M_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
M_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];





%%[*] 外力+曲げ剛性
Q_global = (Qf_global - Qe_global);                                                                     %% 内力＋外力項
Q_global_1 = (Qf_global_1 - Qe_global_1);         

%%[*] 節点共有の境界条件
Q_global(i_vec_marge) = 2*Q_global(i_vec_marge);

%%[*] 節点値固定の境界条件
Q_GLOBAL = [    Q_global;
                Q_global_1];
Q_GLOBAL([ i_vec N_q_all+i_vec_all]) = [];




%%[*] 減衰行列
Qd_mat = Qd_global;
Qd_mat_1 = Qd_global_1;
%%[*] 節点共有の境界条件
%%[*] 1枚目
U_M_mat = 0*Qd_mat;
U_M_mat(i_vec_marge,:) = Qd_mat_1(i_vec_marge,:);
Qd_mat(i_vec_marge,i_vec_marge) = Qd_mat(i_vec_marge,i_vec_marge) + Qd_mat_1(i_vec_marge,i_vec_marge);

%%[*] 2枚目
L_M_mat = 0*Qd_mat_1;
L_M_mat(:,i_vec_marge) = Qd_mat_1(:,i_vec_marge);

%%[*] 節点値固定の境界条件
Qd_GLOBAL = [ 	Qd_mat  	U_M_mat;
                L_M_mat    	Qd_mat_1];
            
Qd_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
Qd_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];





%%[*] 膜剛性行列の勾配行列
dq_Qe_mat = dq_Qe_global;
dq_Qe_mat_1 = dq_Qe_global_1;
%%[*] 節点共有の境界条件
%%[*] 1枚目
U_M_mat = 0*dq_Qe_mat;
U_M_mat(i_vec_marge,:) = dq_Qe_mat_1(i_vec_marge,:);
dq_Qe_mat(i_vec_marge,i_vec_marge) = dq_Qe_mat(i_vec_marge,i_vec_marge) + dq_Qe_mat_1(i_vec_marge,i_vec_marge);

%%[*] 2枚目
L_M_mat = 0*dq_Qe_mat_1;
L_M_mat(:,i_vec_marge) = dq_Qe_mat_1(:,i_vec_marge);

%%[*] 節点値固定の境界条件
dq_Qe_GLOBAL = [ 	dq_Qe_mat  	U_M_mat;
                    L_M_mat    	dq_Qe_mat_1];
            
dq_Qe_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
dq_Qe_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];




eye_mat = speye( 2*N_q_all);
eye_mat([ i_vec N_q_all+i_vec_all],:) = [];
eye_mat(:,[ i_vec N_q_all+i_vec_all]) = [];
zero_mat = sparse( 2*N_q_all-length( [ i_vec N_q_all+i_vec_all]),  2*N_q_all-length( [ i_vec N_q_all+i_vec_all]));




%% [2] 加速度dtt_q算出

if stage == 0                                                                                           %% 1回目の計算値を使いまわす．
    C_damp = (flag_theta_a == 0)*2 + ~(flag_theta_a == 0)*1;                                          	%% シート両方に減衰があるときは半陰解法
    D_matrix = [ eye_mat                                    zero_mat;
                 (Qd_GLOBAL + C_damp*d_t/2*dq_Qe_GLOBAL)  	M_GLOBAL];                                  %% 減衰0で単位行列になる．
     out1.D_matrix = D_matrix;
else
    D_matrix = out1.D_matrix;
end
X2_matrix = [ zero_mat      eye_mat;
              zero_mat      zero_mat];
A_mat1 = D_matrix - alpha_v*d_t*X2_matrix;
A_mat2 = D_matrix + (1 - alpha_v)*d_t*X2_matrix;




%%[2-0] 固定部のノードでは dt_q = 0, dtt_q = 0
not_i_vec = (1:N_q_all);
not_i_vec(i_vec) = [];

not_i_vec_all = (1:N_q_all);
not_i_vec_all(i_vec_all) = [];


%% [3] 状態ベクトル更新

i_vec_X = [ not_i_vec N_q_all+not_i_vec_all 2*N_q_all+not_i_vec 3*N_q_all+not_i_vec_all];

if stage == 0                                                                                           %% 1回目の計算値を使いまわす．
    out1.A1_A2_Xn = A_mat1\( A_mat2*X_vec(i_vec_X) );
end
out_0 = out1.A1_A2_Xn + A_mat1\( d_t*[ zero_mat(:,1);
                                       Q_GLOBAL]);
                    
out = X_vec;
out(i_vec_X) = out_0;
out(N_q_all+[ i_vec_marge 2*N_q_all+i_vec_marge]) = out([ i_vec_marge 2*N_q_all+i_vec_marge]);


end