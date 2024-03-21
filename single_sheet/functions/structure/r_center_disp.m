function [ X_center_disp, Z_center_disp] = r_center_disp( data)

N_element = data.N_element;
Nx = data.Nx;
Ny = data.Ny;
N_qi = data.N_qi;
N_q = data.N_q;
nodes = data.nodes;
h_X_vec = data.h_X_vec;


%% [1] ノード変位データ抽出 (span方向中央)

if mode( Ny, 2)
    %%[*] スパン方向が奇数要素数の時 (各ノードでの位置を平均化)
    
    %% r_X
    %%[*] 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )   
    idx_center_element = ceil( Ny/2):Ny:N_element-ceil( Ny/2);
    node_num = reshape( nodes(idx_center_element,:).', [], 1);

    i_vec = repmat( ( N_qi*(node_num - 1)+1 ), [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( node_num) 1]);
    i_vec = reshape( i_vec.',1,[]);
    %%%
    %%% ^  Y
    %%% | q_(4) ---- q_(3)
    %%% |  |           |
    %%% | q_(1) ---- q_(2)
    %%% O -------------------> X 
    %%%
    i_vec_1elem_X = i_vec(1:N_q:end);                       %% １要素における第1ノードのX座標
    i_vec_2elem_X = i_vec(N_qi+1:N_q:end);              	%% １要素における第2ノードのX座標
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% １要素における第3ノードのX座標
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% １要素における第4ノードのX座標
    
    X_center_disp1 = h_X_vec([ i_vec_1elem_X i_vec_2elem_X(end)],:);
    X_center_disp2 = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);

    X_center_disp =(X_center_disp1 + X_center_disp2)/2;
    
    %% r_Z
    %%[*] 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )   
    i_vec_1elem_Z = i_vec(3:N_q:end);                       %% １要素における第1ノードのZ座標
    i_vec_2elem_Z = i_vec(N_qi+3:N_q:end);              	%% １要素における第2ノードのZ座標
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% １要素における第3ノードのZ座標
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% １要素における第4ノードのZ座標
    
    Z_center_disp1 = h_X_vec([ i_vec_1elem_Z i_vec_2elem_Z(end)],:);
    Z_center_disp2 = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
    
    Z_center_disp = (Z_center_disp1 + Z_center_disp2)/2;
else
    %%[*] スパン方向が偶数要素数の時
    
    %% r_X
    %%[*] 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )   
    idx_center_element = Ny/2:Ny:N_element-Ny/2;
    node_num = reshape( nodes(idx_center_element,:).', [], 1);

    i_vec = repmat( ( N_qi*(node_num - 1)+1 ), [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( node_num) 1]);
    i_vec = reshape( i_vec.',1,[]);
    %%%
    %%% ^  Y
    %%% | q_(4) ---- q_(3)
    %%% |  |           |
    %%% | q_(1) ---- q_(2)
    %%% O -------------------> X 
    %%%
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% １要素における第3ノードのX座標
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% １要素における第4ノードのX座標
    
    X_center_disp = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);


    %% r_Z
    %%[*] 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )   
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% １要素における第3ノードのZ座標
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% １要素における第4ノードのZ座標
    
    Z_center_disp = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
end


