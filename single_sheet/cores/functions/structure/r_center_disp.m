function [ X_center_disp, Z_center_disp] = r_center_disp( data)

N_element = data.N_element;
Nx = data.Nx;
Ny = data.Ny;
N_qi = data.N_qi;
N_q = data.N_q;
nodes = data.nodes;
h_X_vec = data.h_X_vec;


%% [1] �m�[�h�ψʃf�[�^���o (span��������)

if mode( Ny, 2)
    %%[*] �X�p����������v�f���̎� (�e�m�[�h�ł̈ʒu�𕽋ω�)
    
    %% r_X
    %%[*] 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %%[*] 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )   
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
    i_vec_1elem_X = i_vec(1:N_q:end);                       %% �P�v�f�ɂ������1�m�[�h��X���W
    i_vec_2elem_X = i_vec(N_qi+1:N_q:end);              	%% �P�v�f�ɂ������2�m�[�h��X���W
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% �P�v�f�ɂ������3�m�[�h��X���W
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% �P�v�f�ɂ������4�m�[�h��X���W
    
    X_center_disp1 = h_X_vec([ i_vec_1elem_X i_vec_2elem_X(end)],:);
    X_center_disp2 = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);

    X_center_disp =(X_center_disp1 + X_center_disp2)/2;
    
    %% r_Z
    %%[*] 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %%[*] 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )   
    i_vec_1elem_Z = i_vec(3:N_q:end);                       %% �P�v�f�ɂ������1�m�[�h��Z���W
    i_vec_2elem_Z = i_vec(N_qi+3:N_q:end);              	%% �P�v�f�ɂ������2�m�[�h��Z���W
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% �P�v�f�ɂ������3�m�[�h��Z���W
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% �P�v�f�ɂ������4�m�[�h��Z���W
    
    Z_center_disp1 = h_X_vec([ i_vec_1elem_Z i_vec_2elem_Z(end)],:);
    Z_center_disp2 = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
    
    Z_center_disp = (Z_center_disp1 + Z_center_disp2)/2;
else
    %%[*] �X�p�������������v�f���̎�
    
    %% r_X
    %%[*] 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %%[*] 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )   
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
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% �P�v�f�ɂ������3�m�[�h��X���W
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% �P�v�f�ɂ������4�m�[�h��X���W
    
    X_center_disp = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);


    %% r_Z
    %%[*] 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %%[*] 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )   
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% �P�v�f�ɂ������3�m�[�h��Z���W
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% �P�v�f�ɂ������4�m�[�h��Z���W
    
    Z_center_disp = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
end


