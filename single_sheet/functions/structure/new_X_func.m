function out = new_X_func( X_vec, M_global, Qf_global, dq_Qe_global, Qe_global, Qd_global, var_param)


%% [0] �ϐ����o
node_r_0 = var_param.node_r_0;
node_dxr_0 = var_param.node_dxr_0;
node_dyr_0 = var_param.node_dyr_0;
coordinates = var_param.coordinates;
d_t = var_param.d_t;
alpha_v = var_param.alpha_v;

N_qi = var_param.N_qi; 
N_q_all = var_param.N_q_all;






%% [1] ���E����
%%[1-0] �ψʋ��E����
if ~isempty( node_r_0)
    i_r = repmat( ( N_qi*(node_r_0 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0) 1]);        %% �ψʍS����������m�[�h�ɑΉ�����x,y�ψʐ����ԍ�(z=0 [m])
    i_r = reshape(i_r.',1,[]);
else
    i_r = [];
end

%%[1-1] ���z���E���� (x����)
if ~isempty( node_dxr_0)
    i_dx_r = repmat( ( N_qi*(node_dxr_0 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0) 1]);	%% dx_r = [1 0 0]^T�D
    i_dx_r = reshape(i_dx_r.',1,[]);
else
    i_dx_r = [];
end

%%[1-2] ���z���E���� (y����)
if ~isempty( node_dyr_0)
    i_dy_r = repmat( ( N_qi*(node_dyr_0 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0) 1]);	%% dy_r = [0 1 0]^T�D
    i_dy_r = reshape(i_dy_r.',1,[]);
else
    i_dy_r = [];
end

i_vec = [ i_r i_dx_r i_dy_r];

%% [2] �����xdtt_q�Z�o

Q_global = (Qf_global - Qe_global);                                                                     %% ���́{�O�͍�


M_global(i_vec,:) = [];
M_global(:,i_vec) = [];
Q_global(i_vec) = []; 

Qd_global(i_vec,:) = [];
Qd_global(:,i_vec) = [];
dq_Qe_global(i_vec,:) = [];
dq_Qe_global(:,i_vec) = [];

eye_mat = speye( N_q_all-length( i_vec));
zero_mat = sparse( N_q_all-length( i_vec),  N_q_all-length( i_vec));
D_matrix = [ eye_mat                                    zero_mat;
             M_global\(Qd_global + d_t*dq_Qe_global)  	eye_mat];                                           %% ����0�ŒP�ʍs��ɂȂ�D
X2_matrix = [ zero_mat      eye_mat;
              zero_mat      zero_mat];
A_mat1 = D_matrix - alpha_v*d_t*X2_matrix;
A_mat2 = D_matrix + (1 - alpha_v)*d_t*X2_matrix;




%%[2-0] �Œ蕔�̃m�[�h�ł� dt_q = 0, dtt_q = 0
n_i_vec = (1:N_q_all);
n_i_vec(i_vec) = [];


%% [3] ��ԃx�N�g���X�V

out_0 = A_mat1\(A_mat2*X_vec([ n_i_vec N_q_all+n_i_vec]) ...
                           ...
                            + d_t*[ zero_mat(:,1);
                                    M_global\Q_global]);
                    
out =  X_vec;
out([ n_i_vec N_q_all+n_i_vec]) = out_0;


end