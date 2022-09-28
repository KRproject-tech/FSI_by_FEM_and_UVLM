function [ out, out1] = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, stage, out1)


%% [0] �ϐ����o
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


flag_theta_a = prod( abs( theta_a_vec));                                    %% �����ꂩ�̃V�[�g�̌�����0�̎��C0�l���Ƃ�D


%%[*] 1����
M_global = m_global_struct.M_global;                                        %% ���ʍs��     : M (���ʁ{���g�̃V�[�g����̕t�����ʌ��� [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_1*d_t^2 q_1)
M_global1 = m_global_struct.M_global1;                                      %% ���ʍs��     : M ( 2���ڂ̃V�[�g����̕t�����ʌ���: [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_2*d_t^2 q_2 )
Qf_global = qf_global_struct.Qf_global;                                     %% �̐ϗ�       : Q_f
dq_Qe_global = dq_qe_global_struct.dq_Qe_global;                            %% dq_Qe(q(n))
Qe_global = qe_global_struct.Qe_global;                                     %% ������       : Q_e
Qd_global = qd_global_struct.Qd_global;                                     %% ������       : Q_d

%%[*] 2����
M_global_1 = m_global_struct.M_global_1;                                 	%% ���ʍs��     : M (���ʁ{���g�̃V�[�g����̕t�����ʌ��� [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_2*d_t^2 q_2)
M_global1_1 = m_global_struct.M_global1_1;                                 	%% ���ʍs��     : M (1���ڂ̃V�[�g����̕t�����ʌ��� [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_1*d_t^2 q_1)
Qf_global_1 = qf_global_struct.Qf_global_1;                              	%% �̐ϗ�       : Q_f
dq_Qe_global_1 = dq_qe_global_struct.dq_Qe_global_1;                       	%% dq_Qe(q(n))
Qe_global_1 = qe_global_struct.Qe_global_1;                              	%% ������       : Q_e
Qd_global_1 = qd_global_struct.Qd_global_1;                              	%% ������       : Q_d


%% [1] ���E����

%%[*] 0�l�Œ� (�V�[�g1)
%%[1-0] �ψʋ��E����
if ~isempty( node_r_0)
    i_r = repmat( ( N_qi*(node_r_0 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0) 1]);                  %% �ψʍS����������m�[�h�ɑΉ�����x,y�ψʐ����ԍ�(z=0 [m])
    i_r = reshape(i_r.',1,[]);
else
    i_r = [];
end

%%[1-1] ���z���E���� (x����)
if ~isempty( node_dxr_0)
    i_dx_r = repmat( ( N_qi*(node_dxr_0 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0) 1]);           %% dx_r = [1 0 0]^T�D
    i_dx_r = reshape(i_dx_r.',1,[]);
else
    i_dx_r = [];
end

%%[1-2] ���z���E���� (y����)
if ~isempty( node_dyr_0)
    i_dy_r = repmat( ( N_qi*(node_dyr_0 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0) 1]);           %% dy_r = [0 1 0]^T�D
    i_dy_r = reshape(i_dy_r.',1,[]);
else
    i_dy_r = [];
end





%%[*] 0�l�Œ� (�V�[�g2)
%%[1-0] �ψʋ��E����
if ~isempty( node_r_0_1)
    i_r_1 = repmat( ( N_qi*(node_r_0_1 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0_1) 1]);                  %% �ψʍS����������m�[�h�ɑΉ�����x,y�ψʐ����ԍ�(z=0 [m])
    i_r_1 = reshape(i_r_1.',1,[]);
else
    i_r_1 = [];
end

%%[1-1] ���z���E���� (x����)
if ~isempty( node_dxr_0_1)
    i_dx_r_1 = repmat( ( N_qi*(node_dxr_0_1 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0_1) 1]);           %% dx_r = [1 0 0]^T�D
    i_dx_r_1 = reshape(i_dx_r_1.',1,[]);
else
    i_dx_r_1 = [];
end

%%[1-2] ���z���E���� (y����)
if ~isempty( node_dyr_0_1)
    i_dy_r_1 = repmat( ( N_qi*(node_dyr_0_1 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0_1) 1]);           %% dy_r = [0 1 0]^T�D
    i_dy_r_1 = reshape(i_dy_r_1.',1,[]);
else
    i_dy_r_1 = [];
end





%%[*] �ߓ_�l���L
%%[1-0] �ψʋ��E����
if ~isempty( node_r_marge)
    i_r_marge = repmat( ( N_qi*(node_r_marge - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_marge) 1]);                  
    i_r_marge = reshape(i_r_marge.',1,[]);
else
    i_r_marge = [];
end

%%[1-1] ���z���E���� (x����)
if ~isempty( node_dxr_marge)
    i_dx_r_marge = repmat( ( N_qi*(node_dxr_marge - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_marge) 1]);           
    i_dx_r_marge = reshape(i_dx_r_marge.',1,[]);
else
    i_dx_r_marge = [];
end

%%[1-2] ���z���E���� (y����)
if ~isempty( node_dyr_marge)
    i_dy_r_marge = repmat( ( N_qi*(node_dyr_marge - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_marge) 1]);           
    i_dy_r_marge = reshape(i_dy_r_marge.',1,[]);
else
    i_dy_r_marge = [];
end


%%[*] �ߓ_�l�̌Œ�
%%[*-0] �V�[�g1����
i_vec = [ i_r i_dx_r i_dy_r];
%%[*-1] �V�[�g2����
i_vec_1 = [ i_r_1 i_dx_r_1 i_dy_r_1];


%%[*] �ߓ_���L    
i_vec_marge = [ i_r_marge i_dx_r_marge i_dy_r_marge];



i_vec_all = union( i_vec_1, i_vec_marge).';
if size( i_vec_all, 2) == 1
    i_vec_all = i_vec_all.';
end



%% �O���[�o���s��̑g���Ƌ��E�����̔��f


 

%%[*] ���ʍs��
M_mat_BC = M_global;
M_mat_BC_1 = M_global_1;
%%[*] �ߓ_���L�̋��E����
%%[*] 1����
U_M_mat = M_global1;
U_M_mat(i_vec_marge,:) = M_mat_BC_1(i_vec_marge,:);
M_mat_BC(i_vec_marge,i_vec_marge) = M_mat_BC(i_vec_marge,i_vec_marge) + M_mat_BC_1(i_vec_marge,i_vec_marge);

%%[*] 2����
L_M_mat = M_global1_1;
L_M_mat(:,i_vec_marge) = M_mat_BC_1(:,i_vec_marge);

%%[*] �ߓ_�l�Œ�̋��E����
M_GLOBAL = [    M_mat_BC  	U_M_mat;
                L_M_mat    	M_mat_BC_1];
            
M_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
M_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];





%%[*] �O��+�Ȃ�����
Q_global = (Qf_global - Qe_global);                                                                     %% ���́{�O�͍�
Q_global_1 = (Qf_global_1 - Qe_global_1);         

%%[*] �ߓ_���L�̋��E����
Q_global(i_vec_marge) = 2*Q_global(i_vec_marge);

%%[*] �ߓ_�l�Œ�̋��E����
Q_GLOBAL = [    Q_global;
                Q_global_1];
Q_GLOBAL([ i_vec N_q_all+i_vec_all]) = [];




%%[*] �����s��
Qd_mat = Qd_global;
Qd_mat_1 = Qd_global_1;
%%[*] �ߓ_���L�̋��E����
%%[*] 1����
U_M_mat = 0*Qd_mat;
U_M_mat(i_vec_marge,:) = Qd_mat_1(i_vec_marge,:);
Qd_mat(i_vec_marge,i_vec_marge) = Qd_mat(i_vec_marge,i_vec_marge) + Qd_mat_1(i_vec_marge,i_vec_marge);

%%[*] 2����
L_M_mat = 0*Qd_mat_1;
L_M_mat(:,i_vec_marge) = Qd_mat_1(:,i_vec_marge);

%%[*] �ߓ_�l�Œ�̋��E����
Qd_GLOBAL = [ 	Qd_mat  	U_M_mat;
                L_M_mat    	Qd_mat_1];
            
Qd_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
Qd_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];





%%[*] �������s��̌��z�s��
dq_Qe_mat = dq_Qe_global;
dq_Qe_mat_1 = dq_Qe_global_1;
%%[*] �ߓ_���L�̋��E����
%%[*] 1����
U_M_mat = 0*dq_Qe_mat;
U_M_mat(i_vec_marge,:) = dq_Qe_mat_1(i_vec_marge,:);
dq_Qe_mat(i_vec_marge,i_vec_marge) = dq_Qe_mat(i_vec_marge,i_vec_marge) + dq_Qe_mat_1(i_vec_marge,i_vec_marge);

%%[*] 2����
L_M_mat = 0*dq_Qe_mat_1;
L_M_mat(:,i_vec_marge) = dq_Qe_mat_1(:,i_vec_marge);

%%[*] �ߓ_�l�Œ�̋��E����
dq_Qe_GLOBAL = [ 	dq_Qe_mat  	U_M_mat;
                    L_M_mat    	dq_Qe_mat_1];
            
dq_Qe_GLOBAL([ i_vec N_q_all+i_vec_all],:) = [];
dq_Qe_GLOBAL(:,[ i_vec N_q_all+i_vec_all]) = [];




eye_mat = speye( 2*N_q_all);
eye_mat([ i_vec N_q_all+i_vec_all],:) = [];
eye_mat(:,[ i_vec N_q_all+i_vec_all]) = [];
zero_mat = sparse( 2*N_q_all-length( [ i_vec N_q_all+i_vec_all]),  2*N_q_all-length( [ i_vec N_q_all+i_vec_all]));




%% [2] �����xdtt_q�Z�o

if stage == 0                                                                                           %% 1��ڂ̌v�Z�l���g���܂킷�D
    C_damp = (flag_theta_a == 0)*2 + ~(flag_theta_a == 0)*1;                                          	%% �V�[�g�����Ɍ���������Ƃ��͔��A��@
    D_matrix = [ eye_mat                                    zero_mat;
                 (Qd_GLOBAL + C_damp*d_t/2*dq_Qe_GLOBAL)  	M_GLOBAL];                                  %% ����0�ŒP�ʍs��ɂȂ�D
     out1.D_matrix = D_matrix;
else
    D_matrix = out1.D_matrix;
end
X2_matrix = [ zero_mat      eye_mat;
              zero_mat      zero_mat];
A_mat1 = D_matrix - alpha_v*d_t*X2_matrix;
A_mat2 = D_matrix + (1 - alpha_v)*d_t*X2_matrix;




%%[2-0] �Œ蕔�̃m�[�h�ł� dt_q = 0, dtt_q = 0
not_i_vec = (1:N_q_all);
not_i_vec(i_vec) = [];

not_i_vec_all = (1:N_q_all);
not_i_vec_all(i_vec_all) = [];


%% [3] ��ԃx�N�g���X�V

i_vec_X = [ not_i_vec N_q_all+not_i_vec_all 2*N_q_all+not_i_vec 3*N_q_all+not_i_vec_all];

if stage == 0                                                                                           %% 1��ڂ̌v�Z�l���g���܂킷�D
    out1.A1_A2_Xn = A_mat1\( A_mat2*X_vec(i_vec_X) );
end
out_0 = out1.A1_A2_Xn + A_mat1\( d_t*[ zero_mat(:,1);
                                       Q_GLOBAL]);
                    
out = X_vec;
out(i_vec_X) = out_0;
out(N_q_all+[ i_vec_marge 2*N_q_all+i_vec_marge]) = out([ i_vec_marge 2*N_q_all+i_vec_marge]);


end