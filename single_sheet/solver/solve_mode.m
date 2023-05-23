%% solve mode

%% ���t�_�@q0 = q(0), dt_q0 = 0
q_vec = h_X_vec(1:N_q_all,1);
dt_q_vec = 0*q_vec;


%% �S�e���x�N�g���̑g�� Qe^(n)
flag_output = 1;                                                            %% �����s��Z�o�̗L����
theta_a = 1;                                                                %% ��q_��^TD��q_��, ��q_��^TD��q_�ȎZ�o�̗L����
generate_stiff_matrices;
theta_a = var_param.theta_a;

%% ���`�������s��

%%[0-0] ������
dq_e_Dp_dq_e_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);  

    int_dq_e_Dp_dq_e = Qd_eps_mat_i(:,:,ii);    
    dq_e_Dp_dq_e_global(i_vec,i_vec) = dq_e_Dp_dq_e_global(i_vec,i_vec) + squeeze( zeta_m*int_dq_e_Dp_dq_e);
    
end   
    
K_e_m_mat =  dq_e_Dp_dq_e_global;

%%[0-1] �Ȃ�����
dq_k_Dp_dq_k_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);  

    int_dq_k_Dp_dq_k = Qd_k_mat_i(:,:,ii);    
    dq_k_Dp_dq_k_global(i_vec,i_vec) = dq_k_Dp_dq_k_global(i_vec,i_vec) + squeeze( eta_m*int_dq_k_Dp_dq_k);
    
end   
    
K_e_k_mat =  dq_k_Dp_dq_k_global;


%% ���E���� (�O���p�x�̓���l�S���͍l�����Ȃ�)

%%[*] 0�l�Œ�
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

%%[2-0] ����������mode_num�̌ŗL�l���Z�o
[ Phi_dq_mat, omega_a2] = eigs( (mu_m*M_global_BC)\(K_e_m_mat_BC + K_e_k_mat_BC), mode_num, 'SM');   

%%[2-1] ���������ŗL�U�����F��^* := L/U_in*�� [-]
omega_a = sqrt( diag( omega_a2));

%%[2-1] ���������ψʌŗL���[�h: ��q(t)�ɑΉ�
Phi_dq_mat_BC = zeros(N_q_all,mode_num);
Phi_dq_mat_BC(not_i_vec,:) = Phi_dq_mat;

%%[2-2] ���������ʒu�ŗL���[�h: q(t) = q0 + ��q(t)�ɑΉ�
Phi_q_mat_BC = repmat( q_vec, [ 1 mode_num]) + Phi_dq_mat_BC;



%% Modal damping ratio

%%[3-0] ���[�h������ [-]
zeta_n = theta_a*omega_a/2;                      



disp( 'Nondimensional natural frequency: ')
disp( omega_a)
disp( 'Modal damping ratio: ')
disp( zeta_n)


