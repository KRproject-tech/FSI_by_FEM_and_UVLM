%% �ϐ����o
X_vec = h_X_vec(:,i_time);
q_vec = X_vec(1:N_q_all,1);
dt_q_vec = X_vec(N_q_all+1:end,1);





%% ���Ԕ��W (Euler�̗\���q�C���q�@)

%%[*] ���̗̓x�N�g���F[p] = p_lift + Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2) 

%%[0] ���̗͂̐��`��� (Mf2_2)
Qf_p_global_t = @( time) (Qf_p_global - old_Qf_p_global)*(time - time_fluid)/d_t_wake + Qf_p_global_a;

Qf_p_global_tv_n =  Qf_p_global_t( time);
h_Qf_p_global_t(:,i_time) = Qf_p_global_tv_n;


%%[1] ( -(��x*dx�� + ��y*dy��)*dt_rc ) 
Qf_p_lift2_global_t = @( time) (Qf_p_lift2_mat_global - old_Qf_p_lift2_mat_global)*(time - time_fluid)/d_t_wake + Qf_p_lift2_mat_global_a;

Qf_p_lift2_global_tv_n = Qf_p_lift2_global_t( time);


%%[2] �t�����ʃ}�g���b�N�X�̐��`��� (Mf1)
Qf_p_mat_global_t = @( time) (Qf_p_mat_global - old_Qf_p_mat_global)*(time - time_fluid)/d_t_wake + Qf_p_mat_global_a;

Qf_p_mat_global_tv_n = Qf_p_mat_global_t( time);

%%[3] ���̗͂̐��`��� (Mf2_1)
Qf_p_mat0_global_t = @( time) (Qf_p_mat0_global - old_Qf_p_mat0_global)*(time - time_fluid)/d_t_wake + Qf_p_mat0_global_a;

Qf_p_mat0_global_tv_n = Qf_p_mat0_global_t( time);


%% ���̗͌v�Z
%%[*] �P�ʖ@���x�N�g���Z�o
generate_dt_n_vec;


%%[*] �R���P�[�V�����_�̕ψʑ��x

rc_vec = Sc_mat_col_global*q_vec;
rc_vec = reshape( rc_vec, 3, []).';

dt_rc_vec = Sc_mat_col_global*dt_q_vec;
dt_rc_vec = reshape( dt_rc_vec, 3, []).';





if ~exist( 'Gamma_wake', 'var')
    
   V_wake_plate = 0; 
   dt_Amat2_Gamma = 0;
   dt_Amat1 = 0;
   Gamma = 0;
end



Qf_p_mat0_global_t_n = Qf_p_mat0_global_tv_n*( sum( (dt_rc_vec - V_in - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i, 2) - dt_Amat1*Gamma);
Qf_p_lift2_global_t_n = Qf_p_lift2_global_tv_n*reshape( dt_rc_vec.', [], 1); 

%% �\���q�v�Z (�V���������̏�Ԃ����߂�)


%%[0] �S�e���x�N�g���̑g�� Qe^(n)
flag_output = 1;                                                            %% �����s��Z�o�̗L����
generate_stiff_matrices;

Qk_global_n = Qk_global;
Qe_global_n = Qe_global;
dq_Qe_global_n = dq_Qe_global;
Qd_global_n = Qd_global;
Qd_theta_global_n = Qd_theta_global;
J_global_1_n = J_global_1;
J_global_2_n = J_global_2;




%%[1] �\���q�v�Z
%%[1-0] ���ʍs��     : M
m_global_struct.M_global = M_global + J_global_1_n - Qf_p_mat_global_tv_n;
%%[1-1] �̐ϗ�       : Q_f
qf_global_struct.Qf_global = Qf_global + Qf_time_global*q_in_norm( time) + Qf_p_global_tv_n + Qf_p_mat0_global_t_n + Qf_p_lift2_global_t_n;
%%[1-2] dq_Qe(q(n))
dq_qe_global_struct.dq_Qe_global = dq_Qe_global_n;
%%[1-3] ������       : Q_e
qe_global_struct.Qe_global = Qe_global_n + Qk_global_n;
%%[1-4] ������       : Q_d
qd_global_struct.Qd_global = Qd_global_n + Qd_theta_global_n + J_global_2_n;


[ X_vec_p, out1] = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, 0, []);

                                

                                
                                
                                
                                

%% �C���q�v�Z (�V���������̍����� Qe^(n+1) �̌��ŉ���)

q_vec = X_vec_p(1:N_q_all,1);
dt_q_vec = X_vec_p(N_q_all+1:end,1);


%%[*] �P�ʖ@���x�N�g���Z�o
generate_dt_n_vec;


%%[*] �R���P�[�V�����_�̕ψʑ��x

dt_rc_vec = Sc_mat_col_global*dt_q_vec;
dt_rc_vec = reshape( dt_rc_vec, 3, []).';

%%%(�C���O) Ma=1, Ua25, Nx=18, Ny=10 �ŃJ�I�X�U��������
%%[*] t+dt�ł̗��̗͂�\��
%%[*] ���̗̓x�N�g���F[p] = p_lift + Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2) 
% Qf_p_global_tv_np1 =  Qf_p_global_t( time + d_t);
%%[1] ( -(��x*dx�� + ��y*dy��)*dt_rc ) 
% Qf_p_lift2_global_tv_np1 = Qf_p_lift2_global_t( time + d_t);
%%[2] �t�����ʃ}�g���b�N�X�̐��`��� (Mf1)
% Qf_p_mat_global_tv_np1 = Qf_p_mat_global_t( time + d_t);�@%% Ma=1, Ua25, Nx=18, Ny=10 �ŃJ�I�X�U��������
%%[3] ���̗͂̐��`��� (Mf2_1)
% Qf_p_mat0_global_tv_np1 = Qf_p_mat0_global_t( time + d_t);

%%%(�C���O) Ma=1, Ua25, Nx=18, Ny=10 �ŃJ�I�X�U��������
% Qf_p_mat0_global_t_np1 = Qf_p_mat0_global_tv_np1*( sum( (dt_rc_vec - V_in - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i, 2) - dt_Amat1*Gamma); 
% Qf_p_lift2_global_t_np1 = Qf_p_lift2_global_tv_np1*reshape( dt_rc_vec.', [], 1); 
Qf_p_mat0_global_t_np1 = Qf_p_mat0_global_tv_n*( sum( (dt_rc_vec - V_in - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i, 2) - dt_Amat1*Gamma);
Qf_p_lift2_global_t_np1 = Qf_p_lift2_global_tv_n*reshape( dt_rc_vec.', [], 1); 



%%[0] �S�e���x�N�g���̑g�� Qe^(n+1)
flag_output = 0;                                                            %% �����s��Z�o�̖�����
generate_stiff_matrices;
Qk_global_np1 = Qk_global;
Qd_theta_global_np1 = Qd_theta_global;



%%[1-0] ���ʍs��     : M
% m_global_struct.M_global = M_global + J_global_1_n - (
% Qf_p_mat_global_tv_n + Qf_p_mat_global_tv_np1 )/2;   %% Ma=1, Ua25, Nx=18, Ny=10 �ŃJ�I�X�U��������
m_global_struct.M_global = M_global + J_global_1_n - Qf_p_mat_global_tv_n;

%%%(�C���O) Ma=1, Ua25, Nx=18, Ny=10 �ŃJ�I�X�U��������
%%[1-1] �̐ϗ�       : Q_f
% qf_global_struct.Qf_global = Qf_global + Qf_time_global*q_in_norm( time) ...
%                                         + ( Qf_p_global_tv_n + Qf_p_global_tv_np1 )/2 ...
%                                         + ( Qf_p_mat0_global_t_n + Qf_p_mat0_global_t_np1 )/2 ...
%                                         + ( Qf_p_lift2_global_t_n + Qf_p_lift2_global_t_np1 )/2;
qf_global_struct.Qf_global = Qf_global + Qf_time_global*q_in_norm( time) ...
                                        + Qf_p_global_tv_n ...
                                        + ( Qf_p_mat0_global_t_n + Qf_p_mat0_global_t_np1 )/2 ...
                                        + ( Qf_p_lift2_global_t_n + Qf_p_lift2_global_t_np1 )/2;
%%[1-2] dq_Qe(q(n))
dq_qe_global_struct.dq_Qe_global = dq_Qe_global_n;
%%[1-3] ������       : Q_e
%%[*] Qe(q(n+1)) = Qe(q(n)) + ��t*dq_Qe(q(n))*dt_q(n+1)
qe_global_struct.Qe_global = Qe_global_n + (Qk_global_n + Qk_global_np1)/2;
%%[1-4] ������       : Q_d
qd_global_struct.Qd_global = Qd_global_n + Qd_theta_global_n + J_global_2_n;


%%[*] Qe(q(n+1)) = Qe(q(n)) + ��t*dq_Qe(q(n))*dt_q(n+1)
new_X_vec = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, 1, out1);

h_X_vec(:,i_time+1) = new_X_vec;                                            %% (Qe^(n)+Qe^(n+1))/2�̌��ŉ����� X(n+1) 


%% ���U�������͉�͂��~
if sum( isnan( X_vec))

    warndlg( 'Divergence!!') 
%     break;
end
