%% [0] �\���p�����[�^�����l���

h_X_vec = zeros(4*N_q_all,length( time_m));
N_time = length( time_m);
N_fluid_time = ceil( N_time/dt_wake_per_dt);

idx_r = reshape( [1:N_qi:N_q_all; 2:N_qi:N_q_all; 3:N_qi:N_q_all], 1, []);    	%% �����`��ł̃m�[�h���W���� [-]
r0_vec = reshape( coordinates.',[],1);
h_X_vec(idx_r,1) = r0_vec;

idx_r1 = N_q_all+idx_r;                                                         %% 2���ڂ̕��ɂ��Ă��C�����l�Ɉʒu���W���� [-]
r0_vec1 = reshape( coordinates_1.',[],1);
h_X_vec(idx_r1,1) = r0_vec1;


idx_dx_r = reshape( [4:N_qi:N_q_all; 5:N_qi:N_q_all; 6:N_qi:N_q_all], 1, []); 	%% �����`��ł�x�����̌��z���� [-]
h_X_vec(idx_dx_r,1) = reshape( dx_r0_vec.', 1, []);

idx_dx_r1 = N_q_all+idx_dx_r;                                                   %% 2���ڂ̕��ɂ��Ă��C�����`��ł�x�����̌��z���� [-]
h_X_vec(idx_dx_r1,1) = reshape( dx_r0_vec_1.', 1, []);


idx_dy_r = reshape( [7:N_qi:N_q_all; 8:N_qi:N_q_all; 9:N_qi:N_q_all], 1, []); 	%% �����`��ł�y�����̌��z���� (dy_r = [0 1 0]^T [-])
idx_dy_r = [ idx_dy_r N_q_all+idx_dy_r];                                        %% 2���ڂ̕��ɂ��Ă��C�����`��ł�y�����̌��z���� (dy_r = [0 1 0]^T [-]) [-]
dy_r0_vec = repmat( [0 1 0].', [ N_node 1]);
h_X_vec(idx_dy_r,1) = repmat( dy_r0_vec, [ 2 1]);

                   



%% [1] ���̃p�����[�^�����l���
h_r_panel_vec = zeros(4*N_element,3,N_fluid_time);                            	%% �p�l���m�[�h�_�i�[ [-]
h_rcol_vec = zeros(N_element,3,N_fluid_time);                                   %% �R���P�[�V�����_�i�[ [-]
h_n_vec = h_rcol_vec;                                                           %% �P�ʖ@���x�N�g�� [-]
h_dt_n_vec = h_rcol_vec;                                                        %% �P�ʖ@���x�N�g�� [-]
h_Amat = zeros(2*N_element,2*N_element,N_fluid_time);                        	%%  A��=b [-]
h_dt_Amat = h_Amat;                                                             %%  A��=b [-]

h_V_surf = zeros(2*N_element,3,N_fluid_time);                                   %% �\�ʗ����x�N�g�� [-]
h_V_wake_end = zeros(Ny+1,3,N_fluid_time);                                      %% �㉏�����x�N�g�� [-]


h_r_end = zeros(Ny+1,3,N_fluid_time);                                           %% �㉏�m�[�h���W�x�N�g�� [-]
h_dp_add = zeros(N_element,N_fluid_time);                                       %% �t�����ʌ��� (��_n+1 - ��_n)/��t [-]
h_dp_add_estimate = zeros(2*N_element,N_fluid_time);                          	%% �t�����ʌ��� [-]
h_dp_lift = h_dp_add_estimate;                                                  %% ��헬�̗͌��� [-]
h_dp_vec = h_dp_add_estimate;                                                   %% ���̗� [-]
h_Qf_p_global_t = zeros(N_q_all,N_time);                                        %% ���̗̓}�g���b�N�X [-]
h_Qf_p_global_t_1 = h_Qf_p_global_t;

old_Qf_p_global = zeros(N_q_all,1);
old_Qf_p_global_1 = old_Qf_p_global;
old_Qf_p_mat_global = zeros(N_q_all,2*N_q_all);
old_Qf_p_mat_global_1 = old_Qf_p_mat_global;
old_Qf_p_mat0_global = zeros(N_q_all,2*N_element);
old_Qf_p_mat0_global_1 = old_Qf_p_mat0_global;
old_Qf_p_lift2_mat_global = zeros(N_q_all,3*N_element);
old_Qf_p_lift2_mat_global_1 = old_Qf_p_lift2_mat_global;

Qf_p_global_a = old_Qf_p_global;
Qf_p_global_a_1 = old_Qf_p_global;
Qf_p_mat_global_a = old_Qf_p_mat_global;
Qf_p_mat_global_a_1 = old_Qf_p_mat_global;
Qf_p_mat0_global_a = old_Qf_p_mat0_global;
Qf_p_mat0_global_a_1 = old_Qf_p_mat0_global;
Qf_p_lift2_mat_global_a = old_Qf_p_lift2_mat_global;
Qf_p_lift2_mat_global_a_1 = old_Qf_p_lift2_mat_global_1;

Qf_p_global = old_Qf_p_global;
Qf_p_global_1 = old_Qf_p_global;
Qf_p_mat_global = old_Qf_p_mat_global;
Qf_p_mat_global_1 = old_Qf_p_mat_global;
Qf_p_mat0_global = old_Qf_p_mat0_global;
Qf_p_mat0_global_1 = old_Qf_p_mat0_global;
Qf_p_lift2_mat_global = old_Qf_p_lift2_mat_global;
Qf_p_lift2_mat_global_1 = old_Qf_p_lift2_mat_global_1;


%% [2] �G�l���M�E�d����������
h_E_inertia = zeros(2,length( time_m));
h_E_em = h_E_inertia;
h_E_ek = h_E_inertia;


h_W_dm = h_E_inertia;
h_W_dk = h_E_inertia;
h_W_em2 = h_E_inertia;
h_W_ek2 = h_E_inertia;


h_W_f_ext = zeros(2,N_fluid_time);
h_W_f_ext_Xdist = zeros(N_fluid_time,Nx);