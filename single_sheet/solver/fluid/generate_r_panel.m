%% �p�l���m�[�h�_�v�Z [m]


%%[0-0] �_�@ [-]
r_panel_vec_1v = Sc_mat_panel_global_1*q_vec;
r_panel_vec_1 = reshape( r_panel_vec_1v, 3, []).';

%%[0-1] �_�A [-]�@(�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
r_panel_vec_2v = Sc_mat_panel_global_2*q_vec;
r_panel_vec_2 = reshape( r_panel_vec_2v, 3, []).';

%%[0-2] �_�B [-] (�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
r_panel_vec_3v = Sc_mat_panel_global_3*q_vec;
r_panel_vec_3 = reshape( r_panel_vec_3v, 3, []).';

%%[0-3] �_�C [-]
r_panel_vec_4v = Sc_mat_panel_global_4*q_vec;
r_panel_vec_4 = reshape( r_panel_vec_4v, 3, []).';



%%[0-4] �_�@ [-]
dt_r_panel_vec_1v = Sc_mat_panel_global_1*dt_q_vec;
dt_r_panel_vec_1 = reshape( dt_r_panel_vec_1v, 3, []).';

%%[0-5] �_�A [-]�@(�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
dt_r_panel_vec_2v = Sc_mat_panel_global_2*dt_q_vec;
dt_r_panel_vec_2 = reshape( dt_r_panel_vec_2v, 3, []).';

%%[0-6] �_�B [-] (�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
dt_r_panel_vec_3v = Sc_mat_panel_global_3*dt_q_vec;
dt_r_panel_vec_3 = reshape( dt_r_panel_vec_3v, 3, []).';

%%[0-7] �_�C [-]
dt_r_panel_vec_4v = Sc_mat_panel_global_4*dt_q_vec;
dt_r_panel_vec_4 = reshape( dt_r_panel_vec_4v, 3, []).';