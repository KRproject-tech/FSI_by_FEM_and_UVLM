%% 変数抽出
X_vec = h_X_vec(:,i_time);
q_vec = X_vec(1:N_q_all,1);

%% 剛性行列の組立
generate_stiff_matrices;

%% 時間発展 (Eulerの予測子修正子法)

%%[0] 流体力の線形補間
Qf_p_global_t = (Qf_p_global - old_Qf_p_global)*(time - time_fluid)/d_t_wake + old_Qf_p_global;
h_Qf_p_global_t(:,i_time) = Qf_p_global_t;

%%[1] 付加質量マトリックスの線形補間
Qf_p_mat_global_t = (Qf_p_mat_global - old_Qf_p_mat_global)*(time - time_fluid)/d_t_wake + old_Qf_p_mat_global;


f_sys = F_sys( X_vec, M_global-Qf_p_mat_global_t, Qf_global+Qf_p_global_t, Qe_global, Qd_global, var_param);
if i_time == 1;
   old_fsys = f_sys;
end

%% 予測子の計算
X_vec_p = X_vec + d_t*( 3/2*f_sys - 1/2*old_fsys );
old_fsys = f_sys;

q_vec = X_vec_p(1:N_q_all,1);

%% 剛性行列の組立
generate_stiff_matrices;

fsys_np1 = F_sys( X_vec_p, M_global-Qf_p_mat_global_t, Qf_global+Qf_p_global_t, Qe_global, Qd_global, var_param);
new_X_vec = X_vec + d_t/2*(fsys_np1 + f_sys);


h_X_vec(:,i_time+1) = new_X_vec;

%% 発散した時は解析を停止
if sum( isnan( X_vec))

    warndlg( 'Divergence!!') 
    break;
end