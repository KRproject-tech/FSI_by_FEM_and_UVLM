%% Οo
X_vec = h_X_vec(:,i_time);
q_vec = X_vec(1:N_q_all,1);

%% «sρΜg§
generate_stiff_matrices;

%% Τ­W (EulerΜ\ͺqC³q@)

%%[0] ¬ΜΝΜό`βΤ
Qf_p_global_t = (Qf_p_global - old_Qf_p_global)*(time - time_fluid)/d_t_wake + old_Qf_p_global;
h_Qf_p_global_t(:,i_time) = Qf_p_global_t;

%%[1] tΑΏΚ}gbNXΜό`βΤ
Qf_p_mat_global_t = (Qf_p_mat_global - old_Qf_p_mat_global)*(time - time_fluid)/d_t_wake + old_Qf_p_mat_global;


f_sys = F_sys( X_vec, M_global-Qf_p_mat_global_t, Qf_global+Qf_p_global_t, Qe_global, Qd_global, var_param);
if i_time == 1;
   old_fsys = f_sys;
end

%% \ͺqΜvZ
X_vec_p = X_vec + d_t*( 3/2*f_sys - 1/2*old_fsys );
old_fsys = f_sys;

q_vec = X_vec_p(1:N_q_all,1);

%% «sρΜg§
generate_stiff_matrices;

fsys_np1 = F_sys( X_vec_p, M_global-Qf_p_mat_global_t, Qf_global+Qf_p_global_t, Qe_global, Qd_global, var_param);
new_X_vec = X_vec + d_t/2*(fsys_np1 + f_sys);


h_X_vec(:,i_time+1) = new_X_vec;

%% ­U΅½ΝπΝπβ~
if sum( isnan( X_vec))

    warndlg( 'Divergence!!') 
    break;
end