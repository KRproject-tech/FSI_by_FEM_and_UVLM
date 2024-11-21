%% 変数抽出 (初期形状)
X_vec = h_X_vec(:,1);
q_vec_0 = X_vec(1:N_q_all,1);



%% Rigid motion

theta_pitch_vec = theta_pitch_time( time);
dt_theta_pitch_vec = dt_theta_pitch_time( time);

R_pitch_mat = sparse( R_pitch( theta_pitch_vec));
str_mat = repmat( 'R_pitch_mat,', [ 1 3*N_node]);
eval( [ 'R_pitch_mat_global = blkdiag( ', str_mat(1:end-1), ');']);

dt_R_pitch_mat = sparse( dt_R_pitch( theta_pitch_vec, dt_theta_pitch_vec));
str_mat = repmat( 'dt_R_pitch_mat,', [ 1 3*N_node]);
eval( [ 'dt_R_pitch_mat_global = blkdiag( ', str_mat(1:end-1), ');']);

q_vec_offset = zeros(N_qi,1);
q_vec_offset(1) = L_pitch_center;
q_vec_offset_global = repmat( q_vec_offset, [ N_node 1]);




%%%
%%% q_G(t) = R(θ(t)) q_G(0) 
%%% dt_q_G(t) = dt_θ(t) dθ_R(θ(t)) q_G(0) 
%%%
q_vec = R_pitch_mat_global*(q_vec_0 - q_vec_offset_global) + q_vec_offset_global;
dt_q_vec = dt_R_pitch_mat_global*(q_vec_0 - q_vec_offset_global);


%% Update variables

new_X_vec = X_vec;

new_X_vec(1:N_q_all,1) = q_vec;
new_X_vec(N_q_all+1:end,1) = dt_q_vec;


h_X_vec(:,i_time+1) = new_X_vec;                                            %% (Qe^(n)+Qe^(n+1))/2の元で解いた X(n+1) 

%% 流体力算用にX_vecを更新
X_vec = new_X_vec;