function [ v_wake, q1234_wake_mat] = V_wake_func( r_vec, r_vec_1, r_vec_2, r_vec_3, r_vec_4, Gamma, var_param, fine_flag)



q1234_wake_mat = generate_q1234_mat( r_vec, r_vec_1, r_vec_2, r_vec_3, r_vec_4, var_param, fine_flag);   
N_r_vec = size( r_vec, 1);
q_gamma_wake = q1234_wake_mat.*( ones(N_r_vec,1)*kron( Gamma.', ones(1,3)) );    
v_wake = [ sum( q_gamma_wake(:,1:3:end), 2) sum( q_gamma_wake(:,2:3:end), 2) sum( q_gamma_wake(:,3:3:end), 2)];


end