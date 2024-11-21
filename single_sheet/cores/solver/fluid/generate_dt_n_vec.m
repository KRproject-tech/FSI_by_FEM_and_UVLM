%% 単位法線ベクトル算出

%%[1-0] n_vec
r13_vec_v = Sc_mat_31*q_vec;
r13_vec = reshape( r13_vec_v, 3, []).';

r42_vec_v = Sc_mat_24*q_vec;
r42_vec = reshape( r42_vec_v, 3, []).';


r13_cross_r42 = cross_mat( r13_vec, r42_vec);
n_vec_i = r13_cross_r42./norm_mat( r13_cross_r42); 

%%[1-1] dt_n_vec
dt_r13_vec_v = Sc_mat_31*dt_q_vec;
dt_r13_vec = reshape( dt_r13_vec_v, 3, []).';

dt_r42_vec_v = Sc_mat_24*dt_q_vec;
dt_r42_vec = reshape( dt_r42_vec_v, 3, []).';

%%% dt(r13 x r42) = (dt_r13 x r42) + (r13 x dt_r42)
dt_r13_cross_r42 = cross_mat( dt_r13_vec, r42_vec) + cross_mat( r13_vec, dt_r42_vec);

dt_n_vec_i = dt_r13_cross_r42./norm_mat( r13_cross_r42) - n_vec_i.*inner_mat( dt_r13_cross_r42./norm_mat( r13_cross_r42), n_vec_i);