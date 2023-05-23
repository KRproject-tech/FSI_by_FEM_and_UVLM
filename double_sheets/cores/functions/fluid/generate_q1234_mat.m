function q1234_mat = generate_q1234_mat(  rc_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, var_param, fine_flag)


N_element1 = size( r_panel_vec_1, 1);
N_element2 = size( rc_vec, 1);


r_eps = var_param.r_eps;
Ncore = var_param.Ncore;
eps_v = var_param.eps_v;

    
Length = var_param.Length;
Nx = var_param.Nx;




%% Generation of influence coefficient matrix
%%%
%%%    ^ Y
%%%    | ①-----2--------②
%%% 　 | |               |
%%%    | 1　　　 X　　　　2
%%%    | |               |
%%%    | ④-----4--------③
%%%    |--------------------------> X
%%%

rc_mat = repmat( rc_vec, [ 1 N_element1]);
r_panel_mat_1 = ones(N_element2,1)*reshape( r_panel_vec_1.', 1, []);
r_panel_mat_2 = ones(N_element2,1)*reshape( r_panel_vec_2.', 1, []);
r_panel_mat_3 = ones(N_element2,1)*reshape( r_panel_vec_3.', 1, []);
r_panel_mat_4 = ones(N_element2,1)*reshape( r_panel_vec_4.', 1, []);



%%[1-0] コロケーション点とパネルノード点の間の距離ベクトル [m]
r11_mat = rc_mat - r_panel_mat_1;
r12_mat = rc_mat - r_panel_mat_2;
r13_mat = rc_mat - r_panel_mat_3;
r14_mat = rc_mat - r_panel_mat_4;

r21_mat = rc_mat - r_panel_mat_4;
r22_mat = rc_mat - r_panel_mat_1;
r23_mat = rc_mat - r_panel_mat_2;
r24_mat = rc_mat - r_panel_mat_3;


r01_mat = r11_mat - r21_mat;
r02_mat = r12_mat - r22_mat;
r03_mat = r13_mat - r23_mat;
r04_mat = r14_mat - r24_mat;



%%[1-1] 誘導速度算出 


%%%[*] Vortex core サイズはWakeの各格子長の最大値を基準に設定
norm_r01 = norm_mat( r01_mat);
norm_r02 = norm_mat( r02_mat);
norm_r03 = norm_mat( r03_mat);
norm_r04 = norm_mat( r04_mat);

max_wake_r0 = max( norm_r01(:,1:3:end), norm_r02(:,1:3:end));
max_wake_r1 = max( max_wake_r0, norm_r03(:,1:3:end));
max_wake_r = max( max_wake_r1, norm_r04(:,1:3:end));

max_wake_r = max( max_wake_r, Length/Nx);
    
r_core = (fine_flag == 1)*max_wake_r.*r_eps.fine + ~(fine_flag == 1)*max_wake_r.*r_eps.rough;



r11_cross_r21 = cross_mat( r11_mat, r21_mat);
norm_r11_cross_r21 = norm_mat( r11_cross_r21);
norm_r11 = norm_mat( r11_mat);
norm_r21 = norm_mat( r21_mat);
q1_mat = 1/(4*pi)*r11_cross_r21./( norm_r11_cross_r21.^2 + eps_v).*inner_mat( r01_mat, r11_mat./max( norm_r11, eps) - r21_mat./max( norm_r21, eps) );
% r11_r21 = max( norm_r11, eps).*max( norm_r21, eps);
% q1_mat = 1/(4*pi)*(norm_r11 + norm_r21).*r11_cross_r21./max( inner_mat( r11_mat, r21_mat) + r11_r21, eps^2)./r11_r21;
%%[1-1-0] Vortex core model
h1 = norm_r11_cross_r21./norm_r01;                              %% the perpendicular distance
h1 = h1(:,1:3:end);                                             %% 全成分に対してKvを計算するより速い． 
Kv1 = h1.^2./( h1.^(2*Ncore) + r_core.^(2*Ncore) ).^(1/Ncore);
Kv1 = kron( Kv1, ones(1,3));
q1_mat = Kv1.*q1_mat;


r12_cross_r22 = cross_mat( r12_mat, r22_mat);
norm_r12_cross_r22 = norm_mat( r12_cross_r22);
norm_r12 = norm_mat( r12_mat);
norm_r22 = norm_mat( r22_mat);
q2_mat = 1/(4*pi)*r12_cross_r22./( norm_r12_cross_r22.^2 + eps_v).*inner_mat( r02_mat, r12_mat./max( norm_r12, eps) - r22_mat./max( norm_r22, eps) );
% r12_r22 = max( norm_r12, eps).*max( norm_r22, eps);
% q2_mat = 1/(4*pi)*(norm_r12 + norm_r22).*r12_cross_r22./max( inner_mat( r12_mat, r22_mat) + r12_r22, eps^2)./r12_r22;
%%[1-1-1] Vortex core model
h2 = norm_r12_cross_r22./norm_r02;                              %% the perpendicular distance
h2 = h2(:,1:3:end);                                             %% 全成分に対してKvを計算するより速い． 
Kv2 = h2.^2./( h2.^(2*Ncore) + r_core.^(2*Ncore) ).^(1/Ncore);
Kv2 = kron( Kv2, ones(1,3));
q2_mat = Kv2.*q2_mat;



r13_cross_r23 = cross_mat( r13_mat, r23_mat);
norm_r13_cross_r23 = norm_mat( r13_cross_r23);
norm_r13 = norm_mat( r13_mat);
norm_r23 = norm_mat( r23_mat);
q3_mat = 1/(4*pi)*r13_cross_r23./( norm_r13_cross_r23.^2 + eps_v).*inner_mat( r03_mat, r13_mat./max( norm_r13, eps) - r23_mat./max( norm_r23, eps) );
% r13_r23 = max( norm_r13, eps).*max( norm_r23, eps);
% q3_mat = 1/(4*pi)*(norm_r13 + norm_r23).*r13_cross_r23./max( inner_mat( r13_mat, r23_mat) + r13_r23, eps^2)./r13_r23;
%%[1-1-2] Vortex core model
h3 = norm_r13_cross_r23./norm_r03;                              %% the perpendicular distance
h3 = h3(:,1:3:end);                                             %% 全成分に対してKvを計算するより速い．     
Kv3 = h3.^2./( h3.^(2*Ncore) + r_core.^(2*Ncore) ).^(1/Ncore);
Kv3 = kron( Kv3, ones(1,3));
q3_mat = Kv3.*q3_mat;



r14_cross_r24 = cross_mat( r14_mat, r24_mat);
norm_r14_cross_r24 = norm_mat( r14_cross_r24);
norm_r14 = norm_mat( r14_mat);
norm_r24 = norm_mat( r24_mat);
q4_mat = 1/(4*pi)*r14_cross_r24./( norm_r14_cross_r24.^2 + eps_v).*inner_mat( r04_mat, r14_mat./max( norm_r14, eps) - r24_mat./max( norm_r24, eps) );
% r14_r24 = max( norm_r14, eps).*max( norm_r24, eps);
% q4_mat = 1/(4*pi)*(norm_r14 + norm_r24).*r14_cross_r24./max( inner_mat( r14_mat, r24_mat) + r14_r24, eps^2)./r14_r24;
%%[1-1-3] Vortex core model
h4 = norm_r14_cross_r24./norm_r04;                              %% the perpendicular distance
h4 = h4(:,1:3:end);                                             %% 全成分に対してKvを計算するより速い． 
Kv4 = h4.^2./( h4.^(2*Ncore) + r_core.^(2*Ncore) ).^(1/Ncore);
Kv4 = kron( Kv4, ones(1,3));
q4_mat = Kv4.*q4_mat;



%%[1-2] Generation of influence coefficient matrix
q1234_mat = -(q1_mat + q2_mat + q3_mat + q4_mat);

end