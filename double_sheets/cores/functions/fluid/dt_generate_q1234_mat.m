function dt_q1234_mat = dt_generate_q1234_mat(  rc_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, ...
                                                dt_rc_vec, dt_r_panel_vec_1, dt_r_panel_vec_2, dt_r_panel_vec_3, dt_r_panel_vec_4)


N_element1 = size( r_panel_vec_1, 1);
N_element2 = size( rc_vec, 1);


    





%% Generation of influence coefficient matrix
%%%
%%%    ^ Y
%%%    | �@-----2--------�A
%%% 　 | |               |
%%%    | 1　　　 X　　　　2
%%%    | |               |
%%%    | �C-----4--------�B
%%%    |--------------------------> X
%%%

rc_mat = repmat( rc_vec, [ 1 N_element1]);
r_panel_mat_1 = ones(N_element2,1)*reshape( r_panel_vec_1.', 1, []);
r_panel_mat_2 = ones(N_element2,1)*reshape( r_panel_vec_2.', 1, []);
r_panel_mat_3 = ones(N_element2,1)*reshape( r_panel_vec_3.', 1, []);
r_panel_mat_4 = ones(N_element2,1)*reshape( r_panel_vec_4.', 1, []);

dt_rc_mat = repmat( dt_rc_vec, [ 1 N_element1]);
dt_r_panel_mat_1 = ones(N_element2,1)*reshape( dt_r_panel_vec_1.', 1, []);
dt_r_panel_mat_2 = ones(N_element2,1)*reshape( dt_r_panel_vec_2.', 1, []);
dt_r_panel_mat_3 = ones(N_element2,1)*reshape( dt_r_panel_vec_3.', 1, []);
dt_r_panel_mat_4 = ones(N_element2,1)*reshape( dt_r_panel_vec_4.', 1, []);




%%[1-0] コロケーション点とパネルノード点の間の距離ベクトル [-]
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


%%[1-1] コロケーション点とパネルノード点の間の速度ベクトル [-]
dt_r11_mat = dt_rc_mat - dt_r_panel_mat_1;
dt_r12_mat = dt_rc_mat - dt_r_panel_mat_2;
dt_r13_mat = dt_rc_mat - dt_r_panel_mat_3;
dt_r14_mat = dt_rc_mat - dt_r_panel_mat_4;

dt_r21_mat = dt_rc_mat - dt_r_panel_mat_4;
dt_r22_mat = dt_rc_mat - dt_r_panel_mat_1;
dt_r23_mat = dt_rc_mat - dt_r_panel_mat_2;
dt_r24_mat = dt_rc_mat - dt_r_panel_mat_3;


dt_r01_mat = dt_r11_mat - dt_r21_mat;
dt_r02_mat = dt_r12_mat - dt_r22_mat;
dt_r03_mat = dt_r13_mat - dt_r23_mat;
dt_r04_mat = dt_r14_mat - dt_r24_mat;



%%[1-2] 誘導速度算出 

r11_cross_r21 = cross_mat( r11_mat, r21_mat);
norm_r11_cross_r21 = norm_mat( r11_cross_r21);
dt_r11_cross_r21 = cross_mat( dt_r11_mat, r21_mat) + cross_mat( r11_mat, dt_r21_mat);
dt_r11_cross_r21_norm_r11_cross_r21 = dt_r11_cross_r21./norm_r11_cross_r21.^2 - 2*r11_cross_r21.*inner_mat( dt_r11_cross_r21, r11_cross_r21)./norm_r11_cross_r21.^4;
norm_r11 = norm_mat( r11_mat);
norm_r21 = norm_mat( r21_mat);
r11_per_norm_r11_m_r21_mat_per_norm_r21 = r11_mat./norm_r11 - r21_mat./norm_r21;
r11_cross_r21_per_norm_r11_cross_r21_2 = r11_cross_r21./( norm_r11_cross_r21.^2 );
%%% dt_q1
dt_q1_mat = dt_r11_cross_r21_norm_r11_cross_r21.*inner_mat( r01_mat, r11_per_norm_r11_m_r21_mat_per_norm_r21 ) ...
            + r11_cross_r21_per_norm_r11_cross_r21_2.*inner_mat( dt_r01_mat, r11_per_norm_r11_m_r21_mat_per_norm_r21 ) ...
            + r11_cross_r21_per_norm_r11_cross_r21_2.*inner_mat( r01_mat, dt_r11_mat./norm_r11 - r11_mat.*inner_mat( dt_r11_mat, r11_mat)./norm_r11.^3 ...
                                                                            - dt_r21_mat./norm_r21 + r21_mat.*inner_mat( dt_r21_mat, r21_mat)./norm_r21.^3 );



r12_cross_r22 = cross_mat( r12_mat, r22_mat);
norm_r12_cross_r22 = norm_mat( r12_cross_r22);
dt_r12_cross_r22 = cross_mat( dt_r12_mat, r22_mat) + cross_mat( r12_mat, dt_r22_mat);
dt_r12_cross_r22_norm_r12_cross_r22 = dt_r12_cross_r22./norm_r12_cross_r22.^2 - 2*r12_cross_r22.*inner_mat( dt_r12_cross_r22, r12_cross_r22)./norm_r12_cross_r22.^4;
norm_r12 = norm_mat( r12_mat);
norm_r22 = norm_mat( r22_mat);
r12_per_norm_r12_m_r22_per_norm_r22 = r12_mat./norm_r12 - r22_mat./norm_r22;
r12_cross_r22_per_norm_r12_cross_r22_2 = r12_cross_r22./( norm_r12_cross_r22.^2 );
%%% dt_q2
dt_q2_mat = dt_r12_cross_r22_norm_r12_cross_r22.*inner_mat( r02_mat, r12_per_norm_r12_m_r22_per_norm_r22 ) ...
            + r12_cross_r22_per_norm_r12_cross_r22_2.*inner_mat( dt_r02_mat, r12_per_norm_r12_m_r22_per_norm_r22 ) ...
            + r12_cross_r22_per_norm_r12_cross_r22_2.*inner_mat( r02_mat, dt_r12_mat./norm_r12 - r12_mat.*inner_mat( dt_r12_mat, r12_mat)./norm_r12.^3 ...
                                                                            - dt_r22_mat./norm_r22 + r22_mat.*inner_mat( dt_r22_mat, r22_mat)./norm_r22.^3 );




r13_cross_r23 = cross_mat( r13_mat, r23_mat);
norm_r13_cross_r23 = norm_mat( r13_cross_r23);
dt_r13_cross_r23 = cross_mat( dt_r13_mat, r23_mat) + cross_mat( r13_mat, dt_r23_mat);
dt_r13_cross_r23_norm_r13_cross_r23 = dt_r13_cross_r23./norm_r13_cross_r23.^2 - 2*r13_cross_r23.*inner_mat( dt_r13_cross_r23, r13_cross_r23)./norm_r13_cross_r23.^4;
norm_r13 = norm_mat( r13_mat);
norm_r23 = norm_mat( r23_mat);
r13_per_norm_r13_m_r23_per_norm_r23 = r13_mat./norm_r13 - r23_mat./norm_r23;
r13_cross_r23_norm_r13_cross_r23_2 = r13_cross_r23./( norm_r13_cross_r23.^2 );
%%% dt_q3
dt_q3_mat = dt_r13_cross_r23_norm_r13_cross_r23.*inner_mat( r03_mat, r13_per_norm_r13_m_r23_per_norm_r23 ) ...
            + r13_cross_r23_norm_r13_cross_r23_2.*inner_mat( dt_r03_mat, r13_per_norm_r13_m_r23_per_norm_r23 ) ...
            + r13_cross_r23_norm_r13_cross_r23_2.*inner_mat( r03_mat, dt_r13_mat./norm_r13 - r13_mat.*inner_mat( dt_r13_mat, r13_mat)./norm_r13.^3 ...
                                                                            - dt_r23_mat./norm_r23 + r23_mat.*inner_mat( dt_r23_mat, r23_mat)./norm_r23.^3 );




r14_cross_r24 = cross_mat( r14_mat, r24_mat);
norm_r14_cross_r24 = norm_mat( r14_cross_r24);
dt_r14_cross_r24 = cross_mat( dt_r14_mat, r24_mat) + cross_mat( r14_mat, dt_r24_mat);
dt_r14_cross_r24_norm_r14_cross_r24 = dt_r14_cross_r24./norm_r14_cross_r24.^2 - 2*r14_cross_r24.*inner_mat( dt_r14_cross_r24, r14_cross_r24)./norm_r14_cross_r24.^4;
norm_r14 = norm_mat( r14_mat);
norm_r24 = norm_mat( r24_mat);
r14_per_norm_r14_m_r24_per_norm_r24 = r14_mat./norm_r14 - r24_mat./norm_r24;
r14_cross_r24_per_norm_r14_cross_r24_2 = r14_cross_r24./( norm_r14_cross_r24.^2 );
%%% dt_q4
dt_q4_mat = dt_r14_cross_r24_norm_r14_cross_r24.*inner_mat( r04_mat, r14_per_norm_r14_m_r24_per_norm_r24 ) ...
            + r14_cross_r24_per_norm_r14_cross_r24_2.*inner_mat( dt_r04_mat, r14_per_norm_r14_m_r24_per_norm_r24 ) ...
            + r14_cross_r24_per_norm_r14_cross_r24_2.*inner_mat( r04_mat, dt_r14_mat./norm_r14 - r14_mat.*inner_mat( dt_r14_mat, r14_mat)./norm_r14.^3 ...
                                                                            - dt_r24_mat./norm_r24 + r24_mat.*inner_mat( dt_r24_mat, r24_mat)./norm_r24.^3 );



%%[1-2] Generation of influence coefficient matrix
dt_q1234_mat = -1/(4*pi)*( dt_q1_mat + dt_q2_mat + dt_q3_mat + dt_q4_mat );

end