function [ out, out2] = dq_n_dx_r_dq_theta_0y( q_i_vec, dx_n_Sc_struct, p_vec)

%% x = 0での計算（柔軟平板の前縁）

N_q = length( q_i_vec);                                             %% 1要素当たりのノードの成分数
lgth_p = length( p_vec);


%% 勾配計算

%%[*] 0成分は除く (Sc*q = [ S1*I S2*I ... S12*I]*q = S1*q1_r + S2*q1_dx_r + ... + S12*q4_dy_r, I∈R^3*3, qi_r,qi_dx_r,qi_dy_r∈R^3)
q_i_vec = q_i_vec(:,1,1,ones(1,lgth_p));
q_i_vec = reshape( q_i_vec, 3, [], 1, lgth_p);                      %% 形状関数の行数に合わせる(積の計算のため)．1方向:座標成分(x,y,z)に対応，2:方向(r,dx_r,dy_r)に対応，3,4方向:ガウス求積点


[ n_vec_norm_n, n_vec, norm_n] = n_vec_norm_n_f( q_i_vec, dx_n_Sc_struct);

[ dq_n_norm_n_dx_r_norm_dx_r_dq_theta_0y_v, dq_theta_0y_v] = dq_n_norm_n_dx_r_norm_dx_r_dq_theta_0y( q_i_vec, n_vec_norm_n, dx_n_Sc_struct, N_q, norm_n);

out = dq_n_norm_n_dx_r_norm_dx_r_dq_theta_0y_v;
out2 = dq_theta_0y_v;

end



%% 法線ベクトル (n = dx_r×dy_r)
function out = n_vec_f( q_vec, dx_n_Sc_struct)

dx_Sc_mat_v_x0 = dx_n_Sc_struct.dx_Sc_mat_v_x0;
dy_Sc_mat_v_x0 = dx_n_Sc_struct.dy_Sc_mat_v_x0;


  
dx_r = mntimes2_fast( dx_Sc_mat_v_x0, q_vec);
dy_r = mntimes2_fast( dy_Sc_mat_v_x0, q_vec);

    
out =  cross_fast( dx_r, dy_r);

end

%% n/||n||^1
function [ out, n_vec, norm_n_3] = n_vec_norm_n_f( q_vec, dx_n_Sc_struct)


n_vec = n_vec_f( q_vec, dx_n_Sc_struct);

%%[*] n/||n||^3は誤り？
%%% (Hui Wan et al., Study of Strain Energy in Deformed Insect Wings, Dynamic Behavior of Materials, 
%%%  Proceedings of the 2011 AnnualConference on Experimental and Applied
%%%  Mechanics, Vol. 1, 6 pages.)
norm_n_3 = norm2( n_vec).^1;

out = n_vec./norm_n_3(ones(3,1),:,:,:);

end


%% norm
function out = norm2( a)

out = sqrt( sum( a.^2, 1));

end


%% 外積
function out = cross_fast( a, b)


out = a([2 3 1],:,:,:).*b([3 1 2],:,:,:) - a([3 1 2],:,:,:).*b([2 3 1],:,:,:);

end


%% dq_(n/||n||)^T*dx_r/||dx_r||*∂q_θ0y
function [ out, out2] = dq_n_norm_n_dx_r_norm_dx_r_dq_theta_0y( q_vec, n_vec_norm_n, dx_n_Sc_struct, N_q, norm_n)

%%[*] 0成分を含んだもの
dx_Sc_mat_v_o = dx_n_Sc_struct.dx_Sc_mat_v_x0_o;
dy_Sc_mat_v_o = dx_n_Sc_struct.dy_Sc_mat_v_x0_o;

%%[*] 行列の積演算の高速化のために0成分を除いたもの (d_xi_S(0,y*))
dx_Sc_mat_v = dx_n_Sc_struct.dx_Sc_mat_v_x0;
dy_Sc_mat_v = dx_n_Sc_struct.dy_Sc_mat_v_x0;


dx_r = mntimes2_fast( dx_Sc_mat_v, q_vec);
dy_r = mntimes2_fast( dy_Sc_mat_v, q_vec);


n_vec_norm_n_mat = n_vec_norm_n(:,ones(1,N_q),:,:);



%%[*] dq_(n/||n||)
dq_n = cross_fast( dx_Sc_mat_v_o, dy_r(:,ones(1,N_q),:,:)) + cross_fast( dx_r(:,ones(1,N_q),:,:), dy_Sc_mat_v_o);
nT_dq_n = sum( n_vec_norm_n_mat.*dq_n, 1);

dq_n_per_norm_n = ( dq_n - n_vec_norm_n_mat.*nT_dq_n(ones(1,3),:,:,:) )./norm_n(ones(1,3),ones(1,N_q),:,:);




%%[*] dx_r/||dx_r||
dx_r_mat = permute( dx_r, [ 2 1 3 4]);
dx_r_mat = dx_r_mat(ones(1,N_q),:,:,:);
norm_dx_r = norm2( dx_r);
norm_dx_r_mat = norm_dx_r(ones(1,N_q),ones(1,3),:,:);
dx_r_mat_per_norm_dx_r = dx_r_mat./norm_dx_r_mat;


%%[*] dq_(n/||n||)^T*dx_r/||dx_r||

dq_n_per_norm_nT_dx_r_mat_per_norm_dx_r = sum( permute( dq_n_per_norm_n, [ 2 1 3 4]).*dx_r_mat_per_norm_dx_r, 2);

%%[*] dq_theta_0y
ex = [ 1 0 0].';
ez = [ 0 0 1].';
ex_ezT_m_ez_exT = ex*ez.' - ez*ex.';
dx_rT_ex_ezT_m_ez_exT = sum( dx_r(:,ones(1,3),:,:).*ex_ezT_m_ez_exT(:,:,ones(1,size( dx_r, 3)),ones(1,size( dx_r, 4))), 1);
dx_rT_ex_ezT_m_ez_exT_dx_Sc_mat = repmat( permute( dx_rT_ex_ezT_m_ez_exT, [ 2 1 3 4]), [ 1 N_q/3 1 1]).*dx_Sc_mat_v;
dx_rT_ex_ezT_m_ez_exT_dx_Sc_mat = reshape( dx_rT_ex_ezT_m_ez_exT_dx_Sc_mat, 1, [], size( dx_r, 3), size( dx_r, 4));

dxr_ex_2_p_dxr_ez_2 = repmat( dx_r(1,:,:,:).^2 + dx_r(3,:,:,:).^2, [ 1 N_q 1 1]);

dq_theta_0y = dx_rT_ex_ezT_m_ez_exT_dx_Sc_mat./dxr_ex_2_p_dxr_ez_2;


out = -dq_n_per_norm_nT_dx_r_mat_per_norm_dx_r(:,ones(1,N_q),:,:).*dq_theta_0y(ones(1,N_q),:,:,:);
out2 = dq_theta_0y;

end

