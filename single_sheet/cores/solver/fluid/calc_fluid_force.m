%% 流体力算出
%%%
%%% �冪/ρf = (Vin + Vw + Vb - dt_r)*(τx*dx_Γ + τy*dy_Γ) + dt_Γ
%%%


%% [0] 定常流体力
%%%
%%% [p]_lift := (Vin + Vw + Vb - dt_r)*(τx*dx_Γ + τy*dy_Γ)
%%% Γ := (A+B)^-1*[ (dt_rc - Vin - Σ_(Ny+1)^Nwake Γ_wake*Σq_wake^T)*n_i ]

%%[0-0] 単位接線ベクトル [-]
r21_vec = r_panel_vec_2 - r_panel_vec_1;
r34_vec = r_panel_vec_3 - r_panel_vec_4;
r14_vec = r_panel_vec_1 - r_panel_vec_4;
r23_vec = r_panel_vec_2 - r_panel_vec_3;

tau_x = (r21_vec + r34_vec)/2;
tau_y = (r14_vec + r23_vec)/2;


d_x_vec = norm_mat( tau_x);
d_y_vec = norm_mat( tau_y);

tau_x = tau_x./d_x_vec;
tau_y = tau_y./d_y_vec;

%%[0-1] 循環勾配の算出: dx_Γ, dy_Γ
d_x_mat = reshape( d_x_vec(:,1), Ny, []).';
d_y_mat = reshape( d_y_vec(:,1), Ny, []).';

%%[*] 差分における先頭要素は，dx_Γ1 = ( (Γ_i - Γ_{i-1})/�凅 )|_{i=1} = Γ_1/�凅とする．
%%% M. Ghommem, Modeling and Analysis for Optimization of Unsteady Aeroelastic
%%% Systems, Doctoral dissertation of Virginia Polytechnic Institute and
%%% State University, p. 138, 2011. 
Gamma_mat = reshape( Gamma, Ny, []).';                                      
dx_Gamma = [ Gamma_mat(1,:);
             diff( Gamma_mat, [], 1)]./d_x_mat;

Gamma_mat2 = [ zeros(Nx,1)  Gamma_mat   zeros(Nx,1)];                       %% y方向は中心差分 (圧力分布の対称性のため)
dy_Gamma = (Gamma_mat2(:,3:end) - Gamma_mat2(:,1:end-2))./(2*d_y_mat);      %% y方向は中心差分
dy_Gamma(:,1) = Gamma_mat(:,1)./d_y_mat(:,1);                               %% 両端は片側差分
dy_Gamma(:,end) = -Gamma_mat(:,end)./d_y_mat(:,end);                        %% 両端は片側差分



%%[0-2] 無次元化流体力算出 [-]
tau_x_dx_Gamma = tau_x.*( reshape( dx_Gamma.', [], 1)*ones(1,3) );
tau_y_dy_Gamma = tau_y.*( reshape( dy_Gamma.', [], 1)*ones(1,3) );


dp_add = (Gamma - old_Gamma)/d_t_wake;                                      %% 付加質量効果
dp_lift = sum( V_surf.*(tau_x_dx_Gamma + tau_y_dy_Gamma), 2);               %% 定常流体力効果
dp_lift1 = sum( V_surf1.*(tau_x_dx_Gamma + tau_y_dy_Gamma), 2);             %% 定常流体力効果
dp_lift2 = -(tau_x_dx_Gamma + tau_y_dy_Gamma);                              %% -(τx*dxΓ + τy*dyΓ)*dt_rc 


dp_vec = dp_lift + dp_add; 
dp_mat = reshape( dp_vec, Ny, []).';

h_dp_add(:,i_wake_time) = dp_add;                                           %% 付加質量効果
h_dp_lift(:,i_wake_time) = dp_lift;                                         %% 定常流体力効果
h_dp_vec(:,i_wake_time) = dp_vec;



%% 付加質量マトリックスの計算




%%[1-1] 付加質量マトリックスの計算(Mf2) [-]
%%% dtΓ := (A+B)^-1*[ n_i^T*Sc_i ]dt^2_q 
%%%         + (A+B)^-1*( [ (-Σ_(i=1)^Nwake Γ_wake*Σdt_q_wake^T)*n_i] 
%%%                         + [(dt_rc - V_wake - Vin)^T*n_i]
%%%                         - dt_A*Γ)
%%%


if ~exist( 'old_dt_q_vec_wake', 'var')

   old_dt_q_vec_wake = dt_q_vec;
end

%%[1-1-0] Σ{Γwake*(dt_q_mat)^T*n}
dt_q1234_wake_mat = dt_generate_q1234_mat( rc_vec, r_wake_1, r_wake_2, r_wake_3, r_wake_4, ...
                                           dt_rc_vec, dt_r_wake_1, dt_r_wake_2, dt_r_wake_3, dt_r_wake_4);
 
Gamma_wake_dt_q1234 = dt_q1234_wake_mat.*( ones(N_element,1)*kron( Gamma_wake.', ones(1,3)) );    
Gamma_wake_dt_q1234 = [ sum( Gamma_wake_dt_q1234(:,1:3:end), 2) sum( Gamma_wake_dt_q1234(:,2:3:end), 2) sum( Gamma_wake_dt_q1234(:,3:3:end), 2)];

Gamma_wake_dt_q1234_n = sum( Gamma_wake_dt_q1234.*n_vec_i, 2);

                                 





%%[1-1-2] dt_A = (dt_q_mat)^T*n + q_mat^T*dt_n
dt_q1234_mat = dt_generate_q1234_mat( rc_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, ...
                                      dt_rc_vec, dt_r_panel_vec_1, dt_r_panel_vec_2, dt_r_panel_vec_3, dt_r_panel_vec_4);
                                  
dt_q_mat_ni = inner_mat( dt_q1234_mat, n_vec_i_mat);
dt_q_mat_ni = dt_q_mat_ni(:,1:3:end);                                      	%% inner_mat関数における3成分だけのコピーは不要．    

dt_n_vec_i_mat = repmat( dt_n_vec_i, [ 1 N_element]);
q_mat_dt_ni = inner_mat( q1234_mat, dt_n_vec_i_mat);
q_mat_dt_ni = q_mat_dt_ni(:,1:3:end);                                      	%% inner_mat関数における3成分だけのコピーは不要．                        

dt_Amat = dt_q_mat_ni + q_mat_dt_ni;
h_dt_Amat(:,:,i_wake_time) = dt_Amat;
h_Amat(:,:,i_wake_time) = A_mat;


dt_Amat1 = dt_q_mat_ni;


q_mat_Gamma_vec = q1234_mat.*( ones(N_element,1)*kron( Gamma.', ones(1,3)) );
q_mat_Gamma_vec = [ sum( q_mat_Gamma_vec(:,1:3:end), 2)  sum( q_mat_Gamma_vec(:,2:3:end), 2) sum( q_mat_Gamma_vec(:,3:3:end), 2)];

dt_Amat2_Gamma = q_mat_Gamma_vec;


%%[1-1-3] M_f2
%%[*] 流体力ベクトル：[p] = Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2) 

% Mf2_vec = (A_mat + B_mat)\(  -Gamma_wake_dt_q1234_n ...
%                                 + sum( (dt_rc_vec - V_in - V_wake_plate).*dt_n_vec_i, 2) ...
%                                 - dt_Amat*Gamma );

% Mf2_vec = (A_mat + B_mat)\(  -Gamma_wake_dt_q1234_n - dt_Amat*Gamma );
Mf2_vec1 = A_mat\(  -Gamma_wake_dt_q1234_n  );                                  %% dt_A = [Σ{ (dt_q_wake)^T*ni + q_wake^T*dt_ni }]は構造モデルに組み込む．


                            
Mf2_mat = inv(A_mat);                            


%%[1-2] 付加質量マトリックスの計算(Mf1) [-]

nvec_Sc_global = zeros(N_element,N_q_all);
for ii = 1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    nvec_Sc_global(ii,i_vec) = nvec_Sc_global(ii,i_vec) + n_vec_i(ii,:)*Sc_mat_col_global(j_vec,i_vec);
end
Mf1_mat = A_mat\nvec_Sc_global;




%% added mass effect

% h_dp_add_estimate(:,i_wake_time) = Mf1_mat*(dt_q_vec - old_dt_q_vec_wake)/d_t_wake ...
%                                     + Mf2_vec;                                %% 付加質量効果
h_dp_add_estimate(:,i_wake_time) = Mf1_mat*(dt_q_vec - old_dt_q_vec_wake)/d_t_wake ...
                                    ...
                                    + Mf2_mat*( sum( (dt_rc_vec - V_in - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i, 2) - dt_Amat1*Gamma)...
                                    + Mf2_vec1;                                 %% 付加質量効果

                                
old_dt_q_vec_wake = dt_q_vec;
old_Amat = A_mat;




%% 1step前の値を更新

% old_Qf_p_global = Qf_p_global;
% old_Qf_p_mat_global = Qf_p_mat_global;
% old_Qf_p_mat0_global = Qf_p_mat0_global;
% old_Qf_p_lift2_mat_global = Qf_p_lift2_mat_global;

%%[*] Kuttaの条件を満たさせるため，1step前の循環値を用いる．
old_Gamma = Gamma;  






%% 流体力行列組み立て

if flag_fluid_bench
    
    Qf_p_global = zeros(N_q_all,1);
    Qf_p_mat_global = sparse(N_q_all,N_q_all);
    Qf_p_mat0_global = zeros(N_q_all,N_element);
    Qf_p_lift2_mat_global = zeros(N_q_all,3*N_element);
    
else
    
    if coupling_flag == 1

        %%[*] 流体力ベクトル：[p] = p_lift + Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2) 
        calc_fluid_force_strong;
    else 

        %%[*] 流体力ベクトル：[p] = p_lift + p_add 
        calc_fluid_force_weak;
    end



    %% [3] グローバル行列組立

    Qf_p_global = zeros(N_q_all,1);
    for ii = 1:N_element

        %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
        %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
        i_vec = i_vec_v{ii};

        Qf_p_global(i_vec,1) = Qf_p_global(i_vec,1) + Qf_p_vec_i(:,ii);
    end

    Qf_p_mat_global = sparse(N_q_all,N_q_all);
    Qf_p_mat0_global = zeros(N_q_all,N_element);
    Qf_p_lift2_mat_global = zeros(N_q_all,3*N_element);
    for ii = 1:N_element

        %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
        %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
        i_vec = i_vec_v{ii};

        Qf_p_mat0_global(i_vec,:) = Qf_p_mat0_global(i_vec,:) + squeeze( Qf_p_mat0_i(:,:,ii));
        Qf_p_lift2_mat_global(i_vec,:) = Qf_p_lift2_mat_global(i_vec,:) + squeeze( Qf_p_lift2_mat_i(:,:,ii));
        Qf_p_mat_global(i_vec,i_vec) = Qf_p_mat_global(i_vec,i_vec) + squeeze( Qf_p_mat_i(:,:,ii));
    end
                                                         
end