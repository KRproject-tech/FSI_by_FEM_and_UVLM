%%[*] 流体力ベクトル：[p] = Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2) 

%% [2-0] 流体力ベクトル (Mf2_2)
%%[*] 流体力の線形補間のためのノード値
dp_nvec_i = zeros(N_element,3,3);        
dp_nvec_i(:,:,2) = ( (dp_lift1 + Mf2_vec1)*ones(1,3) ).*n_vec_i;
dp_nvec_i(Ny+1:end,:,1) = dp_nvec_i(1:end-Ny,:,2);
dp_nvec_i(1:end-Ny,:,3) = dp_nvec_i(Ny+1:end,:,2);
        

Qf_p_vec_i = zeros(N_q,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
            
    int_StF = zeros(N_q,1);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre求積
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] 流体力の線形補間
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 3 2 1]);
        dp_ni_interp = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i(ii,:,:), 3);
        
        
        i_eta_a = 1;
        for eta_a = p_vec           
            
            
            %%[*] 法線ベクトルと同じ向きに圧力が作用するので，符号は正．
            StF = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_ni_interp.';
            int_StF = int_StF + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StF;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_vec_i(:,ii) = int_StF;            %% ii要素の外力行列 (単位面積当たり)
end



%% [2-1] 流体力行列 ( -(τx*dxΓ + τy*dyΓ)*dt_rc )
%%[*] 流体力の線形補間のためのノード値
dp_lift2_nvec_i = zeros(3,3,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
dp_lift2 = permute( dp_lift2, [3 2 1]);
niT_dp_lift2 = mntimes2( n_vec_mat, dp_lift2);

dp_lift2_nvec_i(:,:,:,2) = niT_dp_lift2;
dp_lift2_nvec_i(:,:,Ny+1:end,1) = dp_lift2_nvec_i(:,:,1:end-Ny,2);
dp_lift2_nvec_i(:,:,1:end-Ny,3) = dp_lift2_nvec_i(:,:,Ny+1:end,2);

Qf_p_lift2_mat_i = zeros(N_q,3*N_element,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
    
   
            
    int_StM = zeros(N_q,3*N_element);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre求積
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] 流体力の線形補間
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_lift2_nvec_i_interp = sum( p_interp_vec(ones(1,3),ones(1,3),:,:).*dp_lift2_nvec_i(:,:,ii,:), 4);
        
        i_eta_a = 1;
        dp_lift2_nvec_i_interp_v = zeros(3,3*N_element);
        dp_lift2_nvec_i_interp_v(:,3*ii-2:3*ii) = dp_lift2_nvec_i_interp;
        for eta_a = p_vec
            
            %%[*] 法線ベクトルと同じ向きに圧力が作用するので，符号は正．            
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_lift2_nvec_i_interp_v;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_lift2_mat_i(:,:,ii) = int_StM;            %% ii要素の外力行列 (単位面積当たり)
end


%% [2-2] 流体力行列 (Mf2_1)
%%[*] 流体力の線形補間のためのノード値
dp_add_nvec0_i = zeros(3,N_element,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
Mf2_mat_i = permute( Mf2_mat, [3 2 1]);
niT_Mf2_mat = mntimes2( n_vec_mat, Mf2_mat_i);

dp_add_nvec0_i(:,:,:,2) = niT_Mf2_mat;
dp_add_nvec0_i(:,:,Ny+1:end,1) = dp_add_nvec0_i(:,:,1:end-Ny,2);
dp_add_nvec0_i(:,:,1:end-Ny,3) = dp_add_nvec0_i(:,:,Ny+1:end,2);

Qf_p_mat0_i = zeros(N_q,N_element,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
    
   
            
    int_StM = zeros(N_q,N_element);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre求積
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] 流体力の線形補間
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_add_ni0_interp = sum( p_interp_vec(ones(1,3),ones(1,N_element),:,:).*dp_add_nvec0_i(:,:,ii,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] 法線ベクトルと同じ向きに圧力が作用するので，符号は正．
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni0_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat0_i(:,:,ii) = int_StM;            %% ii要素の外力行列 (単位面積当たり)
end




%% [2-3] 付加質量行列 (Mf1)
%%[*] 流体力の線形補間のためのノード値
dp_add_nvec_i = zeros(3,N_q_all,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
Mf1_mat_i = permute( Mf1_mat, [3 2 1]);
niT_Mf1_mat = mntimes2( n_vec_mat, Mf1_mat_i);

dp_add_nvec_i(:,:,:,2) = niT_Mf1_mat;
dp_add_nvec_i(:,:,Ny+1:end,1) = dp_add_nvec_i(:,:,1:end-Ny,2);
dp_add_nvec_i(:,:,1:end-Ny,3) = dp_add_nvec_i(:,:,Ny+1:end,2);

Qf_p_mat_i = zeros(N_q,N_q,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
    
    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    
            
    int_StM = zeros(N_q);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre求積
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] 流体力の線形補間
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_add_ni_interp = sum( p_interp_vec(ones(1,3),ones(1,N_q),:,:).*dp_add_nvec_i(:,i_vec,ii,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] 法線ベクトルと同じ向きに圧力が作用するので，符号は正．
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat_i(:,:,ii) = int_StM;            %% ii要素の外力行列 (単位面積当たり)
end

