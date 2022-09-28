

%% [0] 体積力 (Qf = ∫_A S(x,y)^t*F dA, F ∈ R^3)

Qf_time_vec_i = zeros(N_q,N_element);
for ii = 1:N_element

    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
            
    int_StF = zeros(N_q,1);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre求積
        i_eta_a = 1;
        for eta_a = p_vec
            
            StF = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*q_in_vec;
            int_StF = int_StF + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StF;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_time_vec_i(:,ii) = int_StF;         %% ii要素の外力行列
end


%% [2] グローバル行列組立

Qf_time_global = zeros(N_q_all,1);
for ii = 1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    Qf_time_global(i_vec,1) = Qf_time_global(i_vec,1) + Qf_time_vec_i(:,ii);
end

