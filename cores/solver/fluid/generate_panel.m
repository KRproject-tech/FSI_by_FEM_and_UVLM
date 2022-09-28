%% Vortex lattice method の格子生成 


%% [0] パネルノード点座標 [-]( r_i = Sc(x_i,y_i)*q_i, i∈{1,...,N_element} )

Sc_mat_v_panel_1 = zeros(3,N_q,N_element);
Sc_mat_v_panel_4 = Sc_mat_v_panel_1;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
            
    %%[*] 要素座標のノード点 [-]
    x_i = dL*[ 1/4  1/4];   %% 点①，④
    y_i = dW*[ 1    0];     %% 点①，④

    Sc_mat_v_panel_1(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_4(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end

%% [1] 平板末端ノード点座標 [m]( r_i = Sc(x_i,y_i)*q_i, i∈{1,...,N_element} )

Sc_mat_v_panel_end_2 = zeros(3,N_q,N_element);
Sc_mat_v_panel_end_3 = Sc_mat_v_panel_end_2;
%%[*] 要素番号はY方向に増加するので，平板末端の要素番号は; N_element-Ny+1~N_element．
for ii = N_element-Ny+1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% ii要素の長さ [-]
    dW = dW_vec(ii);                        %% ii要素の幅 [-]
            
    %%[*] 要素座標のノード点 [-]
    x_i = dL*[ 1  1];   %% 点②，③
    y_i = dW*[ 1  0];   %% 点②，③

    Sc_mat_v_panel_end_2(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_end_3(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end


%% [2] グローバル行列組立

Sc_mat_panel_global_1 = sparse(3*N_element,N_q_all);
Sc_mat_panel_global_4 = Sc_mat_panel_global_1;
for ii = 1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_global_1(j_vec,i_vec) = Sc_mat_panel_global_1(j_vec,i_vec) + squeeze( Sc_mat_v_panel_1(:,:,ii));
    Sc_mat_panel_global_4(j_vec,i_vec) = Sc_mat_panel_global_4(j_vec,i_vec) + squeeze( Sc_mat_v_panel_4(:,:,ii));
end



%%[*] 要素番号はY方向に増加するので，平板末端の要素番号は; N_element-Ny+1~N_element．
Sc_mat_panel_end_global_2 = sparse(3*N_element,N_q_all);
Sc_mat_panel_end_global_3 = Sc_mat_panel_end_global_2;
for ii = N_element-Ny+1:N_element

    %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_end_global_2(j_vec,i_vec) = Sc_mat_panel_end_global_2(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_2(:,:,ii));
    Sc_mat_panel_end_global_3(j_vec,i_vec) = Sc_mat_panel_end_global_3(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_3(:,:,ii));
end


%%[2-0] 点② [m]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
Sc_mat_panel_global_2 = [ Sc_mat_panel_global_1(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_2(end-3*Ny+1:end,:) - Sc_mat_panel_global_1(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_1(end-3*Ny+1:end,:)];

%%[2-1] 点③ [m] (平板末端のパネルノード点は，平板末端ノード座標より補間)
Sc_mat_panel_global_3 = [ Sc_mat_panel_global_4(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_3(end-3*Ny+1:end,:) - Sc_mat_panel_global_4(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_4(end-3*Ny+1:end,:)];

                      
%% [3] コロケーション点計算 [-]
% Sc_mat_col_global = ( Sc_mat_panel_global_1 + Sc_mat_panel_global_2 + Sc_mat_panel_global_3 + Sc_mat_panel_global_4 )/4;

%%[*] 不均等要素長に対応
Sc_mat_col_global = sparse(3*N_element,N_q_all);
for ii = 1:N_element    
     
    j_vec = 3*ii-2:3*ii;                                %% [ rx ry rz]^T
    
    if ii <= N_element-Ny
        dL_i = dL_vec(ii);                              %% ii要素の長さ [-]
        dL_ip1 = dL_vec(ii+Ny);                         %% ii要素の後方の要素の長さ [-]
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + diag( [ (dL_i + dL_ip1)/(3*dL_i + dL_ip1)     1/2     (dL_i + dL_ip1)/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/2 ...
                                           	+ diag( [ 2*dL_i/(3*dL_i + dL_ip1)              1/2     2*dL_i/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) )/2;
    else
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + ( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/4;
    end
end


%% [4] 法線ベクトル算出 [-]

Sc_mat_31 = Sc_mat_panel_global_3 - Sc_mat_panel_global_1;
Sc_mat_24 = Sc_mat_panel_global_2 - Sc_mat_panel_global_4;



