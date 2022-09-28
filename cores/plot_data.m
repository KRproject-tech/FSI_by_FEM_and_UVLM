clc
clear all
close all hidden

%% delete
delete( '*.asv')
delete( '*.log')

%% path
add_pathes

%% initializing
initializing


%% parameter
param_setting
tmp.movie_format = movie_format;
tmp.Snapshot_tmin = Snapshot_tmin;
tmp.Snapshot_tmax = Snapshot_tmax;
tmp.i_snapshot = i_snapshot;
tmp.panel_node_plot = panel_node_plot;
tmp.pressure_interp_plot = pressure_interp_plot;


%% version check
version_check


%% multi threading
maxNumCompThreads( core_num);

%% load
disp( 'Now loading ....')
load ./save/NUM_DATA
movie_format = tmp.movie_format;
Snapshot_tmin = tmp.Snapshot_tmin;
Snapshot_tmax = tmp.Snapshot_tmax;
i_snapshot = tmp.i_snapshot;
panel_node_plot = tmp.panel_node_plot;
pressure_interp_plot = tmp.pressure_interp_plot;



%% plot

i_ax = 1;


%% [0] Plotting the Finite Element Mesh
X = zeros(4,N_element) ;
Y = zeros(4,N_element) ;

for ii = 1:N_element
    X(:,ii) = coordinates(nodes(ii,:),1) ;
    Y(:,ii) = coordinates(nodes(ii,:),2) ;
end



h_fig(1) = figure(1);
set( h_fig(1), 'Position', [100 700 600 300])
h_ax(i_ax) = axes( 'Parent', h_fig(1), 'FontSize', 15);

patch( X, Y, 'w', 'Parent', h_ax(i_ax))
axis( h_ax(i_ax), 'equal')
xlabel( h_ax(i_ax), '{\itx} position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
ylabel( h_ax(i_ax), '{\ity} position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)

set( h_ax(i_ax), 'FontName', 'Times New Roman')

i_ax = i_ax + 1;





%% [1] Plotting the nodes

%%[1-0] plot behavior
h_fig(2) = figure(2);
set( h_fig(2), 'Position', [100 100 1200 500], 'render', 'zbuffer', 'doublebuffer','on')
h_ax(i_ax) = axes( 'Parent', h_fig(2), 'FontSize', 15);



h_txt(1) = text( -Length/2, Width, 0, [ 'Time = ', num2str( 0, '%0.3f'), ' [s]'],...
                'FontSize', 12, 'FontName', 'Times New Roman', 'BackgroundColor', 'g');

h_plot(1) = patch( 0, 0, 0, 'r', 'Parent', h_ax(i_ax));
hold( h_ax(i_ax), 'on')
if panel_node_plot == 1;
    h_plot(2) = plot3( h_ax(i_ax), 0, 0, 0, '.');                          	%% パネルノード点 [m]
    h_plot(3) = plot3( h_ax(i_ax),  0, 0, 0, '*');                       	%% コロケーション点 [m]
    h_plot(4) = quiver3( h_ax(i_ax),  0, 0, 0, 0, 0, 0, 'r');             	%% 流体力ベクトル [Pa]
end
h_plot(5) = quiver3( h_ax(i_ax),  0, 0, 0, 0, 0, 0, 'g');               	%% 表面流速ベクトル [m/s]
h_plot(6) = patch( 0, 0, 0, 0, 'Parent', h_ax(i_ax), 'LineStyle', 'none');	%% Wake (1st sheet)
h_plot(16) = patch( 0, 0, 0, 0, 'Parent', h_ax(i_ax), 'LineStyle', 'none');	%% Wake (2nd sheet)
if pressure_interp_plot == 1;
    h_plot(7) = quiver3( h_ax(i_ax),  0, 0, 0, 0, 0, 0, 'b');             	%% 補間流体力ベクトル [Pa]
end
all_Gamma_wake = cell2mat( reshape( [ h_Gamma_wake; h_Gamma_wake_1], [], 1));
caxis( h_ax(i_ax), [ min( all_Gamma_wake) max( all_Gamma_wake)])
h_col = colorbar;
ylabel( h_col, ' Circulation on wake  \Gamma [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
light
lighting gouraud 
view( h_ax(i_ax), [1 -2 1])
axis( h_ax(i_ax), 'equal')
grid( h_ax(i_ax), 'on')
xlim( h_ax(i_ax), [-Length 5*Length])
ylim( h_ax(i_ax), [-Width 2*Width])
zlim( h_ax(i_ax), [-1.5*Length 1.5*Length])
xlabel( h_ax(i_ax), '{\itX}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
ylabel( h_ax(i_ax), '{\itY}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
zlabel( h_ax(i_ax), '{\itZ}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
set( h_ax(i_ax), 'FontName', 'Times New Roman')

i_ax = i_ax + 1;





%%[1-1] Snapshots
h_fig(3) = figure(3);
set( h_fig(3), 'Position', [700 600 600 400])
h_ax(i_ax) = axes( 'Parent', h_fig(3), 'FontSize', 15);

light
lighting gouraud 
view( h_ax(i_ax), [1 2 1])
axis( h_ax(i_ax), 'equal')
grid( h_ax(i_ax), 'on')
xlim( h_ax(i_ax), [-0.5*Length 1.5*Length])
ylim( h_ax(i_ax), [-1.5*Width 2*Width])
zlim( h_ax(i_ax), [-0.5*Length 0.5*Length])
xlabel( h_ax(i_ax), '{\itX}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
ylabel( h_ax(i_ax), '{\itY}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
zlabel( h_ax(i_ax), '{\itZ}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
set( h_ax(i_ax), 'FontName', 'Times New Roman')

i_ax = i_ax + 1;









movie_data = struct( 'cdata', [], 'colormap', []);
i_time = 1;
i_wake_time = 1;
i_movie = 1;
idx_r = reshape( [1:N_qi:N_q_all; 2:N_qi:N_q_all; 3:N_qi:N_q_all], 1, []);  
idx_r1 = N_q_all+idx_r;
time_m = 0:d_t:time;


%%[1-2] 補間圧力ベクトルの時刻歴の算出
if pressure_interp_plot == 1;
    
    lgth_time = length( time_m(time_m <= Snapshot_tmax))*ceil(d_t/dt_wake_per_dt);
    
    h_r_vec_interp = zeros(3*length( p_vec)*length( p_vec)*2*N_element,lgth_time);
    h_dp_ni_interp = zeros(3*length( p_vec)*length( p_vec)*2*N_element,lgth_time);
    
    for time = time_m(time_m <= Snapshot_tmax)

        if mod( i_time, dt_wake_per_dt) == 0
            
            dp_vec_all = h_dp_vec(:,i_wake_time);
            dp_vec = dp_vec_all(1:end/2);
            dp_vec_1 = dp_vec_all(end/2+1:end);

            
            q_vec = h_X_vec(1:N_q_all,i_time);
            q_vec_1 = h_X_vec(N_q_all+1:2*N_q_all,i_time);
            
            %%[*-0] (シート1)   
            r_vec_interp = zeros(3,length( p_vec),length( p_vec),N_element);
            dp_ni_interp = zeros(3,length( p_vec),length( p_vec),N_element);
            %%[*-1] (シート2)   
            r_vec_interp_1 = zeros(3,length( p_vec),length( p_vec),N_element);
            dp_ni_interp_1 = zeros(3,length( p_vec),length( p_vec),N_element);

            %%[*] 流体力の線形補間のためのノード値
            n_vec_i = ( dp_vec*ones(1,3) ).*h_n_vec(:,:,i_wake_time);
            n_vec_i_1 = ( dp_vec_1*ones(1,3) ).*h_n_vec_1(:,:,i_wake_time);
            %%[*-0] (シート1)        
            dp_nvec_i = zeros(N_element,3,3);        
            dp_nvec_i(:,:,2) = n_vec_i;
            dp_nvec_i(Ny+1:end,:,1) = dp_nvec_i(1:end-Ny,:,2);
            dp_nvec_i(1:end-Ny,:,3) = dp_nvec_i(Ny+1:end,:,2);
            %%[*-1] (シート2)
            dp_nvec_i_1 = zeros(N_element,3,3);        
            dp_nvec_i_1(:,:,2) = n_vec_i_1;
            dp_nvec_i_1(Ny+1:end,:,1) = dp_nvec_i_1(1:end-Ny,:,2);
            dp_nvec_i_1(1:end-Ny,:,3) = dp_nvec_i_1(Ny+1:end,:,2);

            for ii = 1:N_element
                dL = dL_vec(ii);                        %% ii要素の長さ [-]
                dW = dW_vec(ii);                        %% ii要素の幅 [-]

                %% 1ノード当たり9成分 ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
                %% 1要素当たり36成分　( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
                i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
                i_vec = reshape(i_vec.',1,[]);
                q_i_vec = q_vec(i_vec);
                q_i_vec_1 = q_vec_1(i_vec);


                i_xi_a = 1;
                for xi_a = p_vec                        %% Gauss-Legendre求積

                    x_i = dL*(xi_a + 1)/2;
                    p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
                    p_interp_vec = permute( p_interp_vec, [ 3 2 1]);
                    dp_ni = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i(ii,:,:), 3);
                    dp_ni_1 = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i_1(ii,:,:), 3);

                    i_eta_a = 1;
                    for eta_a = p_vec


                        y_i = dW*(eta_a + 1)/2;
                        %%[*-0] (シート1)
                        r_vec_interp(:,i_xi_a,i_eta_a,ii) = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii)*q_i_vec;
                        dp_ni_interp(:,i_xi_a,i_eta_a,ii) = dp_ni;
                        %%[*-1] (シート2)
                        r_vec_interp_1(:,i_xi_a,i_eta_a,ii) = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii)*q_i_vec_1;
                        dp_ni_interp_1(:,i_xi_a,i_eta_a,ii) = dp_ni_1;

                        i_eta_a = i_eta_a+1;
                    end
                    i_xi_a = i_xi_a+1;
                end
            end
            
            h_r_vec_interp(:,i_wake_time) = [ r_vec_interp(:);
                                              r_vec_interp_1(:)];
            h_dp_ni_interp(:,i_wake_time) = [ dp_ni_interp(:);
                                              dp_ni_interp_1(:)];
            
            i_wake_time = i_wake_time + 1;            
        end

        i_time = i_time + 1;
    end
end




%%[1-3] アニメーション

i_time = 1;
i_wake_time = 1;
for time = time_m(time_m <= Snapshot_tmax)
    
    if mod( i_time, i_snapshot) == 0
        
        disp( [ 'Time = ', num2str( time, '%0.4f'), ' [s]'])

        %%[1-0] ノード変位データ抽出
        r_vec = reshape( h_X_vec(idx_r,i_time), 3, []);
        r_vec_1 = reshape( h_X_vec(idx_r1,i_time), 3, []);

        %%[1-1] 1要素当たりのノードの座標取得
        X = zeros(4,N_element);
        Y = zeros(4,N_element);
        Z = zeros(4,N_element);
        for ii = 1:N_element
            X(:,ii) = r_vec(1,nodes(ii,:));
            Y(:,ii) = r_vec(2,nodes(ii,:));
            Z(:,ii) = r_vec(3,nodes(ii,:));
        end
        
        X_1 = X;
        Y_1 = Y;
        Z_1 = Z;
        for ii = 1:N_element
            X_1(:,ii) = r_vec_1(1,nodes(ii,:));
            Y_1(:,ii) = r_vec_1(2,nodes(ii,:));
            Z_1(:,ii) = r_vec_1(3,nodes(ii,:));
        end

        %%[1-2] plot更新
        set( h_plot(1), 'XData', [ X X_1], 'YData', [ Y Y_1], 'ZData', [ Z Z_1]);
                
        set( h_txt(1), 'String', [ 'Time = ', num2str( time, '%0.3f'), ' [-]']);
        
        %%[1-3] 挙動のスナップショット
        if Snapshot_tmin <= time && time <= Snapshot_tmax
            
            patch( [ X X_1], [ Y Y_1], [ Z Z_1], 'r', 'Parent', h_ax(3), 'FaceAlpha', 0.5);
            Xv = X(:);
            Yv = Y(:);
            Zv = Z(:);
            Rv = sqrt( Xv.^2+ Yv.^2 + Zv.^2);
            [ dummy, idx_max_r] = max( Rv);
%             text( Xv(idx_max_r(1)), Yv(idx_max_r(1)), Zv(idx_max_r(1)), [ 'Time = ', num2str( time, '%0.3f'), ' [-]'],...
%                     'FontSize', 12, 'FontName', 'Times New Roman', 'BackgroundColor', 'g', 'Parent', h_ax(3));
        end
        
        movie_data(i_movie) = getframe( h_fig(2));
        movie_time(i_movie) = time; %#ok<AGROW>
        
        i_movie = i_movie+1;
    end
    
    
    

        
    if mod( i_time, dt_wake_per_dt) == 0 && i_wake_time <= length( h_r_wake)

        %%[1-2-0-0] パネルノード点 [-] (翼表面)
        x_node = h_r_panel_vec(:,1,i_wake_time);
        y_node = h_r_panel_vec(:,2,i_wake_time);
        z_node = h_r_panel_vec(:,3,i_wake_time);       
        %%[1-2-0-1] パネルノード点 [-] (翼裏面)
        x_node_1 = h_r_panel_vec_1(:,1,i_wake_time);
        y_node_1 = h_r_panel_vec_1(:,2,i_wake_time);
        z_node_1 = h_r_panel_vec_1(:,3,i_wake_time); 
        
        %%[1-2-1-0] コロケーション点 [-] (翼表面)
        r_col = h_rcol_vec(:,:,i_wake_time);
        %%[1-2-1-1] コロケーション点 [-] (翼裏面)
        r_col_1 = h_rcol_vec_1(:,:,i_wake_time);
        
        %%[1-2-2-0] 流体力ベクトル [Pa]: f = -dp*n*dS (翼表面)
        dp_vec_all = h_dp_vec(:,i_wake_time);
        dp_vec = dp_vec_all(1:end/2);
        dp_vec_1 = dp_vec_all(end/2+1:end);
        n_vec_i = ( dp_vec*ones(1,3) ).*h_n_vec(:,:,i_wake_time);
        %%[1-2-2-1] 流体力ベクトル [Pa]: f = -dp*n*dS (翼裏面)
        n_vec_i_1 = ( dp_vec_1*ones(1,3) ).*h_n_vec_1(:,:,i_wake_time);


        if panel_node_plot == 1;
            set( h_plot(2), 'XData', [ x_node(:); x_node_1(:)], 'YData', [ y_node(:); y_node_1(:)], 'ZData', [ z_node(:); z_node_1(:)]);
            set( h_plot(3), 'XData', [ r_col(:,1); r_col_1(:,1)], 'YData', [ r_col(:,2); r_col_1(:,2)], 'ZData', [ r_col(:,3); r_col_1(:,3)]);
            set( h_plot(4), 'XData', r_col(:,1), 'YData', r_col(:,2), 'ZData', r_col(:,3), 'UData', n_vec_i(:,1), 'VData', n_vec_i(:,2), 'WData', n_vec_i(:,3));
        end
        
        %% 補間流体力ベクトル [Pa]
        if pressure_interp_plot == 1;           
            
            set( h_plot(7), 'XData', h_r_vec_interp(1:3:end,i_wake_time), 'YData', h_r_vec_interp(2:3:end,i_wake_time), 'ZData', h_r_vec_interp(3:3:end,i_wake_time),...
                            'UData', h_dp_ni_interp(1:3:end,i_wake_time), 'VData', h_dp_ni_interp(2:3:end,i_wake_time), 'WData', h_dp_ni_interp(3:3:end,i_wake_time));             	
        end
        
        

        %%[1-2-3] 表面流速ベクトル [m/s]
        V_surf = h_V_surf(:,:,i_wake_time);                        
        V_wake_end = h_V_wake_end(:,:,i_wake_time);
        V_surf_all = [ V_surf;
                       V_wake_end];

        r_end = h_r_end(:,:,i_wake_time);
        r_point_all = [ r_col;
                        r_col_1;
                        r_end];
        set( h_plot(5), 'XData', r_point_all(:,1), 'YData', r_point_all(:,2), 'ZData', r_point_all(:,3), 'UData', V_surf_all(:,1), 'VData', V_surf_all(:,2), 'WData', V_surf_all(:,3));


        %%[1-2-4] Wake plot
        %%[*-0] シート1
        r_wake = h_r_wake{i_wake_time};
        
        N_wake = size( r_wake, 1)/4;
        ii_wake = 1:N_wake;
        r_wake1 = r_wake(ii_wake,:);
        r_wake2 = r_wake(N_wake+ii_wake,:);
        r_wake3 = r_wake(2*N_wake+ii_wake,:);
        r_wake4 = r_wake(3*N_wake+ii_wake,:);

        R_wake = zeros(4*N_wake,3);
        R_wake(1:4:end,:) = r_wake4;
        R_wake(2:4:end,:) = r_wake3;
        R_wake(3:4:end,:) = r_wake2;
        R_wake(4:4:end,:) = r_wake1;

        Xw = reshape( R_wake(:,1), 4, []);
        Yw = reshape( R_wake(:,2), 4, []);
        Zw = reshape( R_wake(:,3), 4, []);
        Gamma_wake = h_Gamma_wake{i_wake_time};
        
       	%%[*-1] シート2  
        r_wake_1 = h_r_wake_1{i_wake_time};
        
        N_wake_1 = size( r_wake_1, 1)/4;
        ii_wake_1 = 1:N_wake_1;
        r_wake1_1 = r_wake_1(ii_wake_1,:);
        r_wake2_1 = r_wake_1(N_wake_1+ii_wake_1,:);
        r_wake3_1 = r_wake_1(2*N_wake_1+ii_wake_1,:);
        r_wake4_1 = r_wake_1(3*N_wake_1+ii_wake_1,:);

        R_wake_1 = zeros(4*N_wake_1,3);
        R_wake_1(1:4:end,:) = r_wake4_1;
        R_wake_1(2:4:end,:) = r_wake3_1;
        R_wake_1(3:4:end,:) = r_wake2_1;
        R_wake_1(4:4:end,:) = r_wake1_1;

        Xw_1 = reshape( R_wake_1(:,1), 4, []);
        Yw_1 = reshape( R_wake_1(:,2), 4, []);
        Zw_1 = reshape( R_wake_1(:,3), 4, []);
        Gamma_wake_1 = h_Gamma_wake_1{i_wake_time};
        
        set( h_plot(6), 'XData', Xw, 'YData', Yw, 'ZData', Zw, 'CData', Gamma_wake);
        set( h_plot(16), 'XData', Xw_1, 'YData', Yw_1, 'ZData', Zw_1, 'CData', Gamma_wake_1);
        
        

        i_wake_time = i_wake_time + 1;
    end
        
    

    
    drawnow
    
    i_time = i_time + 1;
end

i_wake_time = i_wake_time - 1;


%% Velocity field

h_fig(4) = figure(4);
set( h_fig(4), 'Position', [100 100 1200 500], 'render', 'OpenGL')
h_ax(i_ax) = axes( 'Parent', h_fig(4), 'FontSize', 15);





%%[2-0] ノード変位データ抽出
r_vec = reshape( h_X_vec(idx_r,i_time), 3, []);
r_vec_1 = reshape( h_X_vec(idx_r1,i_time), 3, []);

%%[2-1] 1要素当たりのノードの座標取得
X = zeros(4,N_element) ;
Y = zeros(4,N_element) ;
Z = zeros(4,N_element) ;
for ii = 1:N_element
    X(:,ii) = r_vec(1,nodes(ii,:));
    Y(:,ii) = r_vec(2,nodes(ii,:));
    Z(:,ii) = r_vec(3,nodes(ii,:));
end

X_1 = X;
Y_1 = Y;
Z_1 = Z;
for ii = 1:N_element
X_1(:,ii) = r_vec_1(1,nodes(ii,:));
Y_1(:,ii) = r_vec_1(2,nodes(ii,:));
Z_1(:,ii) = r_vec_1(3,nodes(ii,:));
end


%%[2-2] パネルノード点 [m]
%%[2-2-0] 翼表面
r_node = h_r_panel_vec(:,:,i_wake_time);

N_node = size( r_node, 1)/4;
ii_node = 1:N_node;
r_node1 = r_node(ii_node,:);
r_node2 = r_node(N_node+ii_node,:);
r_node3 = r_node(2*N_node+ii_node,:);
r_node4 = r_node(3*N_node+ii_node,:);
%%[2-2-1] 翼裏面
r_node_1 = h_r_panel_vec_1(:,:,i_wake_time);

N_node_1 = size( r_node_1, 1)/4;
ii_node_1 = 1:N_node_1;
r_node1_1 = r_node_1(ii_node_1,:);
r_node2_1 = r_node_1(N_node_1+ii_node_1,:);
r_node3_1 = r_node_1(2*N_node_1+ii_node_1,:);
r_node4_1 = r_node_1(3*N_node_1+ii_node_1,:);

r_node1_all = [ r_node1; 
                r_node1_1];
r_node2_all = [ r_node2; 
                r_node2_1];
r_node3_all = [ r_node3; 
                r_node3_1];
r_node4_all = [ r_node4; 
                r_node4_1];


%%[2-3] Wake node
%%[2-3-0] シート1
r_wake = h_r_wake{i_wake_time};

N_wake = size( r_wake, 1)/4;
ii_wake = 1:N_wake;
r_wake1 = r_wake(ii_wake,:);
r_wake2 = r_wake(N_wake+ii_wake,:);
r_wake3 = r_wake(2*N_wake+ii_wake,:);
r_wake4 = r_wake(3*N_wake+ii_wake,:);
%%[2-3-1] シート2
r_wake_1 = h_r_wake_1{i_wake_time};

N_wake_1 = size( r_wake_1, 1)/4;
ii_wake_1 = 1:N_wake_1;
r_wake1_1 = r_wake_1(ii_wake_1,:);
r_wake2_1 = r_wake_1(N_wake_1+ii_wake_1,:);
r_wake3_1 = r_wake_1(2*N_wake_1+ii_wake_1,:);
r_wake4_1 = r_wake_1(3*N_wake_1+ii_wake_1,:);
%%[2-3-2] 統合
r_wake1_all = [ r_wake1; 
                r_wake1_1];
r_wake2_all = [ r_wake2; 
                r_wake2_1];
r_wake3_all = [ r_wake3; 
                r_wake3_1];
r_wake4_all = [ r_wake4; 
                r_wake4_1];
            
            
%%[2-4] 対象座標 [m]
[ X_mat, Z_mat] = meshgrid( linspace( -Length, 4*Length, 60),  linspace( -2.0*Length, 2.0*Length, 160));

[ N_row, N_col] = size( X_mat);

X_v = X_mat(:);
Z_v = Z_mat(:);
Y_v = Width/2*ones(N_row*N_col,1);

r_xyz = [ X_v Y_v Z_v];

%%[2-5] 流速分布
Gamma_wake = h_Gamma_wake{i_wake_time};
Gamma_wake_1 = h_Gamma_wake_1{i_wake_time};
Gamma_wake_all = [  Gamma_wake;
                    Gamma_wake_1];


V_wake = V_wake_func( r_xyz, r_wake1_all, r_wake2_all, r_wake3_all, r_wake4_all, Gamma_wake_all, var_param, 0);
V_gamma = V_wake_func( r_xyz, r_node1_all, r_node2_all, r_node3_all, r_node4_all, Gamma_all, var_param, 0);
V_in = ones(N_row*N_col,1)*V_in(1,:);

V_xyz = V_wake + V_gamma + V_in;



%%[2-6] 流線

Z0_pos = linspace( -2.0*Length, 2.0*Length, 50);
X0_pos = (-Length + eps)*ones(length( Z0_pos),1);


tri = delaunay( r_xyz(:,1), r_xyz(:,3));
FlowP=TriStream( tri, r_xyz(:,1), r_xyz(:,3), V_xyz(:,1), V_xyz(:,3), X0_pos, Z0_pos);
h_plot_streamline = ptristream( FlowP, 0.5*Width, 'b');
set( h_plot_streamline, 'LineWidth', 1)



%%[2-7] plot
idx_plot = 1:2:N_row*N_col;                                                 %% 流速ベクトルをplotするときはデータ数を間引く．

dp_vec_all = h_dp_vec(:,i_wake_time);
%%[2-7-0] 翼表面 
dp_vec = dp_vec_all(1:end/2);
%%[2-7-1] 翼裏面
dp_vec_1 = dp_vec_all(end/2+1:end);

patch( X, Y, Z, 0, 'CData', dp_vec, 'Parent', h_ax(i_ax));                  %% 翼表面の圧力分布 [-]
hold( h_ax(i_ax), 'on')
patch( X_1, Y_1, Z_1, 0, 'CData', dp_vec_1, 'Parent', h_ax(i_ax));          %% 翼裏面の圧力分布 [-]  
% quiver3( h_ax(i_ax), r_xyz(:,1), r_xyz(:,2), r_xyz(:,3), V_xyz(:,1), V_xyz(:,2), V_xyz(:,3), 'r', 'AutoScaleFactor', 2);          %% 表面流速ベクトル [m/s]
quiver5( r_xyz(idx_plot,1), r_xyz(idx_plot,2), r_xyz(idx_plot,3), V_xyz(idx_plot,1), V_xyz(idx_plot,2), V_xyz(idx_plot,3), 'r');    %% 表面流速ベクトル [m/s]

view( h_ax(i_ax), [0 -1 0])
axis( h_ax(i_ax), 'equal')
grid( h_ax(i_ax), 'on')
xlim( h_ax(i_ax), [-Length 5*Length])
ylim( h_ax(i_ax), [-Width 2*Width])
zlim( h_ax(i_ax), [ -1.5*Length, 1.5*Length])

h_col2 = colorbar;
xlabel( h_ax(i_ax), '{\itX} position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
ylabel( h_ax(i_ax), '{\itY} position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
zlabel( h_ax(i_ax), '{\itZ} position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
ylabel( h_col2, 'Pressure [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
set( h_ax(i_ax), 'FontName', 'Times New Roman')

i_ax = i_ax + 1;




%% エネルギ収支

h_W_inertia = [ zeros(2,1) (h_E_inertia(:,3:end) - h_E_inertia(:,1:end-2))/(2*d_t) nan(2,1)];
h_W_em = [ zeros(2,1) (h_E_em(:,3:end) - h_E_em(:,1:end-2))/(2*d_t) nan(2,1)];
h_W_ek = [ zeros(2,1) (h_E_ek(:,3:end) - h_E_ek(:,1:end-2))/(2*d_t) nan(2,1)];
h_W_total_m = h_W_inertia + h_W_em + h_W_ek + h_W_dm + h_W_dk;
h_W_total_m(:,length( time_m):end) = nan;


h_fig(5) = figure(5);
set( h_fig(5), 'Position', [100 100 1200 500])
h_ax(i_ax) = axes( 'Parent', h_fig(5), 'FontSize', 15);

plot( h_ax(i_ax), time_m(1:end-1), h_W_total_m(1,1:length( time_m)-1), 'b--', time_m(1:end-1), h_W_total_m(2,1:length( time_m)-1), 'r--', 'LineWidth', 3)
xlabel( h_ax(i_ax), 'Time [-]')
ylabel( h_ax(i_ax), 'Work rate [-]') 
ylim( h_ax(i_ax), [ -0.2 0.2])
grid( h_ax(i_ax), 'on')
hold( h_ax(i_ax), 'on')
plot( h_ax(i_ax), time_wake_m(1:end-1), h_W_f_ext(1,1:length( time_wake_m)-1), 'co-', time_wake_m(1:end-1), h_W_f_ext(2,1:length( time_wake_m)-1), 'mo-')
plot( h_ax(i_ax), time_m(1:end-1), h_W_dk(1,1:length( time_m)-1), 'g-', time_m(1:end-1), h_W_dk(2,1:length( time_m)-1), 'g--')
plot( h_ax(i_ax), time_m(1:end-1), h_W_dm(1,1:length( time_m)-1), 'k-', time_m(1:end-1), h_W_dm(2,1:length( time_m)-1), 'k--')


legend( '{\itW_{total,{\rm1}}}', '{\itW_{total,{\rm2}}}', '{\itW_{f,{\rm1}}}', '{\itW_{f,{\rm2}}}',...
        '{\itW_{dk,{\rm1}}}', '{\itW_{dk,{\rm2}}}', '{\itW_{dm,{\rm1}}}', '{\itW_{dm,{\rm2}}}')




i_ax = i_ax + 1;



%% modes

if exist( 'mode_num', 'var')
    

    for i_mode = 1:mode_num
       
        %%[0] 固有ノード変位データ抽出
        Phi_r_vec = reshape( Phi_q_mat_BC(idx_r,i_mode), 3, []);
        Phi_r_vec_1 = reshape( Phi_q_mat_BC(idx_r1,i_mode), 3, []);
        
        Phi_dr_vec = reshape( Phi_dq_mat_BC(idx_r,i_mode), 3, []);
        Phi_dr_vec_1 = reshape( Phi_dq_mat_BC(idx_r1,i_mode), 3, []);
        

        %%[1] 1要素当たりのノードの座標取得
        X = zeros(4,N_element);
        Y = zeros(4,N_element);
        Z = zeros(4,N_element);
        for ii = 1:N_element
            X(:,ii) = Phi_r_vec(1,nodes(ii,:));
            Y(:,ii) = Phi_r_vec(2,nodes(ii,:));
            Z(:,ii) = Phi_r_vec(3,nodes(ii,:));
        end
        
        X_1 = X;
        Y_1 = Y;
        Z_1 = Z;
        for ii = 1:N_element
            X_1(:,ii) = Phi_r_vec_1(1,nodes(ii,:));
            Y_1(:,ii) = Phi_r_vec_1(2,nodes(ii,:));
            Z_1(:,ii) = Phi_r_vec_1(3,nodes(ii,:));
        end
        
        dX = X;
        dY = Y;
        dZ = Z;
        for ii = 1:N_element
            dX(:,ii) = Phi_dr_vec(1,nodes(ii,:));
            dY(:,ii) = Phi_dr_vec(2,nodes(ii,:));
            dZ(:,ii) = Phi_dr_vec(3,nodes(ii,:));
        end
        
        dX_1 = dX;
        dY_1 = dY;
        dZ_1 = dZ;
        for ii = 1:N_element
            dX_1(:,ii) = Phi_dr_vec_1(1,nodes(ii,:));
            dY_1(:,ii) = Phi_dr_vec_1(2,nodes(ii,:));
            dZ_1(:,ii) = Phi_dr_vec_1(3,nodes(ii,:));
        end
        
        norm_Phi_dr = sqrt( dX.^2 + dY.^2 + dZ.^2 );
        norm_Phi_dr_1 = sqrt( dX_1.^2 + dY_1.^2 + dZ_1.^2 );
        
      
        %%[2] modeごとにplot
        i_fig = 5+i_mode;
        
        
       
        
        h_fig_mode(i_mode) = figure(i_fig);
        set( h_fig_mode(i_mode), 'Position', [1300 600 600 400])
        h_ax(i_ax) = axes( 'Parent', h_fig_mode(i_mode), 'FontSize', 15);
       
        patch( [ X X_1], [ Y Y_1], [ Z Z_1], [ norm_Phi_dr norm_Phi_dr_1], 'Parent', h_ax(i_ax));
        
        light
        lighting gouraud 
        view( h_ax(i_ax), [1 -2 1])
        grid( h_ax(i_ax), 'on')
        xlim( h_ax(i_ax), [0 Length])
        ylim( h_ax(i_ax), [0 Width])
        xlabel( h_ax(i_ax), '{\itX}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
        ylabel( h_ax(i_ax), '{\itY}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
        zlabel( h_ax(i_ax), '{\itZ}^* position [-]', 'FontName', 'Times New Roman', 'FontSize', 15)
        set( h_ax(i_ax), 'FontName', 'Times New Roman')   
        
        h_txt_mode(i_mode) = text( -Length/2, Width, 0, [ '\omega^*_n = ', num2str( omega_a(i_mode), '%0.3f'), ' [-]', sprintf( '\r\n'), '\zeta_n = ', num2str( zeta_n(i_mode), '%0.3f'), ' [-]'],...
                                    'FontSize', 12, 'FontName', 'Times New Roman', 'BackgroundColor', 'g');
        
        i_ax = i_ax + 1;       
        
    end
    
end



%% save

fig_name = { 'nodes', 'displacement', 'snapshot', 'Velocity_field', 'work_rate'};
fig_name_mode = 'mode';

for ii = 1:length( h_fig)
    saveas( h_fig(ii), [ './save/fig/', fig_name{ii}, '.fig']) 
end
for ii = 1:length( h_fig_mode)
    saveas( h_fig_mode(ii), [ './save/fig/modes/', fig_name_mode, '_', num2str( ii, '%d'), '.fig']) 
end

%% create movie

if strcmp( 'avi', movie_format)
    movie2avi( movie_data, './save/data.avi')
elseif strcmp( 'wmv', movie_format)
    video.frames = movie_data;
    video.width=size( video.frames(1).cdata, 2);
    video.height=size( video.frames(1).cdata, 1);
    video.times = movie_time;
    
    mmwrite( './save/data.wmv', video)
end



%% Computing time

disp( [ 'Computing time of structure: ', num2str( measure_time_struct, '%0.2f'), ' [s]'])
disp( [ 'Computing time of fluid: ', num2str( measure_time_fluid, '%0.2f'), ' [s]'])


%% Finish
warndlg( 'Finish')



