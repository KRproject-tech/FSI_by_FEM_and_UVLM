%% Wakeパネルノード点上の流速を算出+Wake生成


    
if i_wake_time == 1    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 計算開始直後  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%[2-0] 後縁パネルノード座標取得 [-]
    i_trail = N_element-Ny+1:N_element;
    N_trail = length( i_trail);                                             %% 後縁パネル数 [-]
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);                             %% Y=0 [-]におけるWakeノード


    %%[2-1] 平板上の渦パネルによる後縁パネルノードへの誘導速度計算 [-]
    r_trail_vec = [ r_panel_vec_2_end;
                    r_panel_vec_31_end];       
    V_gamma_wake = V_wake_func( r_trail_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    

    V_wake_2 = V_gamma_wake(1:end-1,:) + V_in(i_trail,:) - dt_r_panel_vec_2(i_trail,:);
    V_wake_31 = V_gamma_wake(end,:) + V_in(i_trail(1),:) - dt_r_panel_vec_3(i_trail(1),:);


    %%[2-2-0] Wakeパネルノード座標生成 [-]    
    r_wake_1 = r_panel_vec_2_end;
    r_wake_4 = r_panel_vec_3_end;
    r_wake_2 = r_panel_vec_2_end + V_wake_2*d_t_wake;
    r_wake_31 = r_panel_vec_31_end + V_wake_31*d_t_wake;
    r_wake_3 = [ r_wake_31;
                 r_wake_2(1:end-1,:)];

    h_r_wake(i_wake_time) = { [ r_wake_1;
                                r_wake_2;
                                r_wake_3;
                                r_wake_4] };
                            
    %%[2-2-1] Wakeノード移流速度 [-]
    dt_r_wake_1 = dt_r_panel_vec_2(i_trail,:);
    dt_r_wake_4 = dt_r_panel_vec_3(i_trail,:);
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    dt_r_wake_3 = [ dt_r_wake_31;                                           %% Y=0 [m]におけるWakeノード�Bの移流速度
                    dt_r_wake_2(1:end-1,:)];                                %% その他のWakeノード�Bの移流速度

    
    %%[2-3] 循環の保存
    Gamma_wake = Gamma_trail;
    
    h_Gamma_wake(i_wake_time) = { Gamma_wake };
    
    h_r_end(:,:,i_wake_time) = [ r_panel_vec_31_end;
                                 r_panel_vec_2_end];
    h_V_wake_end(:,:,i_wake_time) = [ V_wake_31(1,:);
                                      V_wake_2(1:Ny,:)];

elseif i_wake_time > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2回目の計算以降  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%[*] 閾値以降では変形させない (後縁部のwakeの切り取りによる変動の防止)
    idx_R_wake_x_threshold_no_change = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold_no_change);    	%% R_wake_x_threshold_no_change: 変形を許容するWake末端位置の閾値 [m]
    if ~isempty( idx_R_wake_x_threshold_no_change)

        idx_R_wake_x_threshold_no_change = floor( idx_R_wake_x_threshold_no_change(1)/Ny)*Ny + 1;
        idx_R_wake_x_threshold_no_change_31 = floor( idx_R_wake_x_threshold_no_change(1)/Ny) + 1;
    end
    
    
    

    N_wake = size( r_wake_2, 1);                                            %% Wakeパネル数 [-]

    %%[2-4] Wakeノード座標
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);
  
    %%[*] 平板後縁のノード点も含める． [-]
    old_r_wake_2 = [ r_panel_vec_2_end;             
                     r_wake_2];
    old_r_wake_31 = [ r_panel_vec_31_end;
                      r_wake_3(1:Ny:end,:)];                                %% Y=0 [m]におけるWakeノード
    r_wake_23 = [ old_r_wake_2;
                  old_r_wake_31];
              
              
    %%[2-5] 平板後縁の移動速度 [-]    
    ii_trail = 1:N_trail;
    N_wake_trail = N_wake+N_trail; 
    dt_r_panel_vec_2_wake = zeros(N_wake_trail,3);
    dt_r_panel_vec_31_wake = zeros(N_wake_trail/Ny,3);
    dt_r_panel_vec_2_wake(ii_trail,:) = dt_r_panel_vec_2(i_trail,:);
    dt_r_panel_vec_31_wake(1,:) = dt_r_panel_vec_3(i_trail(1),:);
    
    %%[2-6] Wakeパネルノード座標生成 [-]   
    %%[*] 陽的Euler法では不安定
    %%% H. ABEDI, Development of Vortex Filament Method for Aerodynamic
    %%% Loads on Rotor Blades, Doctor thesis of Chalmers University of Technology, p. 33, 2013.  
    
        
    
    %%[*] 【k1】
    %%% RK4 method    
    %%[*] 平板上の渦パネルによるWakeパネルノードへの誘導速度計算 [-]           
    V_gamma_wake_k1 = V_wake_func( r_wake_23, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] WakeパネルによるWakeパネルノードへの誘導速度計算 [-]     
    V_gamma_w_wake_k1 = V_wake_func( r_wake_23, r_wake_1, r_wake_2, r_wake_3, r_wake_4, Gamma_wake, var_param, 0);
         
    %%[*] 平板後縁 + Wakeパネルノードでの流速 [-]
    V_in_wake = ones(N_wake_trail,1)*U_in*[ 1 0 0]; 
    V_wake_2_k1 = V_gamma_wake_k1(1:N_wake_trail,:) + V_gamma_w_wake_k1(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k1 = V_gamma_wake_k1(N_wake_trail+1:end,:) + V_gamma_w_wake_k1(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;

        
        
    %%[*] 【k2】
    r_wake_2_k2 = old_r_wake_2 + V_wake_2_k1*d_t_wake/2;
    r_wake_31_k2 = old_r_wake_31 + V_wake_31_k1*d_t_wake/2;   
    r_wake_3_k2 = zeros(N_wake_trail,3);
    r_wake_3_k2(1:Ny:end,:) = r_wake_31_k2;                                   %% Y=0 [-]におけるWakeノード�B
    idx_r2 = 1:N_wake_trail;
    idx_r2(1:Ny:end) = [];
    r_wake_3_k2(idx_r2,:) = r_wake_2_k2(idx_r2-1,:);                          %% その他のWakeノード�B
    r_wake_1_k2 = [ r_panel_vec_2_end;
                   r_wake_2_k2(1:end-Ny,:)];
    r_wake_4_k2 = [ r_panel_vec_3_end;
                   r_wake_3_k2(1:end-Ny,:)];  
               
    %%[*] 循環の保存
    %%% Kuttaの条件を満たすように放出渦を決める．
    Gamma_wake = [ Gamma_trail;
                   Gamma_wake];
             

    %%[*] 平板上の渦パネルによるWakeパネルノードへの誘導速度計算 [-]                       
    r_wake_23_k2 = [ r_wake_2_k2;
                    r_wake_31_k2];
    V_gamma_wake_k2 = V_wake_func( r_wake_23_k2, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] WakeパネルによるWakeパネルノードへの誘導速度計算 [-]     
    V_gamma_w_wake_k2 = V_wake_func( r_wake_23_k2, r_wake_1_k2, r_wake_2_k2, r_wake_3_k2, r_wake_4_k2, Gamma_wake, var_param, 0);
                 
    %%[*] 平板後縁 + Wakeパネルノードでの流速 [-]
    V_wake_2_k2 = V_gamma_wake_k2(1:N_wake_trail,:) + V_gamma_w_wake_k2(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k2 = V_gamma_wake_k2(N_wake_trail+1:end,:) + V_gamma_w_wake_k2(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] 【k3】    
    %%[*] 平板後縁 + Wakeパネルノード [-]
    r_wake_2_k3 = old_r_wake_2 + V_wake_2_k2*d_t_wake/2;
    r_wake_31_k3 = old_r_wake_31 + V_wake_31_k2*d_t_wake/2;
    
    r_wake_3_k3 = zeros(N_wake_trail,3);
    r_wake_3_k3(1:Ny:end,:) = r_wake_31_k3;                                       %% Y=0 [-]におけるWakeノード�B
    r_wake_3_k3(idx_r2,:) = r_wake_2_k3(idx_r2-1,:);                              %% その他のWakeノード�B
    r_wake_1_k3 = [ r_panel_vec_2_end;
                    r_wake_2_k3(1:end-Ny,:)];
    r_wake_4_k3 = [ r_panel_vec_3_end;
                    r_wake_3_k3(1:end-Ny,:)];        
                

    r_wake_23_k3 = [ r_wake_2_k3;
                    r_wake_31_k3];
    V_gamma_wake_k3 = V_wake_func( r_wake_23_k3, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] WakeパネルによるWakeパネルノードへの誘導速度計算 [-]     
    V_gamma_w_wake_k3 = V_wake_func( r_wake_23_k3, r_wake_1_k3, r_wake_2_k3, r_wake_3_k3, r_wake_4_k3, Gamma_wake, var_param, 0);
                
    %%[*] 平板後縁 + Wakeパネルノードでの流速 [-]
    V_wake_2_k3 = V_gamma_wake_k3(1:N_wake_trail,:) + V_gamma_w_wake_k3(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k3 = V_gamma_wake_k3(N_wake_trail+1:end,:) + V_gamma_w_wake_k3(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] 【k4】   
    %%[*] 平板後縁 + Wakeパネルノード [-]
    r_wake_2_k4 = old_r_wake_2 + V_wake_2_k3*d_t_wake;
    r_wake_31_k4 = old_r_wake_31 + V_wake_31_k3*d_t_wake;
    
    r_wake_3_k4 = zeros(N_wake_trail,3);
    r_wake_3_k4(1:Ny:end,:) = r_wake_31_k4;                                       %% Y=0 [-]におけるWakeノード�B
    r_wake_3_k4(idx_r2,:) = r_wake_2_k4(idx_r2-1,:);                              %% その他のWakeノード�B
    r_wake_1_k4 = [ r_panel_vec_2_end;
                    r_wake_2_k4(1:end-Ny,:)];
    r_wake_4_k4 = [ r_panel_vec_3_end;
                    r_wake_3_k4(1:end-Ny,:)];        
                

    r_wake_23_k4 = [ r_wake_2_k4;
                    r_wake_31_k4];
    V_gamma_wake_k4 = V_wake_func( r_wake_23_k4, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] WakeパネルによるWakeパネルノードへの誘導速度計算 [-]     
    V_gamma_w_wake_k4 = V_wake_func( r_wake_23_k4, r_wake_1_k4, r_wake_2_k4, r_wake_3_k4, r_wake_4_k4, Gamma_wake, var_param, 0);
                 
    %%[*] 平板後縁 + Wakeパネルノードでの流速 [-]
    V_wake_2_k4 = V_gamma_wake_k4(1:N_wake_trail,:) + V_gamma_w_wake_k4(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k4 = V_gamma_wake_k4(N_wake_trail+1:end,:) + V_gamma_w_wake_k4(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
    
    
     
    
    %%[*] 時間発展 (RK4)
    V_wake_2 = (V_wake_2_k1 + 2*V_wake_2_k2 + 2*V_wake_2_k3 + V_wake_2_k4)/6;
    V_wake_31 = (V_wake_31_k1 + 2*V_wake_31_k2 + 2*V_wake_31_k3 + V_wake_31_k4)/6;
    
    
    %%[*] 閾値以降では変形させない (後縁部のwakeの切り取りによる変動の防止)
    if ~isempty( idx_R_wake_x_threshold_no_change)
        
        idx_2 = idx_R_wake_x_threshold_no_change:N_wake_trail;
        idx_31 = idx_R_wake_x_threshold_no_change_31:size( old_r_wake_31, 1);
        
        V_wake_2(idx_2,:) = V_in_wake(idx_2,:);
        V_wake_31(idx_31,:) = V_in_wake(idx_31,:);
    end   
    
    
    

    r_wake_2 = old_r_wake_2 + V_wake_2*d_t_wake;
    r_wake_31 = old_r_wake_31 + V_wake_31*d_t_wake;
    
    r_wake_3 = zeros(N_wake_trail,3);
    r_wake_3(1:Ny:end,:) = r_wake_31;                                       %% Y=0 [m]におけるWakeノード�B
    r_wake_3(idx_r2,:) = r_wake_2(idx_r2-1,:);                              %% その他のWakeノード�B
    r_wake_1 = [ r_panel_vec_2_end;
             	 r_wake_2(1:end-Ny,:)];
    r_wake_4 = [ r_panel_vec_3_end;
           	     r_wake_3(1:end-Ny,:)];        
                
                
                
                

    h_r_wake(i_wake_time) = { [ r_wake_1;
                                r_wake_2;
                                r_wake_3;
                                r_wake_4] }; 
                            
    %%[*] Wakeノード移流速度 [-]
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    
    dt_r_wake_3 = zeros(N_wake_trail,3);
    dt_r_wake_3(1:Ny:end,:) = dt_r_wake_31;                                	%% Y=0 [m]におけるWakeノード�Bの移流速度
    dt_r_wake_3(idx_r2,:) = dt_r_wake_2(idx_r2-1,:);                      	%% その他のWakeノード�Bの移流速度
    dt_r_wake_1 = [ dt_r_panel_vec_2(i_trail,:);
                    dt_r_wake_2(1:end-Ny,:)];
    dt_r_wake_4 = [ dt_r_panel_vec_3(i_trail,:);
                    dt_r_wake_3(1:end-Ny,:)];                           
                            
                            
    h_Gamma_wake(i_wake_time) = { Gamma_wake };

    h_r_end(:,:,i_wake_time) = [ r_panel_vec_31_end;
                                 r_panel_vec_2_end];
    h_V_wake_end(:,:,i_wake_time) = [ dt_r_wake_31(1,:);
                                      dt_r_wake_2(1:Ny,:)];

end




%% Wake末端削除(計算コスト削減)





idx_R_wake_x_threshold = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold);     %% R_wake_x_threshold: Wake末端位置の閾値 [m]
if ~isempty( idx_R_wake_x_threshold)
    idx_R_wake_x_threshold = floor( idx_R_wake_x_threshold(1)/Ny)*Ny + 1;

    r_wake_1(idx_R_wake_x_threshold:end,:) = [];
    r_wake_2(idx_R_wake_x_threshold:end,:) = [];
    r_wake_3(idx_R_wake_x_threshold:end,:) = [];
    r_wake_4(idx_R_wake_x_threshold:end,:) = [];
    
    dt_r_wake_1(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_2(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_3(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_4(idx_R_wake_x_threshold:end,:) = [];

    Gamma_wake(idx_R_wake_x_threshold:end) = [];
end




