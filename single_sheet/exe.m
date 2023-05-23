clc
clear all
close all hidden

%% delete
delete( '*.asv')
delete( '*.log')

%% path
add_pathes

%% parameter
param_setting

%% version check
version_check


%% multi threading
maxNumCompThreads( core_num);


%% exe

%%[0] Definition of shape functions
generate_shape_function;

%%[1] Mesh generation
%%[1-0] FEMメッシュ生成
generate_elements;
%%[1-1] Vortex lattice method の格子生成
generate_panel;

%%[1-2] Matrix generation
generate_matrices;

%%[1-3] 外乱入力
generate_Qf_time_mat;




%%[2] 時間発展 [-]

time_m = 0:d_t:End_Time;
initial_values;


%%[3] solve modes
solve_mode;


%%[4] 計算速度確認
if speed_check == 1
    time_m = time_m(1:10);
    profile on;
end

tic
i_time = 1;
i_time_cnt = 1;
i_wake_time = 1;
time_fluid = 0;
d_t_wake = d_t*dt_wake_per_dt;
time_wake_m = 0;
fluid_compute_flag = 1;
time = 0;
while time <= time_m(end) || ~fluid_compute_flag
    
    time = i_time*d_t;
    
    disp( [ 'Time = ', num2str( time, '%0.4f'), ' [-]'])

    %%[5] 構造解析
    solve_structure;    
   
    
    %%[6] 流体解析       
    if mod( i_time, dt_wake_per_dt) == 1       
        drawnow                                         	%% 処理中のフリーズ防止 
        
        if fluid_compute_flag
                        
            %% 1step前の値を更新
            old_Qf_p_global = Qf_p_global;
            old_Qf_p_mat_global = Qf_p_mat_global;
            old_Qf_p_mat0_global = Qf_p_mat0_global;
            old_Qf_p_lift2_mat_global = Qf_p_lift2_mat_global;
        
            solve_fluid;    
            i_wake_time = i_wake_time + 1;   
        else
            
            %% 1step前の値を更新
            Qf_p_global_a = Qf_p_global;
            Qf_p_mat_global_a = Qf_p_mat_global;
            Qf_p_mat0_global_a = Qf_p_mat0_global;
            Qf_p_lift2_mat_global_a = Qf_p_lift2_mat_global;
                           
            time_fluid = time;                                  %% 流体パラメータ更新時刻 [-] (流体力の時間方向補間用)
        end
                

        disp( [ 'Ma = ', num2str( Ma), ' [-], Ua = ', num2str( Ua), ' [-], theta_a = ', num2str( theta_a), ' [-], C_theta_a = ', num2str( C_theta_a), ' [-], J_a = ', num2str( J_a), ' [-]'])
    end
    
    
    %%[7] エネルギ収支の評価
    solve_energy;
    
    
    %%[8] 解析データ中途保存
    if mod( i_time, 500) == 0 && ~fluid_compute_flag
        save ./save/NUM_DATA 
    end
    
    %%[9] 反復計算
    if mod( i_time, dt_wake_per_dt) == 1      
        
        if fluid_compute_flag
            
            i_time = i_time - i_time_cnt;            
            fluid_compute_flag = 0;
        else
            
            i_time_cnt = 0;
            fluid_compute_flag = 1;
        end
    end
    
    i_time = i_time + 1;
    i_time_cnt = i_time_cnt + 1;
end
toc

%% save 

%%[*] 計算速度確認
if speed_check == 1
    profile viewer;
else
    save ./save/NUM_DATA
end