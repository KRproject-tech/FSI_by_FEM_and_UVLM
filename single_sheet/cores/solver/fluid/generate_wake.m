%% Wake�p�l���m�[�h�_��̗������Z�o+Wake����


    
if i_wake_time == 1    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% �v�Z�J�n����  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%[2-0] �㉏�p�l���m�[�h���W�擾 [-]
    i_trail = N_element-Ny+1:N_element;
    N_trail = length( i_trail);                                             %% �㉏�p�l���� [-]
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);                             %% Y=0 [-]�ɂ�����Wake�m�[�h


    %%[2-1] ����̉Q�p�l���ɂ��㉏�p�l���m�[�h�ւ̗U�����x�v�Z [-]
    r_trail_vec = [ r_panel_vec_2_end;
                    r_panel_vec_31_end];       
    V_gamma_wake = V_wake_func( r_trail_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    

    V_wake_2 = V_gamma_wake(1:end-1,:) + V_in(i_trail,:) - dt_r_panel_vec_2(i_trail,:);
    V_wake_31 = V_gamma_wake(end,:) + V_in(i_trail(1),:) - dt_r_panel_vec_3(i_trail(1),:);


    %%[2-2-0] Wake�p�l���m�[�h���W���� [-]    
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
                            
    %%[2-2-1] Wake�m�[�h�ڗ����x [-]
    dt_r_wake_1 = dt_r_panel_vec_2(i_trail,:);
    dt_r_wake_4 = dt_r_panel_vec_3(i_trail,:);
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    dt_r_wake_3 = [ dt_r_wake_31;                                           %% Y=0 [m]�ɂ�����Wake�m�[�h�B�̈ڗ����x
                    dt_r_wake_2(1:end-1,:)];                                %% ���̑���Wake�m�[�h�B�̈ڗ����x

    
    %%[2-3] �z�̕ۑ�
    Gamma_wake = Gamma_trail;
    
    h_Gamma_wake(i_wake_time) = { Gamma_wake };
    
    h_r_end(:,:,i_wake_time) = [ r_panel_vec_31_end;
                                 r_panel_vec_2_end];
    h_V_wake_end(:,:,i_wake_time) = [ V_wake_31(1,:);
                                      V_wake_2(1:Ny,:)];

elseif i_wake_time > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2��ڂ̌v�Z�ȍ~  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%[*] 臒l�ȍ~�ł͕ό`�����Ȃ� (�㉏����wake�̐؂���ɂ��ϓ��̖h�~)
    idx_R_wake_x_threshold_no_change = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold_no_change);    	%% R_wake_x_threshold_no_change: �ό`�����e����Wake���[�ʒu��臒l [m]
    if ~isempty( idx_R_wake_x_threshold_no_change)

        idx_R_wake_x_threshold_no_change = floor( idx_R_wake_x_threshold_no_change(1)/Ny)*Ny + 1;
        idx_R_wake_x_threshold_no_change_31 = floor( idx_R_wake_x_threshold_no_change(1)/Ny) + 1;
    end
    
    
    

    N_wake = size( r_wake_2, 1);                                            %% Wake�p�l���� [-]

    %%[2-4] Wake�m�[�h���W
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);
  
    %%[*] ���㉏�̃m�[�h�_���܂߂�D [-]
    old_r_wake_2 = [ r_panel_vec_2_end;             
                     r_wake_2];
    old_r_wake_31 = [ r_panel_vec_31_end;
                      r_wake_3(1:Ny:end,:)];                                %% Y=0 [m]�ɂ�����Wake�m�[�h
    r_wake_23 = [ old_r_wake_2;
                  old_r_wake_31];
              
              
    %%[2-5] ���㉏�̈ړ����x [-]    
    ii_trail = 1:N_trail;
    N_wake_trail = N_wake+N_trail; 
    dt_r_panel_vec_2_wake = zeros(N_wake_trail,3);
    dt_r_panel_vec_31_wake = zeros(N_wake_trail/Ny,3);
    dt_r_panel_vec_2_wake(ii_trail,:) = dt_r_panel_vec_2(i_trail,:);
    dt_r_panel_vec_31_wake(1,:) = dt_r_panel_vec_3(i_trail(1),:);
    
    %%[2-6] Wake�p�l���m�[�h���W���� [-]   
    %%[*] �z�IEuler�@�ł͕s����
    %%% H. ABEDI, Development of Vortex Filament Method for Aerodynamic
    %%% Loads on Rotor Blades, Doctor thesis of Chalmers University of Technology, p. 33, 2013.  
    
        
    
    %%[*] �yk1�z
    %%% RK4 method    
    %%[*] ����̉Q�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]           
    V_gamma_wake_k1 = V_wake_func( r_wake_23, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wake�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]     
    V_gamma_w_wake_k1 = V_wake_func( r_wake_23, r_wake_1, r_wake_2, r_wake_3, r_wake_4, Gamma_wake, var_param, 0);
         
    %%[*] ���㉏ + Wake�p�l���m�[�h�ł̗��� [-]
    V_in_wake = ones(N_wake_trail,1)*U_in*[ 1 0 0]; 
    V_wake_2_k1 = V_gamma_wake_k1(1:N_wake_trail,:) + V_gamma_w_wake_k1(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k1 = V_gamma_wake_k1(N_wake_trail+1:end,:) + V_gamma_w_wake_k1(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;

        
        
    %%[*] �yk2�z
    r_wake_2_k2 = old_r_wake_2 + V_wake_2_k1*d_t_wake/2;
    r_wake_31_k2 = old_r_wake_31 + V_wake_31_k1*d_t_wake/2;   
    r_wake_3_k2 = zeros(N_wake_trail,3);
    r_wake_3_k2(1:Ny:end,:) = r_wake_31_k2;                                   %% Y=0 [-]�ɂ�����Wake�m�[�h�B
    idx_r2 = 1:N_wake_trail;
    idx_r2(1:Ny:end) = [];
    r_wake_3_k2(idx_r2,:) = r_wake_2_k2(idx_r2-1,:);                          %% ���̑���Wake�m�[�h�B
    r_wake_1_k2 = [ r_panel_vec_2_end;
                   r_wake_2_k2(1:end-Ny,:)];
    r_wake_4_k2 = [ r_panel_vec_3_end;
                   r_wake_3_k2(1:end-Ny,:)];  
               
    %%[*] �z�̕ۑ�
    %%% Kutta�̏����𖞂����悤�ɕ��o�Q�����߂�D
    Gamma_wake = [ Gamma_trail;
                   Gamma_wake];
             

    %%[*] ����̉Q�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]                       
    r_wake_23_k2 = [ r_wake_2_k2;
                    r_wake_31_k2];
    V_gamma_wake_k2 = V_wake_func( r_wake_23_k2, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wake�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]     
    V_gamma_w_wake_k2 = V_wake_func( r_wake_23_k2, r_wake_1_k2, r_wake_2_k2, r_wake_3_k2, r_wake_4_k2, Gamma_wake, var_param, 0);
                 
    %%[*] ���㉏ + Wake�p�l���m�[�h�ł̗��� [-]
    V_wake_2_k2 = V_gamma_wake_k2(1:N_wake_trail,:) + V_gamma_w_wake_k2(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k2 = V_gamma_wake_k2(N_wake_trail+1:end,:) + V_gamma_w_wake_k2(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] �yk3�z    
    %%[*] ���㉏ + Wake�p�l���m�[�h [-]
    r_wake_2_k3 = old_r_wake_2 + V_wake_2_k2*d_t_wake/2;
    r_wake_31_k3 = old_r_wake_31 + V_wake_31_k2*d_t_wake/2;
    
    r_wake_3_k3 = zeros(N_wake_trail,3);
    r_wake_3_k3(1:Ny:end,:) = r_wake_31_k3;                                       %% Y=0 [-]�ɂ�����Wake�m�[�h�B
    r_wake_3_k3(idx_r2,:) = r_wake_2_k3(idx_r2-1,:);                              %% ���̑���Wake�m�[�h�B
    r_wake_1_k3 = [ r_panel_vec_2_end;
                    r_wake_2_k3(1:end-Ny,:)];
    r_wake_4_k3 = [ r_panel_vec_3_end;
                    r_wake_3_k3(1:end-Ny,:)];        
                

    r_wake_23_k3 = [ r_wake_2_k3;
                    r_wake_31_k3];
    V_gamma_wake_k3 = V_wake_func( r_wake_23_k3, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wake�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]     
    V_gamma_w_wake_k3 = V_wake_func( r_wake_23_k3, r_wake_1_k3, r_wake_2_k3, r_wake_3_k3, r_wake_4_k3, Gamma_wake, var_param, 0);
                
    %%[*] ���㉏ + Wake�p�l���m�[�h�ł̗��� [-]
    V_wake_2_k3 = V_gamma_wake_k3(1:N_wake_trail,:) + V_gamma_w_wake_k3(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k3 = V_gamma_wake_k3(N_wake_trail+1:end,:) + V_gamma_w_wake_k3(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] �yk4�z   
    %%[*] ���㉏ + Wake�p�l���m�[�h [-]
    r_wake_2_k4 = old_r_wake_2 + V_wake_2_k3*d_t_wake;
    r_wake_31_k4 = old_r_wake_31 + V_wake_31_k3*d_t_wake;
    
    r_wake_3_k4 = zeros(N_wake_trail,3);
    r_wake_3_k4(1:Ny:end,:) = r_wake_31_k4;                                       %% Y=0 [-]�ɂ�����Wake�m�[�h�B
    r_wake_3_k4(idx_r2,:) = r_wake_2_k4(idx_r2-1,:);                              %% ���̑���Wake�m�[�h�B
    r_wake_1_k4 = [ r_panel_vec_2_end;
                    r_wake_2_k4(1:end-Ny,:)];
    r_wake_4_k4 = [ r_panel_vec_3_end;
                    r_wake_3_k4(1:end-Ny,:)];        
                

    r_wake_23_k4 = [ r_wake_2_k4;
                    r_wake_31_k4];
    V_gamma_wake_k4 = V_wake_func( r_wake_23_k4, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wake�p�l���ɂ��Wake�p�l���m�[�h�ւ̗U�����x�v�Z [-]     
    V_gamma_w_wake_k4 = V_wake_func( r_wake_23_k4, r_wake_1_k4, r_wake_2_k4, r_wake_3_k4, r_wake_4_k4, Gamma_wake, var_param, 0);
                 
    %%[*] ���㉏ + Wake�p�l���m�[�h�ł̗��� [-]
    V_wake_2_k4 = V_gamma_wake_k4(1:N_wake_trail,:) + V_gamma_w_wake_k4(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k4 = V_gamma_wake_k4(N_wake_trail+1:end,:) + V_gamma_w_wake_k4(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
    
    
     
    
    %%[*] ���Ԕ��W (RK4)
    V_wake_2 = (V_wake_2_k1 + 2*V_wake_2_k2 + 2*V_wake_2_k3 + V_wake_2_k4)/6;
    V_wake_31 = (V_wake_31_k1 + 2*V_wake_31_k2 + 2*V_wake_31_k3 + V_wake_31_k4)/6;
    
    
    %%[*] 臒l�ȍ~�ł͕ό`�����Ȃ� (�㉏����wake�̐؂���ɂ��ϓ��̖h�~)
    if ~isempty( idx_R_wake_x_threshold_no_change)
        
        idx_2 = idx_R_wake_x_threshold_no_change:N_wake_trail;
        idx_31 = idx_R_wake_x_threshold_no_change_31:size( old_r_wake_31, 1);
        
        V_wake_2(idx_2,:) = V_in_wake(idx_2,:);
        V_wake_31(idx_31,:) = V_in_wake(idx_31,:);
    end   
    
    
    

    r_wake_2 = old_r_wake_2 + V_wake_2*d_t_wake;
    r_wake_31 = old_r_wake_31 + V_wake_31*d_t_wake;
    
    r_wake_3 = zeros(N_wake_trail,3);
    r_wake_3(1:Ny:end,:) = r_wake_31;                                       %% Y=0 [m]�ɂ�����Wake�m�[�h�B
    r_wake_3(idx_r2,:) = r_wake_2(idx_r2-1,:);                              %% ���̑���Wake�m�[�h�B
    r_wake_1 = [ r_panel_vec_2_end;
             	 r_wake_2(1:end-Ny,:)];
    r_wake_4 = [ r_panel_vec_3_end;
           	     r_wake_3(1:end-Ny,:)];        
                
                
                
                

    h_r_wake(i_wake_time) = { [ r_wake_1;
                                r_wake_2;
                                r_wake_3;
                                r_wake_4] }; 
                            
    %%[*] Wake�m�[�h�ڗ����x [-]
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    
    dt_r_wake_3 = zeros(N_wake_trail,3);
    dt_r_wake_3(1:Ny:end,:) = dt_r_wake_31;                                	%% Y=0 [m]�ɂ�����Wake�m�[�h�B�̈ڗ����x
    dt_r_wake_3(idx_r2,:) = dt_r_wake_2(idx_r2-1,:);                      	%% ���̑���Wake�m�[�h�B�̈ڗ����x
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




%% Wake���[�폜(�v�Z�R�X�g�팸)





idx_R_wake_x_threshold = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold);     %% R_wake_x_threshold: Wake���[�ʒu��臒l [m]
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



