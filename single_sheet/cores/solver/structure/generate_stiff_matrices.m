%% �����s��̑g��

%%[*] Gauss-Legendre���ς̏d�݂��s�񉻁F1:1~N_q, 2:1~N_q, 3:�̂̐ϕ����_, 4:�ł̐ϕ����_, 5:��q�Ȃ������v�Z���鎞�̐����ԍ�
ww_mat = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ N_q 1 1 1]);
ww_mat2 = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ N_q N_q 1 1]);
ww_mat3 = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ 1 1 1 1]);
ww_mat4 = repmat( permute( w_vec, [ 3 4 1 2]), [ N_q N_q 1 1]);



%%[*] ������Gauss-Legendre���ς̐ϕ����_�����R�s�[�F1:1~N_q, 2:1~N_q, 3:�̂̐ϕ����_, 4:�ł̐ϕ����_, 5:��q�Ȃ������v�Z���鎞�̐����ԍ�
Dp_mat2 = repmat( Dp_mat, [ 1 1 length( p_vec) length( p_vec)]);


%% [0] �����Ђ��� (Qe�� = t��_A (��q_��)^t*Dm*�� dA �� R^36)
if flag_output == 1     
    
    int_e_Dp_e_vec = zeros(1,N_element);

    Qe_eps_vec_i = zeros(N_q,N_element);
    dq_Qe_eps_mat_i = zeros(N_q,N_q,N_element);                           	%% ��_A (dqdq_��^T*D*�� + dq_��^T*D*dq_��) dA 
    Qd_eps_mat_i = zeros(N_q,N_q,N_element);                            	%% ��_A dq_��^T*D*dq_�� dA (dq_��*dt_q = dt_��)    

 
    for ii = 1:N_element

        dL = dL_vec(ii);                        %% ii�v�f�̒��� [-]
        dW = dW_vec(ii);                        %% ii�v�f�̕� [-]


        %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
        %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
        i_vec = i_vec_v{ii};
        q_i_vec = q_vec(i_vec);

        dx_n_sc_struct = dx_n_Sc_struct(ii);       


        [ dq_eps_Dp_eps_vec, dq_eps_Dp_dq_eps_mat, dqdq_eps_Dp_eps_P_dq_eps_Dp_dq_eps_mat, eps_Dp_eps] = dq_eps_Dm_eps( q_i_vec, Dp_mat2, dx_n_sc_struct, p_vec);
        int_dq_e_Dp_e = dL*dW/4*ww_mat.*dq_eps_Dp_eps_vec;
        int_dq_e_Dp_e = sum( sum( int_dq_e_Dp_e, 4), 3).';        

        Qe_eps_vec_i(:,ii) = int_dq_e_Dp_e;                                 %% ii�v�f�̊O�͍s��

        int_e_Dp_e = dL*dW/4*ww_mat3.*eps_Dp_eps;
        int_e_Dp_e_vec(ii) = sum( sum( int_e_Dp_e, 4), 3).';
        

        int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e = dL*dW/4*ww_mat2.*dqdq_eps_Dp_eps_P_dq_eps_Dp_dq_eps_mat;
        int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e = sum( sum( int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e, 4), 3).';

        dq_Qe_eps_mat_i(:,:,ii) = int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e;       	%% ii�v�f�̊O�͍s��

        if theta_a ~= 0                                                     %% �\�����������݂���ꍇ�̂݌v�Z (�v�Z�R�X�g�팸)
            int_dq_e_Dp_dq_e = dL*dW/4*ww_mat2.*dq_eps_Dp_dq_eps_mat;
            int_dq_e_Dp_dq_e = sum( sum( int_dq_e_Dp_dq_e, 4), 3).';

            Qd_eps_mat_i(:,:,ii) = int_dq_e_Dp_dq_e;                        %% ii�v�f�̊O�͍s��
        end    


    end
end

%% [1] �Ȃ��Ђ��� (Qe�� = ��_A (��q_��)^t*D��*�� dA �� R^36)
int_k_Dp_k_vec = zeros(1,N_element);

Qe_k_vec_i = zeros(N_q,N_element);
Qd_k_mat_i = zeros(N_q,N_q,N_element);                                      %% ��_A dq_��^T*D*dq_�� dA (dq_��*dt_q = dt_��)


%%[1-1] ��q�Ȃ������v�Z����Ƃ��̔����ʁF1:1~N_q, 2:1~N_q, 3:�̂̐ϕ����_, 4:�ł̐ϕ����_, 5:��q�Ȃ������v�Z���鎞�̐����ԍ�
% eye_N_q = repmat( eye(N_q), [ 1 1 length( p_vec)  length( p_vec)]);
% d_q = 1e-10;
% dq_vec = d_q*eye_N_q;
% dq_vec = permute( dq_vec, [ 5 1 3 4 2]);
% dq_vec = reshape( dq_vec, 3, [], length( p_vec), length( p_vec), N_q);	 %% �`��֐��̍s���ɍ��킹��(�ς̌v�Z�̂���)�D1����:���W����(x,y,z)�ɑΉ��C2:����(r,dx_r,dy_r)�ɑΉ��C3,4����:�K�E�X���ϓ_�C5����:��q�Ȃ������v�Z���鎞�̐����ԍ�



for ii = 1:N_element
    
    dL = dL_vec(ii);                %% ii�v�f�̒��� [-]
    dW = dW_vec(ii);                %% ii�v�f�̕� [-]
    
    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);    
    
    dx_n_sc_struct = dx_n_Sc_struct(ii);            

    [ dq_k_Dp_k_vec, dq_k_Dp_dq_k_mat, k_Dp_k] = dq_k_Dk_k_FAST( q_i_vec, Dp_mat2, dx_n_sc_struct, p_vec, theta_a);
    
    int_k_Dp_k = dL*dW/4*ww_mat3.*k_Dp_k;
    int_k_Dp_k = sum( sum( int_k_Dp_k, 4), 3).'; 
    int_k_Dp_k_vec(ii) = int_k_Dp_k;
    
    int_dq_k_Dp_k = dL*dW/4*ww_mat.*dq_k_Dp_k_vec;
    int_dq_k_Dp_k = sum( sum( int_dq_k_Dp_k, 4), 3).';  
    Qe_k_vec_i(:,ii) = int_dq_k_Dp_k;                                       %% ii�v�f�̊O�͍s��
    
    if theta_a ~= 0                                                         %% �\�����������݂���ꍇ�̂݌v�Z (�v�Z�R�X�g�팸)
        int_dq_k_Dp_dq_k = dL*dW/4*ww_mat2.*dq_k_Dp_dq_k_mat;
        int_dq_k_Dp_dq_k = sum( sum( int_dq_k_Dp_dq_k, 4), 3).';
        
        Qd_k_mat_i(:,:,ii) = int_dq_k_Dp_dq_k;
    end 

end





%% [2] ��]�_���p [-] (Qd�� = -��_��W dq_(n/||n||)^T*dx_r/||dx_r||*��q_��0y dy �� R^36x36)
 

Qd_theta_mat_i = zeros(N_q,N_q,N_element);                                  %% �� dq(n/|n|)^T*(dx_r/|dx_r|)*dq_��_0y dy 
dt_theta_0y = zeros(size( element_C_theta));

for ii = element_C_theta

    dW = dW_vec(ii);                %% ii�v�f�̕� [-]

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec); 
    dt_q_i_vec = dt_q_vec(i_vec);   
    

    dx_n_sc_struct = dx_n_Sc_struct(ii);     

    [ q_d_theta_mat, dq_theta_0y_v] = dq_n_dx_r_dq_theta_0y( q_i_vec, dx_n_sc_struct, p_vec);

    int_q_d_theta_mat = dW/2*ww_mat4.*q_d_theta_mat;                        %% x�����͒萔 (x=0 [-])�Cy�����ɐϕ� 
    int_q_d_theta_mat = sum( sum( int_q_d_theta_mat, 4), 3);   
    
    int_dq_theta_0y_v = dW/2*ww_mat4(1,:,:,:).*dq_theta_0y_v;
    int_dq_theta_0y_v = sum( sum( int_dq_theta_0y_v, 4), 3);
    

    Qd_theta_mat_i(:,:,ii) = int_q_d_theta_mat;
    
    dt_theta_0y(ii) = int_dq_theta_0y_v*dt_q_i_vec;
end





%% [3] ��]�����s��̑�� [-] (-��_��W dq_(n/||n||)^T*dx_r/||dx_r||*dt��q_��0y dy �� R^36x36)
 

J2_mat_i = zeros(N_q,N_q,N_element);                                  %% �� dq(n/|n|)^T*(dx_r/|dx_r|)*dq_��_0y dy 


for ii = element_C_theta

    dW = dW_vec(ii);                %% ii�v�f�̕� [-]

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);  
    dt_q_i_vec = dt_q_vec(i_vec);    

    dx_n_sc_struct = dx_n_Sc_struct(ii);     

    J2_mat = dq_n_dx_r_dtdq_theta_0y( q_i_vec, dt_q_i_vec, dx_n_sc_struct, p_vec);

    int_J2_mat = dW/2*ww_mat4.*J2_mat;                        %% x�����͒萔 (x=0 [-])�Cy�����ɐϕ� 
    int_J2_mat = sum( sum( int_J2_mat, 4), 3);           

    J2_mat_i(:,:,ii) = int_J2_mat;
end





%% [4] �O���[�o���s��g��

%%[3-0] �e����
Qe_global = zeros(N_q_all,1);
Qk_global = zeros(N_q_all,1);
for ii = 1:N_element

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};    

    if flag_output == 1 
        Qe_global(i_vec,1) = Qe_global(i_vec,1) + zeta_m*Qe_eps_vec_i(:,ii);
    end
    Qk_global(i_vec,1) = Qk_global(i_vec,1) + eta_m*Qe_k_vec_i(:,ii);
end


%%[4-1] �����s�� + �\������
dq_Qe_global = sparse(N_q_all,N_q_all);
Qd_global = sparse(N_q_all,N_q_all);

dq_Qe_global_vec = zeros(1,length( dq_Qe_eps_mat_i(:)));
Qd_global_vec = zeros(1,length( Qd_eps_mat_i(:)));
if flag_output == 1 
    
    for ii = 1:N_element

        %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
        %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
        i_vec = i_vec_v{ii};        
        
        dq_Qe_global_vec_i = zeta_m*dq_Qe_eps_mat_i(:,:,ii);
        dq_Qe_global_vec((ii - 1)*N_q^2+(1:N_q^2)) = dq_Qe_global_vec_i(:).';
        
%         dq_Qe_global(i_vec,i_vec) = dq_Qe_global(i_vec,i_vec) + squeeze( zeta_m*dq_Qe_eps_mat_i(:,:,ii));

        if theta_a ~= 0                                                         %% �\�����������݂���ꍇ�̂݌v�Z (�v�Z�R�X�g�팸)
            
%             Qd_global(i_vec,i_vec) = Qd_global(i_vec,i_vec) + squeeze( theta_a*zeta_m*Qd_eps_mat_i(:,:,ii) +  theta_a*eta_m*Qd_k_mat_i(:,:,ii));
            
            Qd_global_vec_i = theta_a*zeta_m*Qd_eps_mat_i(:,:,ii) +  theta_a*eta_m*Qd_k_mat_i(:,:,ii);
            Qd_global_vec((ii - 1)*N_q^2+(1:N_q^2)) = Qd_global_vec_i(:).';
        end     
    end

    dq_Qe_global = setsparse( dq_Qe_global, I_vec, J_vec, dq_Qe_global_vec, @plus);
    if theta_a ~= 0                                                         %% �\�����������݂���ꍇ�̂݌v�Z (�v�Z�R�X�g�팸)
        Qd_global = setsparse( Qd_global, I_vec, J_vec, Qd_global_vec, @plus);
    end
    
end


%%[4-2] ��]����
Qd_theta_global = sparse(N_q_all,N_q_all);

    
for ii = element_C_theta

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};

    Qd_theta_global(i_vec,i_vec) = Qd_theta_global(i_vec,i_vec) + squeeze( C_theta_a*Qd_theta_mat_i(:,:,ii));
end


%%[4-3] ��]����: J1*dt^2_q + J2*dt_q 
J_global_1 = sparse(N_q_all,N_q_all);    
for ii = element_C_theta

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii};  

    J_global_1(i_vec,i_vec) = J_global_1(i_vec,i_vec) + squeeze( J_a*Qd_theta_mat_i(:,:,ii));
end

J_global_2 = sparse(N_q_all,N_q_all);    
for ii = element_C_theta

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = i_vec_v{ii}; 

    J_global_2(i_vec,i_vec) = J_global_2(i_vec,i_vec) + squeeze( J_a*J2_mat_i(:,:,ii));
end

    
