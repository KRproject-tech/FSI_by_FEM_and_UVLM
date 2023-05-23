%% Vortex lattice method �̊i�q���� 


%% [0] �p�l���m�[�h�_���W [-]( r_i = Sc(x_i,y_i)*q_i, i��{1,...,N_element} )

Sc_mat_v_panel_1 = zeros(3,N_q,N_element);
Sc_mat_v_panel_4 = Sc_mat_v_panel_1;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% ii�v�f�̒��� [-]
    dW = dW_vec(ii);                        %% ii�v�f�̕� [-]
            
    %%[*] �v�f���W�̃m�[�h�_ [-]
    x_i = dL*[ 1/4  1/4];   %% �_�@�C�C
    y_i = dW*[ 1    0];     %% �_�@�C�C

    Sc_mat_v_panel_1(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_4(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end

%% [1] �����[�m�[�h�_���W [m]( r_i = Sc(x_i,y_i)*q_i, i��{1,...,N_element} )

Sc_mat_v_panel_end_2 = zeros(3,N_q,N_element);
Sc_mat_v_panel_end_3 = Sc_mat_v_panel_end_2;
%%[*] �v�f�ԍ���Y�����ɑ�������̂ŁC�����[�̗v�f�ԍ���; N_element-Ny+1~N_element�D
for ii = N_element-Ny+1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% ii�v�f�̒��� [-]
    dW = dW_vec(ii);                        %% ii�v�f�̕� [-]
            
    %%[*] �v�f���W�̃m�[�h�_ [-]
    x_i = dL*[ 1  1];   %% �_�A�C�B
    y_i = dW*[ 1  0];   %% �_�A�C�B

    Sc_mat_v_panel_end_2(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_end_3(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end


%% [2] �O���[�o���s��g��

Sc_mat_panel_global_1 = sparse(3*N_element,N_q_all);
Sc_mat_panel_global_4 = Sc_mat_panel_global_1;
for ii = 1:N_element

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_global_1(j_vec,i_vec) = Sc_mat_panel_global_1(j_vec,i_vec) + squeeze( Sc_mat_v_panel_1(:,:,ii));
    Sc_mat_panel_global_4(j_vec,i_vec) = Sc_mat_panel_global_4(j_vec,i_vec) + squeeze( Sc_mat_v_panel_4(:,:,ii));
end



%%[*] �v�f�ԍ���Y�����ɑ�������̂ŁC�����[�̗v�f�ԍ���; N_element-Ny+1~N_element�D
Sc_mat_panel_end_global_2 = sparse(3*N_element,N_q_all);
Sc_mat_panel_end_global_3 = Sc_mat_panel_end_global_2;
for ii = N_element-Ny+1:N_element

    %% 1�m�[�h������9���� ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T �� R^9 )
    %% 1�v�f������36�����@( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T �� R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_end_global_2(j_vec,i_vec) = Sc_mat_panel_end_global_2(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_2(:,:,ii));
    Sc_mat_panel_end_global_3(j_vec,i_vec) = Sc_mat_panel_end_global_3(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_3(:,:,ii));
end


%%[2-0] �_�A [m]�@(�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
Sc_mat_panel_global_2 = [ Sc_mat_panel_global_1(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_2(end-3*Ny+1:end,:) - Sc_mat_panel_global_1(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_1(end-3*Ny+1:end,:)];

%%[2-1] �_�B [m] (�����[�̃p�l���m�[�h�_�́C�����[�m�[�h���W�����)
Sc_mat_panel_global_3 = [ Sc_mat_panel_global_4(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_3(end-3*Ny+1:end,:) - Sc_mat_panel_global_4(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_4(end-3*Ny+1:end,:)];

                      
%% [3] �R���P�[�V�����_�v�Z [-]
% Sc_mat_col_global = ( Sc_mat_panel_global_1 + Sc_mat_panel_global_2 + Sc_mat_panel_global_3 + Sc_mat_panel_global_4 )/4;

%%[*] �s�ϓ��v�f���ɑΉ�
Sc_mat_col_global = sparse(3*N_element,N_q_all);
for ii = 1:N_element    
     
    j_vec = 3*ii-2:3*ii;                                %% [ rx ry rz]^T
    
    if ii <= N_element-Ny
        dL_i = dL_vec(ii);                              %% ii�v�f�̒��� [-]
        dL_ip1 = dL_vec(ii+Ny);                         %% ii�v�f�̌���̗v�f�̒��� [-]
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + diag( [ (dL_i + dL_ip1)/(3*dL_i + dL_ip1)     1/2     (dL_i + dL_ip1)/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/2 ...
                                           	+ diag( [ 2*dL_i/(3*dL_i + dL_ip1)              1/2     2*dL_i/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) )/2;
    else
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + ( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/4;
    end
end


%% [4] �@���x�N�g���Z�o [-]

Sc_mat_31 = Sc_mat_panel_global_3 - Sc_mat_panel_global_1;
Sc_mat_24 = Sc_mat_panel_global_2 - Sc_mat_panel_global_4;



