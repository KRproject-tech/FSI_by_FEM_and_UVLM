%% �����󐶐�





%% �v�f����W�_�E���z�x�N�g���ɕϊ� (Nx+1)*(Ny+1)��


%%[0] �v�f����W�_ [ r_x r_y r_z]^T

coordinates_tmp = coordinates;
%%[0-0] 1���� 
coordinates = [ coordinates_tmp(:,1).';                                      	%% r_x0
                coordinates_tmp(:,2).';                                        	%% r_y0
                0*coordinates_tmp(:,1).' + Height/2].';                        	%% r_z0
%%[0-1] 2����
coordinates_1 = [	coordinates_tmp(:,1).';                                    	%% r_x0
                    coordinates_tmp(:,2).';                                    	%% r_y0
                    0*coordinates_tmp(:,1).' - Height/2].';                    	%% r_z0



%%[1-1] ���z�x�N�g�� [ dx_r_x dx_r_y dx_r_z]^T
%%[1-0] 1���� 
dx_r0_vec = [  	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec = dx_r0_vec./norm_mat( dx_r0_vec);                                    %% dx_r (���K��)


%%[1-1] 2���� 
dx_r0_vec_1 = [	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec_1 = dx_r0_vec_1./norm_mat( dx_r0_vec_1);                             	%% dx_r (���K��)










