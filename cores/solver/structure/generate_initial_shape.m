%% 初期状生成





%% 要素上座標点・勾配ベクトルに変換 (Nx+1)*(Ny+1)個


%%[0] 要素上座標点 [ r_x r_y r_z]^T

coordinates_tmp = coordinates;
%%[0-0] 1枚目 
coordinates = [ coordinates_tmp(:,1).';                                      	%% r_x0
                coordinates_tmp(:,2).';                                        	%% r_y0
                0*coordinates_tmp(:,1).' + Height/2].';                        	%% r_z0
%%[0-1] 2枚目
coordinates_1 = [	coordinates_tmp(:,1).';                                    	%% r_x0
                    coordinates_tmp(:,2).';                                    	%% r_y0
                    0*coordinates_tmp(:,1).' - Height/2].';                    	%% r_z0



%%[1-1] 勾配ベクトル [ dx_r_x dx_r_y dx_r_z]^T
%%[1-0] 1枚目 
dx_r0_vec = [  	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec = dx_r0_vec./norm_mat( dx_r0_vec);                                    %% dx_r (正規化)


%%[1-1] 2枚目 
dx_r0_vec_1 = [	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec_1 = dx_r0_vec_1./norm_mat( dx_r0_vec_1);                             	%% dx_r (正規化)










