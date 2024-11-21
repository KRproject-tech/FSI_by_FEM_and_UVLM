%% äOêœ
function out = cross_mat( a, b)

three_N_element = size( a, 2);

i_vec = [ 2:3:three_N_element;
          3:3:three_N_element; 
          1:3:three_N_element];
i_vec = reshape( i_vec, 1, []);      
      
j_vec = [ 3:3:three_N_element;
          1:3:three_N_element; 
          2:3:three_N_element];
j_vec = reshape( j_vec, 1, []);            


out = a(:,i_vec).*b(:,j_vec) - a(:,j_vec).*b(:,i_vec);

end