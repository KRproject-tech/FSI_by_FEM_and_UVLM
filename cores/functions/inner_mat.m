%% ŠOÏ
function out = inner_mat( a, b)



out0 = a(:,1:3:end).*b(:,1:3:end) + a(:,2:3:end).*b(:,2:3:end) + a(:,3:3:end).*b(:,3:3:end);

out = kron( out0, ones(1,3));  


end