%% version
exe_ver = 12.0;

%% version check

display( [ 'Solver version: ', num2str( exe_ver, '%0.1f')])
if ~exist( 'param_ver', 'var')   
    param_ver = nan;
    display( 'Parameter file version: no definition')
else    
    display( [ 'Parameter file version: ', num2str( param_ver, '%0.1f')])
end

if exe_ver ~= param_ver
   
    warndlg( 'This parameter file is not compatible!!')
end

pause( 1.0);