%% Function: Read CM1 output data

function data_struct = read_cm1_output(cm1_file)

    %%% Create data structure

    % Mesh grid
    data_struct.xh = double(ncread(cm1_file, 'xh'));
    data_struct.yh = double(ncread(cm1_file, 'yh'));
    data_struct.time = ncread(cm1_file, 'time');

    % Surface pressure [Pa]
    try
        psfc = ncread(cm1_file, 'psfc');  p_mean = mean(psfc, 'all');
    catch
        psfc = ncread(cm1_file, 'prs');  p_mean = mean(psfc, 'all');
    end

    % Pressure perturbation [Pa]
    % Remove spatial mean at each time
    pp = psfc - p_mean;  pp = pp - mean(pp, [1,2]);
    data_struct.pp = pp;
    
end
