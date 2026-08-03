%% Create database for simulations

% Directory to analyze
maindir = '/path/to/Data/CM1/param_search';

% Field names: Inputs
propVars = {'name'};
basicVars = {'nx', 'dx', 'cm1setup', 'testcase', 'imoist', 'sgsmodel', ...
    'idiss', 'isnd', 'iwnd'};
sfcVars = {'isfcflx', 'sfcmodel', 'tsk0', 'xland0', 'lu0', 'season', ...
    'cecd', 'cnstce', 'cnstcd', 'set_flx', 'cnst_shflx', 'cnst_lhflx', ...
    'set_znt', 'cnst_znt'};
hurrVars = {'hurr_vg', 'hurr_rad', 'hurr_vgpl', 'hurr_rotate'};

% Fields names: Results
outVars = {'S10', 'ust', 'ust_tau', 'inflow', 'hpbl'};

%% Loop over simulation cases

% Summarize field names
allVars = [propVars, basicVars, sfcVars, hurrVars, outVars];
inputVars = [basicVars, sfcVars, hurrVars];

% Initialize an empty structure
Database = struct();
for j = 1:length(allVars)
    Database.(allVars{j}) = [];
end

% All subfolders
subdirs = dir(maindir);
subdirs = subdirs([subdirs.isdir]);
subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));
[~, inds_] = sort({subdirs.name});
subdirs = subdirs(inds_);

for i = 1:length(subdirs)

    % Path to the input file
    dirPath = fullfile(maindir, subdirs(i).name);
    filePath = fullfile(dirPath, 'namelist.input');

    % Check input file
    if ~isfile(filePath)
        warning('Input file not found in %s', subdirs(i).name);
        continue;
    end

    % Input filename
    fileID = fopen(filePath, 'r');
    Database.('name'){end+1} = subdirs(i).name;

    %%% Read input parameters %%%
    j = 1;  varName = inputVars{j};
    pattern = sprintf('%s\\s*=\\s*(-?[\\d\\.]+|\\w+)', varName);

    while ~feof(fileID)
        
        % Regular expression to match the pattern 'varName = value'
        line = fgetl(fileID);
        tokens = regexp(line, pattern, 'tokens');
        if isempty(tokens)
            continue;
        end

        % Input parameters always numeric
        paramValue = tokens{1}{1};
        paramValue = str2double(paramValue);
        Database.(varName)(end+1) = paramValue;
        j = j + 1;
        
        % Finish reading relevant parameters
        if j > length(inputVars)
            break;
        end

        % Update varName
        varName = inputVars{j};
        pattern = sprintf('%s\\s*=\\s*(-?[\\d\\.]+|\\w+)', varName);
    end
    fclose(fileID);

    %%% Read simulation results %%%
    try
        load(fullfile(dirPath, 'evolution.mat'));
        load(fullfile(dirPath, 'diag_profile_6.mat'));
    catch
        for j = 1:length(outVars)
            Database.(outVars{j})(end+1) = NaN;
        end
        warning('Diagnostic files not found in %s', subdirs(i).name);
        continue;
    end
    command = sprintf('grep "diagfrq " %s | awk -F "=" ''{print $2}''', filePath);
    [~, result] = system(command);
    diagfrq = str2double(strtrim(result));

    % Time range for averaging
    istart = time_average(1) * 3600 / diagfrq + 1;
    iend = time_average(2) * 3600 / diagfrq + 1;
    
    % Vertical grids
    zh = scalar_struct.zh;  zf = wlev_struct.zf;
    
    % Obtain 10 m height
    dz = zh(2) - zh(1);  
    zh_mask = (zh <= 0.01 + dz*2/3) & (zh >= 0.01 - dz*2/3);
    zf_mask = (zf <= 0.01 + dz*2/3) & (zf >= 0.01 - dz*2/3);

    % Diagnostic variables
    wsp = scalar_struct.wsp;
    u = scalar_struct.u;  v = scalar_struct.v;
    ufr = wlev_struct.ufr + wlev_struct.ufd;  
    vfr = wlev_struct.vfr + wlev_struct.vfd;
    try
        ufs = wlev_struct.ufs + wlev_struct.ufw;  
        vfs = wlev_struct.vfs + wlev_struct.vfw;
    catch
        ufs = wlev_struct.ufs;  vfs = wlev_struct.vfs;
    end

    % Intermediate results
    ttau = sqrt((ufr + ufs).^2 + (vfr + vfs).^2);
    sf_tau = mean(evo_struct.stau(istart:iend));
    ind_ = find(ttau < 0.05*sf_tau, 1);

    % Analyzed results
    S10 = mean(wsp(zh_mask));
    ust = mean(evo_struct.ust(istart:iend));
    ust_tau = sqrt(ttau(1));
    inflow = -rad2deg(atan(mean(u(zh_mask)) / mean(v(zh_mask))));
    h_pbl = - sf_tau * zf(ind_) / (ttau(ind_) - sf_tau);

    % Write to database
    Database.S10(end+1) = S10;
    Database.ust(end+1) = ust;
    Database.ust_tau(end+1) = ust_tau;
    Database.inflow(end+1) = inflow;
    Database.hpbl(end+1) = h_pbl;
    
    % Reset variables
    clear time_average evo_struct scalar_struct wlev_struct;

end

%% Save database

matFileName = fullfile(maindir, 'param_search.mat');
save(matFileName, 'Database');
fprintf('Database saved to %s\n', matFileName);
