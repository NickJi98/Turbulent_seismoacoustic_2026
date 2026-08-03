%% Analysis of CM1 output: Temporal evolution

% Output directory 
cm1_dir = './';

% Time interval for recording
diagfrq = 60.0;     % diagnostic output

%% Get variable names to analyze temporal evolution

% Make directory to save data
data_dir = fullfile(cm1_dir, 'Data');
if exist(data_dir, 'dir') ~= 7
    mkdir(data_dir);
    disp(['Directory ', data_dir, ' created successfully.']);
else
    disp(['Directory ', data_dir, ' already exists.']);
end

% All NC diagnostic files
filelist = dir(fullfile(cm1_dir, 'cm1out_diag*.nc'));
Nfile = length(filelist);

% Starting file number
i0 = regexp(filelist(1).name, 'cm1out_diag_(\d+).nc', 'tokens', 'once');
i0 = str2double(i0{1});
clear filelist;

% Example file to obtain variable names & vertical grid
ref_file = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', i0), '.nc']);
ref_info = ncinfo(ref_file);

% List of variable names (whose data size is 1)
evo_var = {};
for i = 1:length(ref_info.Variables)
    data_size = ref_info.Variables(i).Size;

    if all(data_size == 1)
        evo_var{end+1} = ref_info.Variables(i).Name;
    end
end

% List of variable names
disp(['Number of variables: ', num2str(length(evo_var))]);
disp(strjoin(evo_var, ', '));

%% Read output data

% Create empty struct 
evo_struct = struct;

% Parallel setup
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

% Read variable
for i = 1:length(evo_var) 
    varname = evo_var{i};
    var_mat = zeros(Nfile, 1);

    parfor j = i:Nfile
        nc_file = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', j+i0-1), '.nc']);
        var_mat(j) = squeeze(ncread(nc_file, varname));
    end

    evo_struct.(varname) = var_mat;
end
disp('Finish reading variables.');

%% Calculate PBL depth

pbl_depth = zeros(Nfile, 1);

parfor j = 1:Nfile
    nc_file = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', j+i0-1), '.nc']);
    
    % Vertical grid & Surface stress
    zf = squeeze(ncread(nc_file, 'zf')) ./ 1e3;
    sf_tau = squeeze(ncread(nc_file, 'stau'));

    % Relevant quantities
    ufr = squeeze(ncread(nc_file, 'ufr'));
    ufd = squeeze(ncread(nc_file, 'ufd'));
    ufs = squeeze(ncread(nc_file, 'ufs'));

    vfr = squeeze(ncread(nc_file, 'vfr'));
    vfd = squeeze(ncread(nc_file, 'vfd'));
    vfs = squeeze(ncread(nc_file, 'vfs'));

    try
        ufw = squeeze(ncread(nc_file, 'ufw'));
        vfw = squeeze(ncread(nc_file, 'vfw'));
    catch
        if j == 1; disp('No wall model variable ufw, vfw. Skipped.'); end
        ufw = zeros(size(ufr));  vfw = zeros(size(vfr));
    end

    % Vertical momentum flux
    ttau = sqrt((ufr + ufd + ufs + ufw).^2 + (vfr + vfd + vfs + vfw).^2);

    % PBL depth (from vertical momentum flux)
    ind_ = find(ttau < 0.05*sf_tau, 1);
    pbl_depth(j) = - sf_tau * zf(ind_) / (ttau(ind_) - sf_tau);

end

% Save PBL depth to struct
evo_struct.pbl_depth = pbl_depth;
disp('Finish calculating PBL depth.');

%% Save data

% Save data to matfile
matfile = 'evolution.mat';
save(fullfile(cm1_dir, 'Data', matfile), "evo_struct");