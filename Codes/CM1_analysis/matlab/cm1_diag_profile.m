%% Analysis of CM1 output: Averaged vertical profile

% Output directory 
cm1_dir = './';

% Time interval for recording
diagfrq = 60.0;     % diagnostic output

% Time range for average [hour]
time_average = [5, 6];

%% Get variable names to analyze vertical profiles

% Make directory to save data
data_dir = fullfile(cm1_dir, 'Data');
if exist(data_dir, 'dir') ~= 7
    mkdir(data_dir);
    disp(['Directory ', data_dir, ' created successfully.']);
else
    disp(['Directory ', data_dir, ' already exists.']);
end

% NC diagnostic file index within time range
istart = time_average(1) * 3600 / diagfrq + 1;
iend = time_average(2) * 3600 / diagfrq + 1;
disp(['Extract evolution of variables from ', num2str(iend-istart+1), ' files.']);

% Example file to obtain variable names & vertical grid
ref_file = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', istart), '.nc']);
ref_info = ncinfo(ref_file);

% Vertical grid
zh = ncread(ref_file, 'zh') ./ 1e3;  Nh = length(zh);
zf = ncread(ref_file, 'zf') ./ 1e3;  Nf = length(zf);

% List of variable names
scalar_var = {};  wlev_var = {};
for i = 1:length(ref_info.Variables)
    data_size = ref_info.Variables(i).Size;

    if prod(data_size) == Nh
        scalar_var{end+1} = ref_info.Variables(i).Name;
    elseif prod(data_size) == Nf
        wlev_var{end+1} = ref_info.Variables(i).Name;
    end
end
disp(['Number of variables: ', num2str(length(scalar_var) + length(wlev_var)), '.']);
disp(strjoin(scalar_var, ', '));  disp(strjoin(wlev_var, ', '));

% Remove variable names for vertical grid
scalar_var = setdiff(scalar_var, {'zh'});
wlev_var = setdiff(wlev_var, {'zf'});

%% Read output data 

% Create empty struct 
scalar_struct = struct;  wlev_struct = struct;
scalar_struct.zh = zh;   wlev_struct.zf = zf;

% Parallel setup
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

% Read variable (scalar)
for i = 1:length(scalar_var) 
    var_mat = zeros(Nh, iend-istart+1);
    varname = scalar_var{i};

    parfor j = 1:iend-istart+1
        nc_filepath = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', istart+j-1), '.nc']);
        var_mat(:, j) = squeeze(ncread(nc_filepath, varname));
    end

    scalar_struct.(varname) = mean(var_mat, 2);
end
disp('Finish reading scalar variables.');

% Read variable (w level)
for i = 1:length(wlev_var) 
    var_mat = zeros(Nf, iend-istart+1);
    varname = wlev_var{i};

    parfor j = 1:iend-istart+1
        nc_filepath = fullfile(cm1_dir, ['cm1out_diag_', sprintf('%06d', istart+j-1), '.nc']);
        var_mat(:, j) = squeeze(ncread(nc_filepath, varname));
    end

    wlev_struct.(varname) = mean(var_mat, 2);
end
disp('Finish reading w-level variables.');

%% Save data

% Save data to matfile
matfile = sprintf('diag_profile_%d.mat', time_average(end));
save(fullfile(cm1_dir, 'Data', matfile), ...
    "scalar_struct", "wlev_struct", "time_average");
