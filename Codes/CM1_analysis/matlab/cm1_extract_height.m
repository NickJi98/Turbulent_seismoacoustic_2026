%% Extract output at specific height level

% Height level (m)
fprintf('Running with hq = %f\n', hq);

% Output directory 
cm1_dir = './';

% Output file indices
istart = 8;  iend = 3607;
Nt = iend-istart+1;

% Variable names
var_list = {'prs', 'u', 'v', 'th'};

%% Horizontal & Vertical grids

disp(['Extract outputs from ', num2str(Nt), ' files.']);

% Example file to obtain vertical grid
ref_file = fullfile(cm1_dir, ['cm1out_', sprintf('%06d', istart), '.nc']);
ref_info = ncinfo(ref_file);

% Horizontal grids
xh = ncread(ref_file, 'xh');  xf = ncread(ref_file, 'xf');
yh = ncread(ref_file, 'yh');  yf = ncread(ref_file, 'yf');

% Vertical grids
zh = ncread(ref_file, 'zh');  Nh = length(zh);
zf = ncread(ref_file, 'zf');  Nf = length(zf);

% Index of query height level
hq = hq / 1e3;
dz = zh(2) - zh(1);
zh_mask = (zh <= hq + dz*2/3) & (zh >= hq - dz*2/3);
zf_mask = (zf <= hq + dz*2/3) & (zf >= hq - dz*2/3);

%% Extract variables

% Parallel setup
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

for var_ = var_list
    
    % Vertical grid of current variable
    varname = var_{1};
    var_mat = squeeze(ncread(ref_file, varname));
    size_ = size(var_mat);

    if size_(3) == Nh
        hq_mask = zh_mask;
        disp(zh(hq_mask));
    elseif size_(3) == Nf
        hq_mask = zf_mask;
        disp(zf(hq_mask))
    end

    % Horizontal grids
    if size_(1) == Nh
        x = xh;  xname = 'xh';
    elseif size_(1) == Nf
        x = xf;  xname = 'xf';
    end

    if size_(2) == Nh
        y = yh;  yname = 'yh';
    elseif size_(2) == Nf
        y = yf;  yname = 'yf';
    end

    % Initialize array
    combine_mat = zeros(size_(1), size_(2), Nt);
    combine_mat(:, :, 1) = mean(var_mat(:, :, hq_mask), 3);

    % Read each output file
    parfor i = 2:Nt
        nc_filepath = fullfile(cm1_dir, ['cm1out_', sprintf('%06d', istart+i-1), '.nc']);
        var_mat = squeeze(ncread(nc_filepath, varname));
        combine_mat(:, :, i) = mean(var_mat(:, :, hq_mask), 3);
    end

    %%% Save to NetCDF file %%%
    outfile = fullfile(cm1_dir, ['cm1out_', varname, '_', num2str(hq*1e3), 'm.nc']);
    ncid = netcdf.create(outfile, 'NETCDF4');
    
    % Define dimensions
    dimid_x = netcdf.defDim(ncid, xname, length(x));
    dimid_y = netcdf.defDim(ncid, yname, length(y));
    dimid_time = netcdf.defDim(ncid, 'time', Nt);
    
    % Define variables
    varid_x = netcdf.defVar(ncid, xname, 'NC_FLOAT', dimid_x);
    varid_y = netcdf.defVar(ncid, yname, 'NC_FLOAT', dimid_y);
    varid_time = netcdf.defVar(ncid, 'time', 'NC_FLOAT', dimid_time);
    varid_var = netcdf.defVar(ncid, varname, 'NC_FLOAT', [dimid_x, dimid_y, dimid_time]);
    
    % Set compression for each variable
    compressionLevel = 8;
    netcdf.defVarDeflate(ncid, varid_x, true, true, compressionLevel);
    netcdf.defVarDeflate(ncid, varid_y, true, true, compressionLevel);
    netcdf.defVarDeflate(ncid, varid_time, true, true, compressionLevel);
    netcdf.defVarDeflate(ncid, varid_var, true, true, compressionLevel);
    netcdf.endDef(ncid);
    
    % Write data to variables
    netcdf.putVar(ncid, varid_x, single(x));
    netcdf.putVar(ncid, varid_y, single(y));
    netcdf.putVar(ncid, varid_time, single(0:Nt-1));
    netcdf.putVar(ncid, varid_var, single(combine_mat));
    
    % Add global attribute
    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'hq', hq);
    
    % Close the NetCDF file
    netcdf.close(ncid);

end
