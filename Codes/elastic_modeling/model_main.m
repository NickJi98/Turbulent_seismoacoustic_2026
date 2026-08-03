%% Seismic modeling from CM1 output (Layered media or halfspace)

% Check deta directory, matlab functions in `init_env.m' before running
run("init_env.m");

% If you do not want netcdf outputs, modify the last section.

%% Pressure source

%%% Data file %%%
cm1_file = fullfile(cm1_dir, 'cm1out_prs.nc');
fprintf('Read in data: %s\n\n', cm1_file);

%%% Read pressure source %%%
src = read_cm1_output(cm1_file);

%% Elastic structure

% Layer properties (rho, vp, vs, thickness)
% Units: g/cm^3, km/s, km/s, km

% Read from file
% elast_prop = readmatrix('vel_model_4layer.csv', 'NumHeaderLines', 1);
elast_prop = readmatrix(vel_file, 'NumHeaderLines', 1);

% Vs30 model for elastic halfspace
% elast_prop = [1.92, 1.45, 0.27, 0.1; ...
%     1.92, 1.45, 0.27, 0];     % Vs30, TA.645A
% elast_prop = [1.906, 1.407, 0.244, 0.1; ...
%     1.906, 1.407, 0.244, 0];  % Vs30, TA.544A

% Print layered model
Nlayer = size(elast_prop, 1) - 1;
row_names = cellstr(num2str((1:Nlayer)'))';
row_names{end+1} = 'Halfspace';
disp('Layered media:');
disp(array2table(elast_prop, ...
    "VariableNames", {'rho (g/cm^3)', 'Vp (km/s)', 'Vs (km/s)', 'Thickness (km)'}, ...
    "RowNames", row_names));

%% Seismic modeling (Layered medium)

% Numerical solution
[uz_layer, ~, ~] = qs_model(src, elast_prop);

%% Seismic modeling (Halfspace)

% Halfspace property
top_hs_prop = elast_prop(1, 1:3);

% Quasi-static response
[uz_top_hs, ~, ~] = calc_disp_halfspace(src, top_hs_prop, 0);

%% Save variable fields

%%% Save full pressure & vertical displ. fields %%%

save_full = 1;
filename = fullfile(data_dir, 'LES_ref.nc');

% Time step
dt = src.time(2) - src.time(1);

% Full domain
if save_full == 1; save_uz(uz_layer, src, filename, 1); end
