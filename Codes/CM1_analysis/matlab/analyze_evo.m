%% Analyze temporal evolution

run("init_env.m");

% Illustrate LES reaching a quasi-steady state after 6-hr simulation.

%% Load data

% Read evolution of variables
load(fullfile(sim_dir, 'evolution.mat'));

% Read input wind speed (hurr_vg)
input_file = fullfile('./CM1_Output', sim_dir, 'Data', 'namelist.input');
command = sprintf('grep "hurr_vg " %s | awk -F "=" ''{print $2}''', input_file);
[status, result] = system(command);
if status == 0
    input_V = str2double(strtrim(result));
else
    error('Failed to read hurr_vg from namelist.input!');
end

% Time interval for recording diagnostic output
diagfrq = 60.0;

%% Common variables

t_hr = evo_struct.time ./ 3600;  Nt = length(t_hr);

%% Plot temporal evolutions

% Parameters
x_lim = [0, max(t_hr)];  h_lim = [0, 2];

%%% Wind Speed %%%
figure('Name', 'Wind', 'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot(t_hr, evo_struct.v10, 'k-');
xlabel('Time (hour)');  ylabel('V10 (m/s)');
xlim(x_lim);
subplot(1,3,2);         plot(t_hr, evo_struct.u10, 'k-');
xlabel('Time (hour)');  ylabel('U10 (m/s)');
xlim(x_lim);
subplot(1,3,3);         plot(t_hr, evo_struct.s10, 'k-');
xlabel('Time (hour)');  ylabel('S10 (m/s)');
xlim(x_lim);

%%% Boundary Layer Height %%%
figure('Name', 'PBL', 'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot(t_hr, evo_struct.pbl_depth, 'k-');
xlabel('Time (hour)');  ylabel('hpbl (km)');
xlim(x_lim);            ylim(h_lim);
subplot(1,3,2);         plot(t_hr, evo_struct.zi./1e3, 'k-');
xlabel('Time (hour)');  ylabel('zi (km)');
xlim(x_lim);            ylim(h_lim);
subplot(1,3,3);         plot(t_hr, evo_struct.zwspmax./1e3, 'k-');
xlabel('Time (hour)');  ylabel('zwspmax (km)');
xlim(x_lim);            ylim(h_lim);

%% Read statistics outputs

% Statistics output
stats_file = fullfile(sim_dir, 'cm1out_stats.nc');
stats_info = ncinfo(stats_file);
stats_var = {};
for i = 1:length(stats_info.Variables)
    data_size = stats_info.Variables(i).Size;

    if data_size > 1
        stats_var{end+1} = stats_info.Variables(i).Name;
    end
end

stats_struct = struct;
for i = 1:length(stats_var)
    varname = stats_var{i};
    stats_struct.(varname) = squeeze(ncread(stats_file, varname));
end

%% Plot statistics

%%% Wind Speed %%%
figure('Name', 'Wind', 'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot(t_hr, stats_struct.sumax, 'k-');  hold on;
plot(t_hr, stats_struct.sumin, 'k--');
xlabel('Time (hour)');  ylabel('Surface U (m/s)');  xlim(x_lim);
subplot(1,3,2);         plot(t_hr, stats_struct.svmax, 'k-');  hold on;
plot(t_hr, stats_struct.svmin, 'k--');
xlabel('Time (hour)');  ylabel('Surface V (m/s)');  xlim(x_lim);
subplot(1,3,3);         plot(t_hr, stats_struct.wsp10max, 'k-');  hold on;
plot(t_hr, stats_struct.wsp10min, 'k-');
xlabel('Time (hour)');  ylabel('S10 (m/s)');  xlim(x_lim);

%%% Pressure Perturb., Max Wind Level, Eddy Diffusivity %%%
figure('Name', 'Wind', 'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot(t_hr, stats_struct.ppmax, 'k-');  hold on;
plot(t_hr, stats_struct.ppmin, 'k--');
xlabel('Time (hour)');  ylabel('Pressure Perturbation (Pa)');  xlim(x_lim);
subplot(1,3,2);         plot(t_hr, stats_struct.zwspmax, 'k-');
xlabel('Time (hour)');  ylabel('Level of Max Wind (m)');  xlim(x_lim);
subplot(1,3,3);         plot(t_hr, stats_struct.kmvmax, 'k-');
xlabel('Time (hour)');  ylabel('Vert. Eddy Diffusivity (m^2/s)');  xlim(x_lim);
