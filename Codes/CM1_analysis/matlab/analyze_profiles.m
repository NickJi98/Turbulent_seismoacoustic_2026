%% Analyze averaged vertical profile

run("init_env.m");

% Visualize vertical profiles of LES outputs.

%% Load data

% Plot base state
if_base = true;

% Load data
load(fullfile(sim_dir, 'diag_profile_6.mat'));
load(fullfile(sim_dir, 'evolution.mat'));

% Time interval for recording diagnostic output
diagfrq = 60.0;

%% Calculate profiles

% Time range for averaging
istart = time_average(1) * 3600 / diagfrq + 1;
iend = time_average(2) * 3600 / diagfrq + 1;

% Vertical grids
zh = scalar_struct.zh;  zf = wlev_struct.zf;

% Usual variables
wsp = scalar_struct.wsp;
upup = wlev_struct.upup_w;  vpvp = wlev_struct.vpvp_w;  
wpwp = wlev_struct.wpwp_w;

ufr = wlev_struct.ufr + wlev_struct.ufd;  
vfr = wlev_struct.vfr + wlev_struct.vfd;
try
    ufs = wlev_struct.ufs + wlev_struct.ufw;  
    vfs = wlev_struct.vfs + wlev_struct.vfw;
catch
    ufs = wlev_struct.ufs;  vfs = wlev_struct.vfs;
end

rtke = wlev_struct.rtke_w;  stke = wlev_struct.stke;

rtau = sqrt(ufr.^2 + vfr.^2);  stau = sqrt(ufs.^2 + vfs.^2);
ttau = sqrt((ufr + ufs).^2 + (vfr + vfs).^2);
ttau_w = interp1(zf, ttau, zh, 'linear');

% Effective eddy diffusivity
dudz = wlev_struct.dudz;  dvdz = wlev_struct.dvdz;
rKeff = rtau ./ sqrt(dudz.^2 + dvdz.^2);
try
    sKeff = wlev_struct.kmv + wlev_struct.kmw;
catch
    sKeff = wlev_struct.kmv;
end

% Effective mixing length
Leff = sqrt((rKeff+sKeff) ./ sqrt(dudz.^2 + dvdz.^2));

%% Read base-state profiles

if if_base
    base_var = [fieldnames(scalar_struct); fieldnames(wlev_struct)];
    base_file = fullfile(sim_dir, 'base_state.nc');
    
    base_struct = struct;
    for i = 1:length(base_var)
        varname = base_var{i};
        base_struct.(varname) = squeeze(ncread(base_file, varname));
    end
end

%% Calculate parameters

% Obtain 10 m height
dz = zh(2) - zh(1);  
zh_mask = (zh <= 0.01 + dz*2/3) & (zh >= 0.01 - dz*2/3);
zf_mask = (zf <= 0.01 + dz*2/3) & (zf >= 0.01 - dz*2/3);

% Inflow angle (10 m height)
inflow_theta = -rad2deg(atan(mean(scalar_struct.u(zh_mask)) / mean(scalar_struct.v(zh_mask))));

% PBL depth (from vertical momentum flux)
sf_tau = mean(evo_struct.stau(istart:iend));
ind_ = find(ttau < 0.05*sf_tau, 1);
h_pbl = - sf_tau * zf(ind_) / (ttau(ind_) - sf_tau);

% Friction velocity
ust = sqrt(ttau(1));

% Read input wind speed (hurr_vg)
input_file = fullfile(sim_dir, 'namelist.input');
command = sprintf('grep "hurr_vg " %s | awk -F "=" ''{print $2}''', input_file);
[~, result] = system(command);
input_V = str2double(strtrim(result));
if isnan(input_V)
    base_v = base_struct.v;
    input_V = base_v(2) - (base_v(2)-base_v(1))/dz * zh(2);
end

% Print parameters
fprintf('\n');
fprintf('Input hurr_vg: %g m/s\n', input_V);
fprintf('Inflow angle:  %.2f deg\n', inflow_theta);
fprintf('Wind at 10 m:  %.2f m/s\n', mean(wsp(zh_mask)));
fprintf('PBL height:    %.2f km\n', h_pbl);
fprintf('Friction vel.: %.2f m/s\n', ust);
fprintf('\n');

%% Plot profiles (Wind Speed)

% Height limit for plotting (km)
h_lim = [0, 2];

%%% Wind Speed %%%
figure('Name', 'Wind', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);
if if_base; plot(base_struct.u, zh, 'k-'); hold on; end
plot(scalar_struct.u, zh, 'r-');
xlabel('<u> (m/s)');    ylabel('z (km)');  ylim(h_lim);
subplot(1,3,2);
if if_base; plot(base_struct.v, zh, 'k-'); hold on; end
plot(scalar_struct.v, zh, 'r-');
scatter(input_V, 0, 400, 'b', 'filled', 'pentagram');
xlabel('<v> (m/s)');    ylabel('z (km)');  ylim(h_lim);
subplot(1,3,3);
if if_base; plot(base_struct.wsp, zh, 'k-'); hold on; end
plot(scalar_struct.wsp, zh, 'r-');  hold on;
scatter(input_V, 0, 400, 'b', 'filled', 'pentagram');
xlabel('<U> (m/s)');    ylabel('z (km)');  ylim(h_lim);
if if_base; legend('Base Profile', 'Final Profile', 'Location', 'Best'); end

%%% Log-Wind Profile %%%
kappa = 0.41;  z0 = 0.3;  d = 0;
z_lw = linspace(0, 0.05, 1e2);
u = (ust/kappa) .* log((z_lw.*1e3-d)./z0);
plot(u, z_lw, 'b--', 'DisplayName', 'Log Wind Profile');

%%% Wind Speed (Normalized) %%%
figure('Name', 'Wind (Normalized)', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot(scalar_struct.u ./ ust, zh ./ h_pbl, 'r-');  
xlabel('<u>/u_*');      ylabel('z/h');  ylim(h_lim ./ h_pbl);
subplot(1,3,2);         plot(scalar_struct.v ./ ust, zh ./ h_pbl, 'r-');  
xlabel('<v>/u_*');      ylabel('z/h');  ylim(h_lim ./ h_pbl);
subplot(1,3,3);         plot(scalar_struct.wsp ./ ust, zh ./ h_pbl, 'r-');  
xlabel('<U>/u_*');      ylabel('z/h');  ylim(h_lim ./ h_pbl);

%% Plot profiles (TKE)

%%% Wind Variance %%%
figure('Name', 'Wind Variance', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);
if if_base; plot(base_struct.upup_w, zf, 'k-'); hold on; end
plot(upup + 2/3.*stke, zf, 'r-');
xlabel("<u'u'> (m^2/s^2)");    ylabel('z (km)');  ylim(h_lim);
subplot(1,3,2);
if if_base; plot(base_struct.vpvp_w, zf, 'k-'); hold on; end
plot(vpvp + 2/3.*stke, zf, 'r-');
xlabel("<v'v'> (m^2/s^2)");    ylabel('z (km)');  ylim(h_lim);
subplot(1,3,3);
if if_base; plot(base_struct.stke, zf, 'k-'); hold on; end
plot(wpwp + 2/3.*stke, zf, 'r-');  hold on;
xlabel("<w'w'>, TKE (m^2/s^2)");    ylabel('z (km)');  
xlim([0, 20]);          ylim(h_lim);
plot(rtke + stke, zf, 'r--');

%%% Wind Variance (Normalized) %%%
figure('Name', 'Wind Variance (Normalized)', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,3,1);         plot((upup+2/3.*stke) ./ ust^2, zf ./ h_pbl, 'r-');  
xlabel("<u'u'>/u_*^2"); ylabel('z/h');  ylim(h_lim ./ h_pbl);
subplot(1,3,2);         plot((vpvp+2/3.*stke) ./ ust^2, zf ./ h_pbl, 'r-');  
xlabel("<v'v'>/u_*^2"); ylabel('z/h');  ylim(h_lim ./ h_pbl);
subplot(1,3,3);         plot((wpwp+2/3.*stke) ./ ust^2, zf ./ h_pbl, 'r-');  hold on;
xlabel("<w'w'>/u_*^2, TKE/u_*^2");    ylabel('z/h');  
xlim([0, 5]);           ylim(h_lim ./ h_pbl);
plot((rtke + stke) ./ ust^2, zf ./ h_pbl, 'r--');

%% Plot profiles (Temperature & Specific humidity)

%%% Potential Temperature %%%
figure('Name', 'theta & qv', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);
if if_base; plot(base_struct.th, zh, 'k-'); hold on; end
plot(scalar_struct.th, zh, 'r-');
xlabel("<\theta> (K)"); ylabel('z (km)');
xlim([295, 315]);       ylim([0, 2.5]);
if if_base; legend('Base Profile', 'Final Profile', 'Location', 'Best'); end

subplot(1,2,2);
try
    if if_base; plot(base_struct.qv.*1e3, zh, 'k-'); hold on; end
    plot(scalar_struct.qv.*1e3, zh, 'r-');
    xlabel("<q_v> (g/kg)"); ylabel('z (km)');
    xlim([10, 22]);       ylim([0, 2.5]);
    if if_base; legend('Base Profile', 'Final Profile', 'Location', 'Best'); end
catch
    warning('No moisture included in simulation');
end

%% Plot profiles (Subgrid Ratio)

%%% Subgrid TKE & Ratio %%%
figure('Name', 'Subgrid TKE', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);         plot(rtke + stke, zf, 'k-');  hold on;
plot(rtke, zf, 'r--');  plot(stke, zf, 'b--');
xlabel("TKE (m^2/s^2)");ylabel('z (km)');
xlim([0, 20]);          ylim(h_lim);
legend('Total', 'Resolved', 'Subgrid');
subplot(1,2,2);         plot(stke ./ (stke+rtke), zf, 'k-');
xlabel("Ratio (Subgrid / Total)");  ylabel('z (km)');
xlim([0, 0.5]);         ylim(h_lim);

%%% Subgrid Vertical Flux & Ratio %%%
figure('Name', 'Subgrid Vertical Flux', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);         plot(ttau, zf, 'k-');  hold on;
plot(rtau, zf, 'r--');  plot(stau, zf, 'b--');
xlabel("|\tau|/\rho (m^2/s^2)");ylabel('z (km)');
xlim([0, 4]);           ylim(h_lim);
legend('Total', 'Resolved', 'Subgrid');
subplot(1,2,2);         plot(stau ./ ttau, zf, 'k-');
xlabel("Ratio (Subgrid / Total)");  ylabel('z (km)');
xlim([0, 0.5]);         ylim(h_lim);

%% Plot profiles (Vertical Momentum Flux)

figure('Name', 'Vertical Flux', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);         plot(ttau ./ ust^2, zf, 'k-');  hold on;
xlabel("|\tau| / \rhou_*^2"); ylabel('z (km)');
xlim([0, 1.1]);         ylim(h_lim);
subplot(1,2,2);         plot(ttau_w ./ wsp.^2, zh, 'k-');
xlabel("|\tau| / \rho<U>^2");   ylabel('z (km)');  ylim(h_lim);

%% Plot profiles (Vertical Eddy Diffusivity)

figure('Name', 'Eddy Diffusivity', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);         plot(rKeff+sKeff, zf, 'k-');  hold on;
plot(rKeff, zf, 'r--'); plot(sKeff, zf, 'b--');
xlabel("K_{eff} (m^2/s)"); ylabel('z (km)');  ylim(h_lim);
legend('Total', 'Resolved', 'Subgrid');
subplot(1,2,2);         plot(rKeff+sKeff, zf, 'k-');  hold on;
xlabel("K_{eff} (m^2/s)"); ylabel('z (km)');  ylim(h_lim);

%%% Theoretical K_eff %%%
h_ref = h_pbl;  kappa = 0.41;  zf_ = zf(zf <= h_ref);
K_theo = 1e3*kappa*ust .* zf_ .* (1-zf_./h_ref).^2;  plot(K_theo, zf_, 'm--');
K_theo = 1e3*kappa*ust .* zf_ .* (1-zf_./h_ref).^4;  plot(K_theo, zf_, 'm-.');
legend({'K', '$$\kappa u_{\ast} z(1-z/h)^2$$', ...
    '$$\kappa u_{\ast} z(1-z/h)^4$$'}, 'Interpreter', 'latex');

%% Plot profiles (Mixing Length)

figure('Name', 'Mixing Length', ...
    'Position', [0, screen(4)/2, screen(3), screen(4)/2]);
subplot(1,2,1);         plot(Leff, zf, 'k-');  hold on;
xlabel("l_{eff} (m)");  ylabel('z (km)');  ylim(h_lim);
subplot(1,2,2);         plot(Leff, zf, 'k-');  hold on;
xlabel("l_{eff} (m)");  ylabel('z (km)');  ylim(h_lim);

%%% Theoretical L_eff %%%
Linf = 75*1e-3;  kappa = 0.41;
L_theo = 1e3*(Linf^(-1) + (kappa.*zf).^(-1)).^(-1/1);  plot(L_theo, zf, 'r--');
L_theo = 1e3*(Linf^(-2) + (kappa.*zf).^(-2)).^(-1/2);  plot(L_theo, zf, 'r-.');
legend({'$l$', '$${l^{-1} = (\kappa z)^{-1} + 75^{-1}}$$', ...
    '$$l^{-2} = (\kappa z)^{-2} + 75^{-2}$$'}, 'Interpreter', 'latex');

%% Save profiles

save('profiles.mat', 'ust', 'h_pbl', ...
    'zf', 'upup', 'vpvp', 'wpwp', 'stke', 'rtke', 'rtau', 'ttau', 'stau', ...
    'rKeff', 'sKeff', 'Leff');
