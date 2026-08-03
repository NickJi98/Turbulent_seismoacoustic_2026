%% Seismic modeling from CM1 output (Layered media or halfspace)

run("init_env.m");

% Loop over each layer for +1% Vs perturbation. 
% If you do not want netcdf outputs, just comment out the final line.

%% Pressure source

%%% Data file %%%
cm1_file = fullfile(cm1_dir, 'cm1out_prs.nc');
fprintf('Read in data: %s\n\n', cm1_file);

%%% Read pressure source %%%
src = read_cm1_output(cm1_file);

%% Elastic property relations

% Brocher (2005), Eq. 9
vp_func = @(vs) polyval([-0.0251, 0.2683, -0.8206, 2.0947, 0.9409], vs);

% Brocher (2005), Eq. 2
rho_func = @(vp) 1.75 .* vp.^0.25;

% Reference model
% Layer properties (rho, vp, vs, thickness)
% Units: g/cm^3, km/s, km/s, km
elast_prop = readmatrix(vel_file, 'NumHeaderLines', 1);
Nlayer = size(elast_prop, 1) - 1;

%% Perturbation of each layer

for i = 1:Nlayer+1

    % Perturbed elastic properties (1%)
    vs_ = elast_prop(i, 3) * 1.01;
    vp_ = vp_func(vs_);  rho_ = rho_func(vp_);

    % Update elastic structure
    elast_prop_ = elast_prop;  elast_prop_(i, 1:3) = [rho_, vp_, vs_];

    % Print layered model
    row_names = cellstr(num2str((1:Nlayer)'))';
    row_names{end+1} = 'Halfspace';
    disp('Layered media:');
    disp(array2table(elast_prop_, ...
        "VariableNames", {'rho (g/cm^3)', 'Vp (km/s)', 'Vs (km/s)', 'Thickness (km)'}, ...
        "RowNames", row_names));

    % Seismic modeling (Layered medium)
    [uz_layer, ~, ~] = qs_model(src, elast_prop_);

    % Save vertical displ. field
    filename_ = fullfile(data_dir, sprintf('LES_kernel_%d.nc', i));
    save_uz(uz_layer, src, filename_, 1);

end
