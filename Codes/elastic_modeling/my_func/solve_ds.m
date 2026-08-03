%% Function: Solve displacement-stress vector (P-SV system)

% Return output vectors at different depths specified by zq
% Output at surface z = 0 is always recorded

function output = solve_ds(src, elast_prop, zq)

    %%% Insert depth query points into layered model %%%
    model_prop = insert_query_points(elast_prop, zq);

    %%% Constant values %%%
    % Number of layers (INCLUDE halfspace)
    Nlayer = size(model_prop, 1);

    % Layers to record output at the top
    irec = (model_prop(:, end) == 1);

    % Calculate Lame parameters
    % Unit: GPa = 1e9 Pa
    mu = model_prop(:, 1) .* model_prop(:, 3).^2;
    lambda = model_prop(:, 1) .* model_prop(:, 2).^2 - 2.*mu;
    sigma = lambda + 2*mu;

    % Layer thickness
    % Unit: km
    h = model_prop(:, 4);

    %%% Mesh grid %%%
    % Spatial axes [km]
    Nx = length(src.xh);  dx = src.xh(2)-src.xh(1);
    Ny = length(src.yh);  dy = src.yh(2)-src.yh(1);

    %%% FFT parameters %%%
    % Wavenumber samples [rad/km]
    kx = 2*pi .* [0:Nx/2 (-Nx/2+1):-1]' ./ (Nx*dx);
    ky = 2*pi .* [0:Ny/2 (-Ny/2+1):-1]  ./ (Ny*dy);

    % Radial wavenumber samples [rad/km]
    Nr = max(Nx, Ny) * 10;  dr = min(dx, dy) / sqrt(3);
    % kr = 2*pi .* [1:Nr/2 (-Nr/2+1):-1]' ./ (Nr*dr);
    kr = 2*pi .* (1:Nr/2)' ./ (Nr*dr);

    %%% Initial homogeneous solution %%%
    lambda0 = lambda(end);  mu0 = mu(end);
    
    ds1 = [1/(2*mu0) ./ kr,     1/(2*mu0) ./ abs(kr), ...
           sign(kr),            ones(size(kr))];
    
    ds2 = [sign(kr).*(lambda0+2*mu0)./ kr.^2 ./(2*mu0*(lambda0+mu0)), ...
           -1./ kr.^2 ./(2*(lambda0+mu0)), ...
           1 ./ kr,     zeros(size(kr))];

    %%% Propagator method %%%
    % Initialize arrays
    ds1_surf = zeros([size(ds1) Nlayer]);
    ds2_surf = zeros([size(ds2) Nlayer]);
    
    % Outer loop over non-zero wavenumber
    parfor j = 1:length(kr)
    
        % Wavenumber for current loop
        kj = kr(j);
    
        % Vector at current wavenumber
        ds1j = ds1(j, :)';  ds2j = ds2(j, :)';
        
        % Inner loop over layers (including halfspace with h(end) = 0)
        for i = Nlayer:-1:1
            
            % ODE system
            Ak = [0, -kj, 1/mu(i), 0; ...
                kj*lambda(i)/sigma(i), 0, 0, 1/sigma(i); ...
                4*mu(i)*(lambda(i)+mu(i))/sigma(i)*kj^2, 0, 0, -kj*lambda(i)/sigma(i); ...
                0, 0, kj, 0];
    
            % Propagator matrix
            ds1j = expm(Ak*h(i)) * ds1j;  ds2j = expm(Ak*h(i)) * ds2j;

            % Record output
            ds1_surf(j,:,i) = ds1j;  ds2_surf(j,:,i) = ds2j;
        end
    end

    %%% Output struct %%%
    output.xh = src.xh;  output.yh = src.yh;  output.dx = dx;  output.dy = dy;
    output.kx = kx;  output.ky = ky;  output.kr = kr;
    output.zq = model_prop(irec, 5);  output.prop = model_prop(irec, 1:3);
    output.ds1 = ds1_surf(:,:,irec);  output.ds2 = ds2_surf(:,:,irec);
end

%% Function: Insert depth query points into layered model

function model_prop = insert_query_points(elast_prop, zq)

    %%% Ensure the last row has zero thickness (halfspace) %%%
    elast_prop(end, 4) = 0;
    
    %%% Combine depth sample points %%%
    % Depth of interfaces [km]
    z_layer = cumsum(elast_prop(1:end-1, 4));

    % Insert depth query points
    z_new = sort(union(zq(zq<z_layer(end) & zq>0), z_layer));


    %%% Initialize new layered model %%%
    model_prop = zeros(length(z_new)+1, 6);

    % 4th column: Layer thickness
    model_prop(2:end-1, 4) = diff(z_new);  model_prop(1, 4) = z_new(1);

    % 5th column: Depth of layer top
    model_prop(2:end, 5) = cumsum(model_prop(1:end-1, 4));

    % 6th (last) column: If to record output at the top
    model_prop(2:end, 6) = ismember(z_new, zq);  model_prop(1, 6) = 1;
    
    % Bottom halfspace
    model_prop(end, 1:4) = elast_prop(end, 1:4);

    % Assign each layer properties
    ind = discretize(z_new, union(0,z_layer), 'IncludedEdge', 'right');
    model_prop(1:end-1, 1:3) = elast_prop(ind, 1:3);

end
