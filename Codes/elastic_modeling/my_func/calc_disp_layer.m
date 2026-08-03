%% Function: Calculate displacement in layered medium (P-SV system)

function output = calc_disp_layer(src, sol_ds)

    %%% Read input struct %%%
    Nx = length(sol_ds.xh);  Ny = length(sol_ds.yh);
    dx = sol_ds.dx;  dy = sol_ds.dy;
    kx = sol_ds.kx;  ky = sol_ds.ky;  kr = sol_ds.kr;
    ds1_surf = sol_ds.ds1(:,:,1);  ds2_surf = sol_ds.ds2(:,:,1);


    %%% Convert to Cartesian coordinate %%%
    % Wavenumber grid [rad/km]
    [Kx, Ky] = meshgrid(kx, ky);  Kx = Kx';  Ky = Ky';
    Kr = sqrt(Kx.^2 + Ky.^2);     clear Kx Ky;
    
    % Initialize arrays
    ds1_surf_xy = zeros(Nx, Ny, 4);  ds2_surf_xy = zeros(Nx, Ny, 4);
    
    % Interpolation
    parfor j = 1:4
        ds1_surf_xy(:,:,j) = interp1(kr, ds1_surf(:,j), Kr, 'linear', 0);
        ds2_surf_xy(:,:,j) = interp1(kr, ds2_surf(:,j), Kr, 'linear', 0);
    end

    % Fix k = 0 component
    Kr(1, 1) = Inf;


    %%% Solve linear system at each time %%%
    % Number of time steps
    Nt = size(src.pp, 3);

    % Normal traction (Positive for tensile)
    fk_pp = -fft(fft(src.pp,[],1),[],2) .* dx*dy;

    % Matrix of the linear system (pagewise)
    A2 = cat(4, ds1_surf_xy(:,:,3:4), ds2_surf_xy(:,:,3:4));

    % Initialize output arrays
    uz = zeros(Nx,Ny,Nt);  ux = zeros(Nx,Ny,Nt);  uy = zeros(Nx,Ny,Nt);

    for it = 1:Nt

        % Vector of the linear system (pagewise)
        b2 = zeros(1, Nx, Ny, 2);
        b2(1,:,:,:) = cat(3, zeros(Nx, Ny), fk_pp(:,:,it));

        % Solve linear system (pagewise)
        c = pagemldivide(permute(A2, [3,4,1,2]), permute(b2, [4,1,2,3]));
        c = squeeze(permute(c, [3,4,1,2]));

        % Remove NaN (REMOVE k = 0 component)
        c(1,1,:) = 0;

        % Surface displacement
        fk_uz = c(:,:,1) .* ds1_surf_xy(:,:,2) + c(:,:,2) .* ds2_surf_xy(:,:,2);
        fk_ur = c(:,:,1) .* ds1_surf_xy(:,:,1) + c(:,:,2) .* ds2_surf_xy(:,:,1);

        % Solve horizontal displacement
        fk_ux = 1j.*fk_ur.*kx./Kr;  fk_uy = 1j.*fk_ur.*ky./Kr;

        ux(:,:,it) = real(ifft(ifft(fk_ux,[],1),[],2)) ./ (dx*dy);
        uy(:,:,it) = real(ifft(ifft(fk_uy,[],1),[],2)) ./ (dx*dy);
        uz(:,:,it) = real(ifft(ifft(fk_uz,[],1),[],2)) ./ (dx*dy);

    end


    %%% Output struct %%%
    output.x = sol_ds.xh;  output.y = sol_ds.yh;
    output.ux = ux;  output.uy = uy;  output.uz = uz;

end