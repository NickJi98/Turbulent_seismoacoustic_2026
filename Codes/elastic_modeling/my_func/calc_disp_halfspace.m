%% Function: Calculate displacement in elastic halfspace

function [uz, ux, uy] = calc_disp_halfspace(src, elast_prop, apply_window)

    %%% Constant values
    % Elastic properties
    rho = elast_prop(1);  vp = elast_prop(2);  vs = elast_prop(3);

    % Shear modulus [GPa] & Poisson's ratio
    mu = rho * vs^2;  nu = ((vp/vs)^2-2) / ((vp/vs)^2-1) / 2;

    % Constants for static Green's function [1/GPa]
    const_z = 2*(1-nu)/(4*pi*mu);  const_xy = -(1-2*nu)/(4*pi*mu);

    %%% Mesh grid
    % Spatial axes [km]
    Nx = length(src.xh);  dx = src.xh(2)-src.xh(1);
    Ny = length(src.yh);  dy = src.yh(2)-src.yh(1);
    dS = dx * dy;

    %%% FFT parameters
    % Wavenumber samples [1/km]
    kx = [0:Nx/2 (-Nx/2+1):-1] ./ (Nx*dx);
    ky = [0:Ny/2 (-Ny/2+1):-1] ./ (Ny*dy);

    %%% Window function
    if apply_window

        % Hanning window
        % Hann = @(i,N) 0.5 .* (1-cos(2*pi/N .* (i-0.5)));
        % Wx = Hann(1:Nx,Nx)';  Wy = Hann(1:Ny,Ny);
        % Wt = Hann(reshape(1:Nt,[1,1,Nt]),Nt);

        % Cosine window
        edge = 0.05;
        Wx = tukeywin(Nx, edge*2);  Wy = tukeywin(Ny, edge*2)';
        
        % Apply window function
        Cx = Nx / sum(Wx.^2);  Cy = Ny / sum(Wy.^2);
        src.pp = src.pp .* Wx .* Wy .* sqrt(Cx*Cy);
        disp('Apply window function.');
    end

    %%% FK static Green's function
    % Vertical
    Rz = 1 ./ sqrt(kx'.^2 + ky.^2);
    % Horizontal
    Rx = -1j .* kx' .* Rz.^2;  Ry = -1j .* ky .* Rz.^2;
    % Singularity
    dSk = 1 / (Nx*Ny*dS);  Rz(1, 1) = 2 / sqrt(dSk/pi) / (dx*dy);
    Rx(1, :) = 0;          Ry(1, :) = 0;

    %%% Convolution with source
    % Pressure perturbation
    fk_pp = fft(fft(src.pp,[],1),[],2);  fk_pp(1, 1) = 0;
    % Vertical displacement [μm]
    fk_uz = fk_pp .* Rz;
    uz = real(ifft(ifft(fk_uz,[],1),[],2)) .* const_z;
    % Horizontal displacement [μm]
    fk_ux = fk_pp .* Rx;  fk_uy = fk_pp .* Ry;
    ux = real(ifft(ifft(fk_ux,[],1),[],2)) .* const_xy;
    uy = real(ifft(ifft(fk_uy,[],1),[],2)) .* const_xy;

end