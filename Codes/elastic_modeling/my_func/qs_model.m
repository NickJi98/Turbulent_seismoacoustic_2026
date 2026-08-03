%% Function: Quasi-static seismic modeling in layered medium (P-SV system)

function [uz, ux, uy] = qs_model(src, elast_prop)

    % Solve displacement-stress vector toward the surface
    ds_pm = solve_ds(src, elast_prop, 0);
    
    % Numerical solution
    sol_pm = calc_disp_layer(src, ds_pm);

    % Return modeling results
    % (For vertical displacement, positive downward)
    uz = -sol_pm.uz;  ux = sol_pm.ux;  uy = sol_pm.uy;

end