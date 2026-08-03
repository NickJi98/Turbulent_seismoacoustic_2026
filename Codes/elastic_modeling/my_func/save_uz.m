%% Function: Save results to NetCDF file

function [] = save_uz(uz, src, filename, is_layer)

    xr = src.xh;  yr = src.yh;  pp = src.pp;
    Nx = length(xr);  Ny = length(yr);  Nt = size(pp, 3);
    data_type = 'single';  deflate_level = 5;

    ncfile = filename;
    
    nccreate(ncfile, 'pp', 'Dimensions', {'x', Nx, 'y', Ny, 't', Nt}, ...
        'Datatype', data_type, 'DeflateLevel', deflate_level);
    ncwrite(ncfile, 'pp', pp);

    if is_layer
        nccreate(ncfile, 'uz_layer', 'Dimensions', {'x', Nx, 'y', Ny, 't', Nt}, ...
            'Datatype', data_type, 'DeflateLevel', deflate_level);
        ncwrite(ncfile, 'uz_layer', uz);
    else
        nccreate(ncfile, 'uz_hs', 'Dimensions', {'x', Nx, 'y', Ny, 't', Nt}, ...
            'Datatype', data_type, 'DeflateLevel', deflate_level);
        ncwrite(ncfile, 'uz_hs', uz);
    end

    nccreate(ncfile, 'xr', 'Dimensions', {'x', Nx}, 'Datatype', data_type);
    ncwrite(ncfile, 'xr', xr);
    
    nccreate(ncfile, 'yr', 'Dimensions', {'y', Ny}, 'Datatype', data_type);
    ncwrite(ncfile, 'yr', yr);

end
    