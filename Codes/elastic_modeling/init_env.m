%% Setting global variables

% Add modeling functions
addpath('./my_func');

% Data directory
data_dir = './outputs/best_case';
cm1_dir = './tmp';

% Get the screen size for plotting
screen = get(0, 'ScreenSize');

% Colorbar: red-white-blue
mcolor = colormap('jet');  mcolor = mcolor(end-16:-1:16, :);

% You may obtain the function `slanCM' from MATLAB file exchange for more colormaps.
% https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap

% You can also just keep it as `mcolor = colormap('jet');'