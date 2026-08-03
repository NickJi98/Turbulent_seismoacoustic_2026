%% Setting global variables

% Data directory
data_dir = '../../../Data';

% Simulation result directory
sim_dir = fullfile(data_dir, 'CM1/outputs');

% Get the screen size for plotting
screen = get(0, 'ScreenSize');

% Colorbar: red-white-blue
mcolor = slanCM('bwr');  mcolor = mcolor(end-16:-1:16, :);

% You may obtain the function `slanCM' from MATLAB file exchange.
% You can also just change it to `mcolor = colormap('jet');'
