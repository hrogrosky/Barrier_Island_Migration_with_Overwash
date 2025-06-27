%% Area
% Calculate cells < 0 , cells = 0, cells > 0
neg_cells = sum(sum(H(:,:,1) < 0));
zero_cells = sum(sum(H(:,:,1) == 0));
pos_cells = sum(sum(H(:,:,1) > 0));

% Volume
% Multiple # of cells > 0 by 4x4x0.1 m^3
volume = pos_cells*4*4*0.1;

format longG
[neg_cells, zero_cells, pos_cells, volume]