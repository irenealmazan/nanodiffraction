%% Rocking curve 11, HXN beamtimes line 471
%%{
lst = [45569:1:45629]; 
th_middle = 84.776;
th_start = th_middle - 1.5;
th_end = th_middle + 1.5;
delta_th = 3.0/numel(lst);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns

th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];


del = -17.8;
gam = 11.6;%10.65;
twoTheta = 21.2;
detdist = 0.500; % in meters



ROIxstart = 1;%50;
ROIystart = 132;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

delta_x_microns = 0.1; % step in the sample frame in microns for supergrid


innerpts = 50;
outerpts = 24;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

Grain_label = 'Grain4';
