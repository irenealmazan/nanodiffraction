
%% Rocking curve 10, HXN beamtimes line 465
%lst1_theo = [45346:1:45366];
%%{
lst = [45547:1:45567]; 
th_middle = 84.452;
th_start = th_middle - 0.5;
th_end = th_middle + 0.5;
delta_th = 1.0/(numel(lst)-1);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns


th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];

del = -16.1;
gam = 13.1;%10.65;
twoTheta = 21.2;
detdist = 0.35; % in meters

ROIxstart = 89;%50;
ROIystart = 6;%81;
ROIxsize = 262;%200;
ROIysize = 506;%200;


innerpts = 60;
outerpts = 20;
innerpts_zeropad = 10;
outerpts_zeropad = 10;


delta_x_microns = 0.1; % step in the sample frame in microns for supergrid
