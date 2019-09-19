%% Rocking curve 12, HXN beamtimes line 495
%%{
lst = [45640:1:45650]; 
th_middle = 82.605;
th_end = th_middle + 0.5;
th_start = th_middle - 0.5;
delta_th = 1.0/(numel(lst)-1);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns

th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];

del = -17.0;
gam = 12.50;%10.65;
twoTheta = 21.2;
detdist = 0.500; % in meters

ROIxstart = 6;%50;
ROIystart = 210;%81;
ROIxsize = 295;%200;
ROIysize = 299;%200;


delta_x_microns = 0.1; % step in the sample frame in microns for supergrid



innerpts = 40;
outerpts = 24;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

Grain_label = 'Grain5';
