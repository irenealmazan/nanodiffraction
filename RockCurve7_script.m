
%% Rocking curve 7, HXN beamtimes line 385 - 390
%%{
lst = [45477:1:45487 45449:1:45459 45460:1:45470]; 
th_middle = 84.576;
th_start = th_middle - 0.75;
th_end = th_middle + 0.75;
delta_th = 1.5/(numel(lst)-1);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns


th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];

del = -18.1;
gam = 10.80;%10.65;
twoTheta = 21.2;
detdist = 0.500; % in meters



ROIxstart = 1;%50;
ROIystart = 132;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

innerpts = 41;
outerpts = 21;
innerpts_zeropad = 40;
outerpts_zeropad = 30;

delta_x_microns = 0.1; % step in the sample frame in microns for supergrid

Grain_label = 'Grain2';

