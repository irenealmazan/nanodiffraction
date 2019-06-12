% this scripts compares the centroids maps, and the strain and tilt maps
% for the rock curves 9, 10 and their combination

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

%% Rock curve 9

filename = 'results/data_scan_zeropad';

lst9 = [45526:1:45546]; 
th_end = 84.426;
th_start = th_end - 1.0;
delta_th = 1.0/numel(lst9);


del = -16.1;
gam = 13.1;%10.65;
twoTheta = 21.2;
detdist = 0.35; % in meters


ROIxstart = 89;%50;
ROIystart = 6;%81;
ROIxsize = 262;%200;
ROIysize = 506;%200;

innerpts = 60; % columns
outerpts = 20; % rows
innerpts_zeropad = 10;
outerpts_zeropad = 10;

dat = load([name_file num2str(lst9(1)) '_' num2str(lst9(end)) '.mat'],'dat1');

dat9 = dat.dat1;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst9)-1);
angs = [th_start:th_step:th_end];

eval('Init_parameters');




filename_toload_910 = 'results/data_scan_alignXRF0';

filename_toload_910 = 'results/data_scan_alignXRF0';
