% this is the master script for the data in Hruszkewycz_2018Q1, CdTe taken
% at HXN (NSLS2) in March 2018

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


% Rocking curve 1, HXN beamtimes line 135
lst1 = [45193:1:45233]; 
th_start = 81.25;
th_end = 83.75;
XRFchan = 7; % Cu fluorescence channel 
innerpts = 40;
outerpts = 40;
innerpts_zeropad = 40;
outerpts_zeropad = 40;
th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

%{
% Rocking curve 2, HXN beamtimes line 208 - not good
% lst1 = [45250:1:45264]; 
% th_start = 82.65;
% th_end = 83.75;

% Rocking curve 3, HXN beamtimes line 295
% lst1 = [45361:1:45363 45365:1:45375]; 
% th_start = 82.176;
% th_end = th_start+0.05*numel(lst1);
% innerpts = 30;
% XRFchan = 7; % Cu fluorescence channel PROBLEM WITH FLUORESCENCE CHANNEL
% 
% 
% th_step = 0.05;%(th_end-th_start)/(numel(lst1)-1);
% angs = [th_start:th_step:th_end];


% Rocking curve 4, HXN beamtimes line 352
% lst1 = [45417:1:45437]; 
% th_start = 82.38;
% th_end = 83.38;
% innerpts = 41;
% outerpts = 21;
% XRFchan = 7; % Cu fluorescence channel 


% th_step = (th_end-th_start)/(numel(lst1)-1);
% angs = [th_start:th_step:th_end];

% test load scan in one single scan
%scanid = 45417;
%[scandata,merlimgs,hfig1,hfig2] = ND_read_data.loadscan_HXN(datapath,scanid,XRFchan,'showmerlin',1,'innerpts',innerpts,'outerpts',outerpts);
%[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'thetalist',angs,'innerpts',40,'showmerlin',0);
%}
% read rocking curve
skip = 1;

if ~skip
    [dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'do_padding',1);
    save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
else
    load(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat']);
end

% Display:
[mask,dat1] = ND_analysis.calculateMask(dat1,0.1);
[rock_curve,thetalist] = ND_display_data.displayRockCurveMaps(lst1,dat1,'do_mask',1,'figNum',3000);
ND_display_data.displayRockCurveShift(lst1,dat1,'figNum',4000);
ND_display_data.displayRockCurveLine(lst1,dat1,[10*ones(40,1),[1:1:40]'],'figNum',5010)
%ND_display_data.displayRockCurveLine(lst1,dat1,[[10:2:30]',19*ones(11,1)],'figNum',5000)
dat1 = ND_analysis.computeCentroids_rockCurve(dat1,1);
ND_display_data.display2Dmap(dat1.Xcentroids,'figNum',10,'figTitle',['X centroids scans ' num2str(lst1(1)) ' to ' num2str(lst1(end))]);
[struct_centroidShift] = ND_analysis.computeCentroidShift(dat1,1);

[strain_struct] = ND_analysis.calculateStrain(dat1,1);



%% Rocking curve 2, HXN beamtimes line 135

