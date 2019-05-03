% this is the master script for the data in Hruszkewycz_2018Q1, CdTe taken
% at HXN (NSLS2) in March 2018

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';



%Rocking curve 1, HXN beamtimes line 135
lst1 = [45193:1:45233];
th_start = 81.25;
th_end = 83.75;
th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];
XRFchan = 7; % Cu fluorescence channel 

% test load scan in one single scan
%scanid = 45216;
%[scandata,merlimgs,hfig1,hfig2] = ND_read_data.loadscan_HXN(datapath,scanid,XRFchan,'showmerlin',1,'innerpts',40);

% read rocking curve
skip = 1;

if ~skip
    [dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'thetalist',angs,'innerpts',40,'showmerlin',0);
    save('results/data_scan45193_45233.mat','dat1');
else
    load('results/data_scan45193_45233.mat');
end

% Display:
ND_display_data.displayRockCurve(lst1,dat1,'figNum',3000);
ND_display_data.displayRockCurveShift(lst1,dat1,'figNum',4000);

[struct_centroidShift] = ND_analysis.computeCentroidShift(dat1,1);

%[strain_struct] = ND_analysis.calculateStrain(dat1,1);