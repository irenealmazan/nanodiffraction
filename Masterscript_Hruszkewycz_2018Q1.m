% this is the master script for the data in Hruszkewycz_2018Q1, CdTe taken
% at HXN (NSLS2) in March 2018


% test alignment

%RockCurve1_script;
RockCurve9_10_script;
%RockCurve11_script;
RockCurve12_script;
RockCurve7_script;


return;

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


%% Rocking curve 1, HXN beamtimes line 135
%{
lst1 = [45193:1:45233]; 
th_start = 81.25;
th_end = 83.75;

ROIxstart = 1;%50;
ROIystart = 133;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

innerpts = 40;
outerpts = 40;
innerpts_zeropad = 45;
outerpts_zeropad = 45;


inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];


 [dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');

%}

%% Rocking curve 2, HXN beamtimes line 289
%{
lst1 = [45325:1:45345]; 
th_start = 82.14;
th_end = 83.14;

ROIxstart = 1;%50;
ROIystart = 1;%81;
ROIxsize = 404;%200;
ROIysize = 404;%200;

innerpts = 30;
outerpts = 30;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];


 [dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 3, HXN beamtimes line 289
%{
lst1_theo = [45346:1:45366];
lst1 = [45346:1:45358]; 
th_start = 82.16;
delta_th = 1.0/numel(lst1);
th_end = th_start*delta_th*numel(lst1);

ROIxstart = 1;%50;
ROIystart = 1;%81;
ROIxsize = 404;%200;
ROIysize = 404%200;

innerpts = 30;
outerpts = 30;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 4, HXN beamtimes line 321 NON COMPLETE
%{
%lst1_theo = [45346:1:45366];
lst1 = [45391:1:45398]; 
th_start = 81.38;
delta_th = 1.0/numel(lst1);
th_end = 82.38;

ROIxstart = 1;%50;
ROIystart = 133;%81;
ROIxsize = 404;%200;
ROIysize = 382;%200;

innerpts = 40;
outerpts = 40;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 5, HXN beamtimes line 323 NON COMPLETE
%lst1_theo = [45346:1:45366];
%{
lst1 = [45400:1:45402 45405:1:45407 45417:1:45420];  
th_start = 81.38;
delta_th = 1.0/numel(lst1);
th_end = 82.38;

ROIxstart = 1;%50;
ROIystart = 133;%81;
ROIxsize = 404;%200;
ROIysize = 382;%200;

innerpts = 40;
outerpts = 40;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 6, HXN beamtimes line 352
%lst1_theo = [45346:1:45366];
%{
lst1 = [45417:1:45437]; 
th_start = 81.38;
delta_th = 1.0/numel(lst1);
th_end = 82.38;

ROIxstart = 1;%50;
ROIystart = 132;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

innerpts = 41;
outerpts = 21;
innerpts_zeropad = 11;
outerpts_zeropad = 11;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');

%}
%% Rocking curve 7, HXN beamtimes line 385 - 390
%lst1_theo = [45346:1:45366];
%{
lst1 = [45477:1:45487 45449:1:45459 45460:1:45470]; 
th_end = 85.08;
th_start = 84.08-0.5;
delta_th = 1.0/numel(lst1);


ROIxstart = 1;%50;
ROIystart = 133;%81;
ROIxsize = 404;%200;
ROIysize = 382;%200;

innerpts = 41;
outerpts = 21;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 8, HXN beamtimes line 446
%lst1_theo = [45346:1:45366];
%{
lst1 = [45504:1:45524]; 
th_end = 84.076;
th_start = th_end - 1.0;
delta_th = 1.0/numel(lst1);


ROIxstart = 1;%50;
ROIystart = 133;%81;
ROIxsize = 404;%200;
ROIysize = 382;%200;

innerpts = 31;
outerpts = 11;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 9, HXN beamtimes line 461
%lst1_theo = [45346:1:45366];
%{
lst1 = [45526:1:45546]; 
th_end = 84.426;
th_start = th_end - 1.0;
delta_th = 1.0/numel(lst1);


ROIxstart = 89;%50;
ROIystart = 6;%81;
ROIxsize = 262;%200;
ROIysize = 506;%200;

innerpts = 60;
outerpts = 20;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 10, HXN beamtimes line 465
%lst1_theo = [45346:1:45366];
%{
lst1 = [45547:1:45567]; 
th_end = 84.45;
th_start = th_end - 1.0;
delta_th = 1.0/numel(lst1);


ROIxstart = 89;%50;
ROIystart = 6;%81;
ROIxsize = 262;%200;
ROIysize = 506;%200;

innerpts = 60;
outerpts = 20;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 11, HXN beamtimes line 471
%{
lst1 = [45569:1:45629]; 
th_end = 84.78;
th_start = th_end - 3.0;
delta_th = 1.0/numel(lst1);


ROIxstart = 0;%50;
ROIystart = 132;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

innerpts = 50;
outerpts = 24;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

%% Rocking curve 12, HXN beamtimes line 495
%{
lst1 = [45640:1:45650]; 
th_end = 82.61;
th_start = th_end - 1.0;
delta_th = 1.0/numel(lst1);


ROIxstart = 6;%50;
ROIystart = 210;%81;
ROIxsize = 295;%200;
ROIysize = 299;%200;

innerpts = 40;
outerpts = 24;
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';

XRFchan = 'Det1_Cu';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

th_step = (th_end-th_start)/(numel(lst1)-1);
angs = [th_start:th_step:th_end];

[dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
 save(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
%}

return;


 
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
skip = 0;

if ~skip
    [dat1,imapx,imapy,sumim]=ND_read_data.ThetaScan_film(datapath,lst1,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',1);
    save(['results/data_scan_zeropad_test' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat'],'dat1','-v7.3');
else
    load(['results/data_scan_zeropad_test' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat']);
end

% Display:

[mask,dat1] = ND_analysis.calculateMask(dat1,0.1);
[rock_curve,thetalist] = ND_display_data.displayRockCurveMaps(lst1,dat1,'do_mask',1,'figNum',3000);
ND_display_data.displayRockCurveShift(lst1,dat1,'figNum',4000);
ND_display_data.displayRockCurveLine(lst1,dat1,[10*ones(40,1),[1:1:40]'],'figNum',5010)
%ND_display_data.displayRockCurveLine(lst1,dat1,[[10:2:30]',19*ones(11,1)],'figNum',5000)
dat1 = ND_analysis.computeCentroids_rockCurve(dat1);
ND_display_data.display2Dmap(dat1.Xcentroids,'figNum',10,'figTitle',['X centroids scans ' num2str(lst1(1)) ' to ' num2str(lst1(end))]);
ND_display_data.display2Dmap(dat1.Ycentroids,'figNum',11,'figTitle',['Y centroids scans ' num2str(lst1(1)) ' to ' num2str(lst1(end))]);

[struct_centroidShift] = ND_analysis.computeCentroidShift(dat1,ROIxstart,ROIxsize,ROIystart,ROIysize,1);

[strain_struct] = ND_analysis.calculateStrain(dat1,ROIxstart,ROIxsize,ROIystart,ROIysize,1);




