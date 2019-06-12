
%%% This script displays the rocking curves generated in
%%% Masterscript_Hruskewycz

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


% Rock curve 1
clear all; close all;
lst1 = [45193:1:45233]; 
dat = load(['results/data_scan_zeropad' num2str(lst1(1)) '_' num2str(lst1(end)) '.mat']);
dat1 = dat.dat1;

[mask,dat1] = ND_analysis.calculateMask(dat1,0.1);
[dat1.rock_curve,dat1.thetalist] = ND_display_data.displayRockCurveMaps(lst1,dat1,'do_mask',1,'figNum',1000);
ND_display_data.displayRockCurveShift(lst1,dat1,'figNum',4000);



% Rocking curve 6, HXN beamtimes line 352

lst6 = [45417:1:45437]; 

dat = load(['results/data_scan_zeropad' num2str(lst6(1)) '_' num2str(lst6(end)) '.mat']);
dat6 = dat.dat1;

[mask,dat6] = ND_analysis.calculateMask(dat6,0.2);
[dat6.rock_curve,dat6.thetalist] = ND_display_data.displayRockCurveMaps(lst6,dat6,'do_mask',1,'figNum',1000);
ND_display_data.displayRockCurveShift(lst6,dat6,'figNum',4000);


% Rocking curve 7, HXN beamtimes line 385 - 390
clear all; close all;

lst7 = [45477:1:45487 45449:1:45459 45460:1:45470]; 

dat = load(['results/data_scan_zeropad' num2str(lst7(1)) '_' num2str(lst7(end)) '.mat']);
dat7 = dat.dat1;

[mask,dat7] = ND_analysis.calculateMask(dat7,0.3);
[dat7.rock_curve,dat7.thetalist] = ND_display_data.displayRockCurveMaps(lst7,dat7,'do_mask',1,'figNum',1000);
ND_display_data.displayRockCurveShift(lst7,dat7,'figNum',4000);


% Rock curve 8 -> no good
clear all; close all;

lst8 = [45504:1:45524]; 

dat = load(['results/data_scan_zeropad' num2str(lst8(1)) '_' num2str(lst8(end)) '.mat']);
dat8 = dat.dat1;
[mask,dat8] = ND_analysis.calculateMask(dat8,0.1);
[dat8.rock_curve,dat8.thetalist] = ND_display_data.displayRockCurveMaps(lst8,dat8,'do_mask',1,'figNum',8000);
ND_display_data.displayRockCurveShift(lst8,dat8,'figNum',8500);


% Rocking curve 9, HXN beamtimes line 461 this the ptycho scan
clear all; close all;

lst9 = [45526:1:45546]; 

dat = load(['results/data_scan_zeropad' num2str(lst9(1)) '_' num2str(lst9(end)) '.mat']);
dat9 = dat.dat1;

[mask,dat9] = ND_analysis.calculateMask(dat9,0.1);
[dat9.rock_curve,dat9.thetalist] = ND_display_data.displayRockCurveMaps(lst9,dat9,'do_mask',1,'figNum',9000);
ND_display_data.displayRockCurveShift(lst9,dat9,'figNum',9500);

% Rocking curve 10, HXN beamtimes line 465
clear all; close all;

lst10 = [45547:1:45567]; 

dat = load(['results/data_scan_zeropad' num2str(lst10(1)) '_' num2str(lst10(end)) '.mat']);
dat10 = dat.dat1;

[mask,dat10] = ND_analysis.calculateMask(dat10,0.1);
[dat10.rock_curve,dat10.thetalist] = ND_display_data.displayRockCurveMaps(lst10,dat10,'do_mask',1,'figNum',10000);
ND_display_data.displayRockCurveShift(lst10,dat10,'figNum',10500);


% Rocking curve 11, HXN beamtimes line 471

clear all; close all;

lst11 = [45569:1:45629]; 

dat = load(['results/data_scan_zeropad' num2str(lst11(1)) '_' num2str(lst11(end)) '.mat']);
dat11 = dat.dat1;

[mask,dat11] = ND_analysis.calculateMask(dat11,0.1);
[dat11.rock_curve,dat11.thetalist] = ND_display_data.displayRockCurveMaps(lst11,dat11,'do_mask',1,'figNum',11000);
ND_display_data.displayRockCurveShift(lst11,dat11,'figNum',11500);


%% Rocking curve 12, HXN beamtimes line 495
%%{
clear all; close all;

lst12 = [45640:1:45650]; 

dat = load(['results/data_scan_zeropad' num2str(lst12(1)) '_' num2str(lst12(end)) '.mat']);
dat12 = dat.dat1;

[mask,dat12] = ND_analysis.calculateMask(dat12,0.1);
[dat12.rock_curve,dat12.thetalist] = ND_display_data.displayRockCurveMaps(lst12,dat12,'do_mask',1,'figNum',12000);
ND_display_data.displayRockCurveShift(lst12,dat12,'figNum',4000);


