% This is the script to make figure 5 where we show for grain 1 the
% following
%(a) Full diffraction intensity map with rings (contours), 
%(b) azimuthal tilt, 
%(c) outofplane-tilt, 
%(c) total tilt, 
%(d) strain/dspacing

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


% what to do?
flag_read_HXN_parameters = 1;
do_all_analysis = 1;


Grain_array = [1];%[1];

flag_read_HXN_parameters = 1;

if do_all_analysis == 1


eval(['RockCurve' num2str(Grain_array(1)) '_script']);

load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

[thetalist,rock_curve] = ND_display_data.displayRockCurve(dat,'figNum',900);

[~,dat] = ND_data_processing.calculateMask(dat,0.01);

if Grain_array(1) == 9
dat.mask(5:18,20:65) = 0;
end

[centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',dat.mask,'do_plot',0);
dat.Xcentroids = centroid_struct.Xcentroids;
dat.Ycentroids = centroid_struct.Ycentroids;
dat.Thcentroids = centroid_struct.Thcentroids;

[struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(dat,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);

contour_values_up = [2:1:100]*0.01;%[2:1:90]*0.01;%[2:1:6 7:2:22 22:5:30 35:10:90]*0.01; %grains 1, 7 and 9%[7:2:20 22:5:30 35:10:90]*0.01; %Grain 11 and 12% 
contour_values_down = [1:1:99]*0.01;%[1:1:89]*0.01;%[1:1:5 6:2:20 20:5:25 30:10:85]*0.01;% grains 1, 7 and 9% [6:2:18 20:5:25 30:10:85]*0.01; %Grain 11 and 12%
[distr_struct_dspace,mask_struct] = ND_analysis.computeStrainOrTiltContours(dat,struct_centroidShift,'dspace',contour_values_up,contour_values_down);
[distr_struct_tilt_tot,~] = ND_analysis.computeStrainOrTiltContours(dat,struct_centroidShift,'tilt_tot',contour_values_up,contour_values_down);



grain_struct.map2D_SumInt = dat.map2D_SumInt;
grain_struct.mask0 = dat.mask;
grain_struct.mask_rock = dat.mask_rock;
grain_struct.Xval = dat.xsuperGrid_s(1,:);
grain_struct.Yval = dat.ysuperGrid_s(:,1);
grain_struct.Xcentroids = centroid_struct.Xcentroids;
grain_struct.Ycentroids = centroid_struct.Ycentroids;
grain_struct.Thcentroids = centroid_struct.Thcentroids;
grain_struct.strain = struct_centroidShift.strain;
grain_struct.dspace = struct_centroidShift.dspace;
grain_struct.tilt_x = struct_centroidShift.tilt_x;
grain_struct.tilt_y = struct_centroidShift.tilt_y;
grain_struct.tilt_tot = struct_centroidShift.tilt_tot;
grain_struct.mask_struct = mask_struct;
grain_struct.sigma_dspace_distr = distr_struct_dspace.sigma;
grain_struct.dspace_distr = distr_struct_dspace.distr;
grain_struct.sigma_tilt_tot_distr = distr_struct_tilt_tot.sigma;
grain_struct.tilt_tot_distr = distr_struct_tilt_tot.distr;
grain_struct.contour_values_up = contour_values_up;
grain_struct.contour_values_down = contour_values_down;
grain_struct.contours_to_plot = [2 5 20 50];% grains 1, 7 and 9%[3 5 10 50];% Grain 11 and 12%[2 5 7 12 17]; grains 1, 7 and 9
grain_struct.color_array = ['r','y','g','k'];%['r','y','g','m','k']; % grains 1, 7 and 9
grain_struct.ylim_for_dspace = [3.75 3.78];%Grain 11%[3.74 3.8]; %Grain 9%[3.7 3.78]; %Grain 1%;%
grain_struct.ylim_for_tilt_tot = [0 0.2]; %Grain 1%[0 0.5]; Grain 9
grain_struct.ylim_for_strain = [-1 1]*3e-3;%[3.78 3.86]; Grain 1
grain_struct.ylim_for_tilt_x = [-.4 .6];%[3.78 3.86]; Grain 1
grain_struct.ylim_for_tilt_y = [-.4 .6];%[3.78 3.86]; Grain 1
grain_struct.thetalist = thetalist;
grain_struct.rock_curve = rock_curve;

grain_struct.window_for_maps = [1,size(dat.xsuperGrid_s,1),1,size(dat.xsuperGrid_s,2)];
grain_struct.extra_shift = [0 0];

else
   load(['./results_paper/Figure5_grain' num2str(Grain_array(1)) '_struct.mat'],'grain_struct');
  
end

fig_num = 100;
hfig = ND_paper_figures.display2DmapContoursfig5(grain_struct,{'map2D_SumInt','tilt_x','tilt_y','tilt_tot','dspace'},'extra_shift',grain_struct.extra_shift,'window',grain_struct.window_for_maps,'contours_to_plot',grain_struct.contours_to_plot,'color_array',grain_struct.color_array,'size_figure',[6 60 472 691],'figNum',fig_num);
ND_paper_figures.display2DmapHistogramfig5(grain_struct,{'dspace_distr','sigma_dspace_distr'},'ylabel','D [Angstroms]','ylim',grain_struct.ylim_for_dspace,'figNum',1);
ND_paper_figures.display2DmapHistogramfig5(grain_struct,{'tilt_tot_distr','sigma_tilt_tot_distr'},'ylabel','Total tilt [degrees]','ylim',grain_struct.ylim_for_tilt_tot,'figNum',2);

% histograms:
numbin = 30;

figNum = 10;

grain_struct.numbins_dspace = [grain_struct.ylim_for_dspace(1):(grain_struct.ylim_for_dspace(2)-grain_struct.ylim_for_dspace(1))/numbin:grain_struct.ylim_for_dspace(2)];
[hstruct_dspace,sigma_dspace] = ND_paper_figures.displayHistogramRings(grain_struct,distr_struct_dspace,'numBins',grain_struct.numbins_dspace,'xlim',grain_struct.ylim_for_dspace,'figTitle','D','size_figure',[6 60 472 691],'figNum',figNum);

figNum = 9;
grain_struct.numbins_tilt_tot = [grain_struct.ylim_for_tilt_tot(1):(grain_struct.ylim_for_tilt_tot(2)-grain_struct.ylim_for_tilt_tot(1))/numbin:grain_struct.ylim_for_tilt_tot(2)];
[hstruct_tilt_tot,sigma_tilt_tot] = ND_paper_figures.displayHistogramRings(grain_struct,distr_struct_tilt_tot,'numBins',grain_struct.numbins_tilt_tot,'xlim',grain_struct.ylim_for_tilt_tot,'figTitle','Tilt','size_figure',[6 60 472 691],'figNum',figNum);
 
grain_struct.hstruct_dspace = hstruct_dspace;
grain_struct.hstruct_tilt_tot = hstruct_tilt_tot;
 
save(['./results_paper/Figure5_grain' num2str(Grain_array(1)) '_struct.mat'],'grain_struct','-v7.3');

savefig(figure(fig_num),['./results_paper/Figure5_grain' num2str(Grain_array(1)) '.fig']);
savefig(figure(1),['./results_paper/Figure5_grain' num2str(Grain_array(1)) 'dspace_hist.fig']);
savefig(figure(2),['./results_paper/Figure5_grain' num2str(Grain_array(1)) 'tilt_hist.fig']);

savefig(figure(10),['./results_paper/Figure5_grain' num2str(Grain_array(1)) 'dspace_hist.fig']);
savefig(figure(11),['./results_paper/Figure5_grain' num2str(Grain_array(1)) 'tilt_hist.fig']);
  
 
