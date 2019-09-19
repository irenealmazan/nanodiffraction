% This scripts prepares the Work diagram in Figure 4: 
% summed CCD image with maping of the azhimutal and the radial components
% centroids in (x,y) and (theta,azhimut) components 
% all this in grain 1.
% It also plots maps such as the integrated diffracted intensity, the
% ThCentroid, tilt_x, tilt_y and dspace that are displayed for completitute but which
% are not stored as figures in Figure4/
% It also calculates the error in the centroids positions introduced by the
% fringes and which is discussed in the answer to the referees



clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

folder_name = 'results_paper/Figure4/';

filename = 'results/data_scan';
filename_toload = 'results/data_scan_onlyread_nopad';


% what to do?
flag_read_HXN_parameters = 1;

Grain_array = [1];%[1];


load(['results_paper/grain' num2str(Grain_array(1)) '_struct.mat']);

%mask = percent of the diffracted intensity to show!!
percent = 0.05;
mask = grain_struct.map2D_SumInt > percent*max(max(grain_struct.map2D_SumInt));
grain_struct.mask0 = mask;

if Grain_array == 11
    mask(22:29,437:550) = 0;
end

[mask_nan] = ND_data_processing.turnMaskInNan(double(mask));

grain_struct.mask_nan = mask_nan;

figure_size = [-97 33 1378 672];

grain_struct.ylim_for_radial_proj = [min(min(grain_struct.radial_proj.*grain_struct.mask_nan)) max(max(grain_struct.radial_proj.*grain_struct.mask_nan))];%[20.5 20.8];

fields_to_plot = {'radial_proj'};
figNum = 60;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,fields_to_plot,...
'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);

%%{
namefig = [folder_name 'Radial_map_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}

fields_to_plot = {'azi_proj'};

grain_struct.ylim_for_azi_proj = [min(min(grain_struct.azi_proj.*grain_struct.mask_nan)) max(max(grain_struct.azi_proj.*grain_struct.mask_nan))];%[37.0 38];%[20.5 20.8];

figNum = 61;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,fields_to_plot,...
'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);

%%{
namefig = [folder_name 'Azhimutal_map_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}

%{

ND_display_data.display2Dmap(grain_struct.map2D_SumInt_shift,'figNum',6,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Integrated diffracted intensity','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

savefig(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpng');
%}

fields_to_plot = {'Xcentroid'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [min(min(grain_struct.Xcentroid.*grain_struct.mask_nan)) max(max(grain_struct.Xcentroid.*grain_struct.mask_nan))];%[190 210];
figNum = 7;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);

%%{
namefig = [folder_name 'Xcentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}


%{
ND_display_data.display2Dmap(grain_struct.Xcentroid_shift.*grain_struct.mask_nan_shift,'figNum',7,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Column Centroids','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([130 220]);

savefig(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1))],'-dpng');
%}

fields_to_plot = {'Ycentroid'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) =[min(min(grain_struct.Ycentroid.*grain_struct.mask_nan)) max(max(grain_struct.Ycentroid.*grain_struct.mask_nan))];%[180 210];
figNum = 8;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);

%%{
namefig = [folder_name 'Ycentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}

%{
ND_display_data.display2Dmap(grain_struct.Ycentroid_shift.*grain_struct.mask_nan_shift,'figNum',8,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Rows Centroids ','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
%caxis([170 210]);
caxis([130 220]);

savefig(figure(8),['results_paper/Figure3/Ycentroid_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(8),['results_paper/Figure3/Ycentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(8),['results_paper/Figure3/Ycentroids_grain ' num2str(Grain_array(1))],'-dpng');

%}



fields_to_plot = {'Thcentroid'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [min(min(grain_struct.Thcentroid.*grain_struct.mask_nan)) max(max(grain_struct.Thcentroid.*grain_struct.mask_nan))];%[82.5 82.9];
figNum = 9;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);


%%{
namefig = [folder_name 'Thcentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}

    
    %%%%% integrated CCD

flag_read_HXN_parameters = 1;

eval(['RockCurve' num2str(Grain_array(1)) '_script']);

load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

%dat_supergrid = dat;
 

%  ccd or each angle:
%pixel = [25,26];%-extra_shift;
%pixel = [36,227];%-extra_shift;
pixel = [35,197];%

XRF_maps_to_show = [21];%[2:10:numel(dat.scan)];

counter = 1;
clear dat_supergrid_toshow;
for kk = XRF_maps_to_show
    dat_supergrid_toshow.scan(counter) = dat.scan(kk);
    counter = counter + 1;
end
 dat_supergrid_toshow.ii = dat.ii;

 for kkkk = 1:numel(dat.scan)
 dat_supergrid_toshow.thetalist(kkkk) = dat.scan(kkkk).theta;
 end

 dat_supergrid_toshow.Xcentroids = grain_struct.Xcentroid;
 dat_supergrid_toshow.Ycentroids = grain_struct.Ycentroid;

 figNum = 20;
[figNum_final] = ND_paper_figures.displayCCD(dat_supergrid_toshow,pixel,ROIystart,ROIxstart,ROIysize,ROIxsize,figNum);

rock_curve_pixel = ND_paper_figures.displayCCDsingletheta(dat_supergrid_toshow,pixel,4);


namefig = [folder_name 'sum_CCD' num2str(Grain_array(1)) '_pixel' num2str(pixel(1)) 'and' num2str(pixel(2))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%%%%%% error in the CoM position introduced by the fringes

% integrated diffracted intensity contained in the Bragg reflection
test_im = dat_supergrid_toshow.ii(pixel(1)).jj(pixel(2)).im;
mask_centroid = test_im > 5e-3*max(max(test_im)); % corresponds to a percentage of 0.5
integrated_int_Bragg = sum(sum(test_im.*mask_centroid));

% total integrated diffracted intensity in the summed CCD
total_int = sum(sum(test_im));

% percentage:

int_Bragg_percent = 100*integrated_int_Bragg/total_int;

[tempycen, tempxcen] = ND_analysis.computeCentroids_single_CCD(dat_supergrid_toshow,pixel,'mask',mask_centroid);



% do the same for the full CCD (including fringes)
test_im = dat_supergrid_toshow.ii(pixel(1)).jj(pixel(2)).im;
mask_centroid_2 = test_im > 1e-4*max(max(test_im)); % corresponds to a percentage of 0.5
integrated_int_Bragg = sum(sum(test_im.*mask_centroid_2));

% total integrated diffracted intensity in the summed CCD
total_int = sum(sum(test_im));

% percentage:

int_Bragg_percent = 100*integrated_int_Bragg/total_int;

[tempycen_2, tempxcen_2] = ND_analysis.computeCentroids_single_CCD(dat_supergrid_toshow,pixel,'mask',mask_centroid_2);
tempycen_error = abs(tempycen_2-tempycen)/tempycen_2;
tempxcen_error = abs(tempxcen_2-tempxcen)/tempxcen_2;

figure(900);
clf;
subplot(121);
imagesc(log10(test_im.*mask_centroid));
colorbar;
axis image;
title({['x_0 = ' num2str(tempxcen)],[ 'y_0 = ' num2str(tempycen)]});

%figure(900);
subplot(122);
imagesc(log10(test_im.*mask_centroid_2));
colorbar;
axis image;
title({['x_0 = ' num2str(tempxcen_2)],[ 'y_0 = ' num2str(tempycen_2)]});



figNum = 900;
namefig = [folder_name 'error_sum_CCD' num2str(Grain_array(1)) '_pixel' num2str(pixel(1)) 'and' num2str(pixel(2))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
    
    
%%%%%%%%%%%% extra plots

   % plot of the integrated diffracted intensity (skip because it doesn't
% belong to figure 4

fields_to_plot = {'map2D_SumInt'};
figNum = 6;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,fields_to_plot,...
'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',5,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);


%{
namefig = [folder_name 'Integrated_diff_int_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}
 



%{
ND_paper_figures.displayCCDrowsColsDistribution(dat_supergrid,pixel,101)
savefig(figure(101),['results_paper/Figure3/CCD_rows_distrib_grain ' num2str(Grain_array(1)) '.fig']);
savefig(figure(102),['results_paper/Figure3/CCD_cols_distrib_grain ' num2str(Grain_array(1)) '.fig']);
%}


% strain and tilt:

fields_to_plot = {'dspace'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [3.8 3.86];
figNum = 11;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);


%{
namefig = [folder_name 'dspace_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}


fields_to_plot = {'tilt_x'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [-0.15 .2];
figNum = 13;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);


%{
namefig = [folder_name 'Tiltx_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}


fields_to_plot = {'tilt_y'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [-0.15 .2];
figNum = 14;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',figure_size,...
    'figNum',figNum,'font',30);

set(gcf,'Position',figure_size);


%{
namefig = [folder_name 'Tilty_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');
%}
