% This scripts prepares the Work diagram: 
%stretching 
%alignment 
% centroids 
%strain and tilt explanations
% all this in grain 1



clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


% what to do?
flag_read_HXN_parameters = 1;

Grain_array = [1];%[1];

flag_read_HXN_parameters = 1;

eval(['RockCurve' num2str(Grain_array(1)) '_script']);

load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

dat_supergrid = dat;

%{
load([filename_toload num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');    

[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('x_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('y_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);

dat_orig = dat;

load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

dat_supergrid = dat;

%%%% Display:

XRF_maps_to_show = [13 29];%[2:10:numel(dat.scan)];


counter = 1;
%clear XRF_aligned_toshow


for kk = XRF_maps_to_show

    XRF_aligned_toshow.scan(counter).XRF = dat_orig.scan(kk).XRF;
    XRF_aligned_toshow.scan(counter).theta = dat_orig.scan(kk).theta;
    XRF_aligned_toshow.scan(counter).z_pos = dat_orig.imapx; 
    XRF_aligned_toshow.scan(counter).y_pos = dat_orig.imapy;
    XRF_aligned_toshow.scan(counter).s_pos = dat_orig.imapx.*tand(dat_orig.scan(kk).theta);
    counter = counter + 1;
end

ND_paper_figures.displayAlignedmaps(XRF_aligned_toshow,'XRF',...
    'figNum',1,'calim',[min(min(XRF_aligned_toshow.scan(1).XRF)) max(max(XRF_aligned_toshow.scan(1).XRF))],...
    'size_figure',[1 52 588 653],...
    'font',20,'numRows',numel(XRF_maps_to_show));

namefig = ['results_paper/Figure3/original_XRF_grain' num2str(Grain_array(1))];
savefig(figure(1),[namefig  '.fig']);
print(figure(1),namefig,'-dpdf');
print('-r600',figure(1),namefig,'-dpng');

counter = 1;
clear XRF_aligned_toshow
for kk = XRF_maps_to_show
   
    XRF_aligned_toshow.scan(counter).XRF = circshift(dat_supergrid.scan(kk).XRF,[0 100]);
    XRF_aligned_toshow.scan(counter).s_pos = dat_supergrid.xsuperGrid_s(1,:);
    XRF_aligned_toshow.scan(counter).y_pos = dat_supergrid.ysuperGrid_s(:,1);
    XRF_aligned_toshow.scan(counter).z_pos = dat_supergrid.xsuperGrid_s(1,:)./tand(dat_supergrid.scan(kk).theta);
    XRF_aligned_toshow.scan(counter).theta = dat_supergrid.scan(kk).theta;

    counter = counter + 1;
end

ND_paper_figures.displayAlignedmaps(XRF_aligned_toshow,'XRF','figNum',2,'size_figure',[6 70 473 481],'numRows',numel(XRF_maps_to_show));

namefig = ['results_paper/Figure3/supergrid_XRF_grain' num2str(Grain_array(1))];
savefig(figure(2),[namefig  '.fig']);
print(figure(2),namefig,'-dpdf');
print('-r600',figure(2),namefig,'-dpng');



counter = 1;
clear dat_supergrid_toshow;
dat_supergrid_toshow.ii = dat_supergrid.ii;


for kk = XRF_maps_to_show
    
    %dat_supergrid_toshow.scan(counter) = dat_supergrid.scan(kk);
    dat_supergrid_toshow.scan(counter).s_pos = dat_supergrid.xsuperGrid_s(1,:);
    dat_supergrid_toshow.scan(counter).y_pos = dat_supergrid.ysuperGrid_s(:,1);
    dat_supergrid_toshow.scan(counter).z_pos = dat_supergrid.xsuperGrid_s(1,:)./tand(dat_supergrid.scan(kk).theta);
    dat_supergrid_toshow.scan(counter).theta = dat_supergrid.scan(kk).theta;

    counter = counter + 1;
end


% diffraction maps for each theta value
figNum = 3;
[thetalist,dat_supergrid_toshow] = ND_paper_figures.displayRockCurveMaps(dat_supergrid_toshow,'extra_shift',[0 0],'size_figure',[6 70 473 481],'figNum',figNum);

namefig = ['results_paper/Figure3/diffraction_maps_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf');
print('-r600',figure(figNum),namefig,'-dpng');

%{

savefig(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1))],'-dpng');
%}

%  ccd or each angle:
pixel = [25,26];%-extra_shift;

counter = 1;
clear dat_supergrid_toshow;
for kk = XRF_maps_to_show
    dat_supergrid_toshow.scan(counter) = dat_supergrid.scan(kk);
    counter = counter + 1;
end
 dat_supergrid_toshow.ii = dat_supergrid.ii;

 for kkkk = 1:numel(dat_supergrid.scan)
 dat_supergrid_toshow.thetalist(kkkk) = dat_supergrid.scan(kkkk).theta;
 end

rock_curve_pixel = ND_paper_figures. displayCCDsingletheta(dat_supergrid_toshow,pixel,4);

figNum = 4;
namefig = ['results_paper/Figure3/ccd_image1_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf');
print('-r600',figure(figNum),namefig,'-dpng');

figNum = 5;
namefig = ['results_paper/Figure3/ccd_image2_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf');
print('-r600',figure(figNum),namefig,'-dpng');


%{
savefig(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1))],'-dpng');


savefig(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1))],'-dpng');
%}
varnames = {['theta'] ['intensity']};
filename = ['results_paper/Figure3/rock_curve_single_pixel_grain' num2str(Grain_array(1)) '.xlsx'];
T = table(dat_supergrid_toshow.thetalist',rock_curve_pixel','VariableNames',varnames);
writetable(T,filename,'Sheet',1,'Range','D1');
%}

% integrated diffracted intensity map and centroid maps, shifted in order
% to avoid the splitting in the figure

[mask0,dat_supergrid] = ND_data_processing.calculateMask(dat_supergrid,0.05);

[mask_nan] = ND_data_processing.turnMaskInNan(double(dat_supergrid.mask));

[centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat_supergrid,'mask',dat_supergrid.mask,'do_plot',0);
dat_supergrid.Xcentroids = centroid_struct.Xcentroids;
dat_supergrid.Ycentroids = centroid_struct.Ycentroids;
dat_supergrid.Thcentroids = centroid_struct.Thcentroids;
    
[struct_centroidShift,angles] = ND_analysis.computeCentroidShiftAndStrain(dat_supergrid,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],1);

%%%% radial and azhimutal maps:

figNum = 11;
figure(figNum);axis image;
namefig = ['results_paper/Figure3/Radial' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


figNum = 12;
figure(figNum);axis image;
namefig = ['results_paper/Figure3/Azimuthal' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');




grain_struct.extra_shift = [0 100];
grain_struct.map2D_SumInt = dat_supergrid.map2D_SumInt;%circshift(dat_supergrid.map2D_SumInt,grain_struct.extra_shift);
grain_struct.Xcentroid = dat_supergrid.Xcentroids; %circshift(dat_supergrid.Xcentroids,grain_struct.extra_shift);
grain_struct.Ycentroid = dat_supergrid.Ycentroids;%circshift(dat_supergrid.Ycentroids,grain_struct.extra_shift);
grain_struct.Thcentroid = dat_supergrid.Thcentroids;%circshift(dat_supergrid.Thcentroids,grain_struct.extra_shift);
grain_struct.mask0 = mask0;%circshift(mask0,grain_struct.extra_shift);
grain_struct.mask_nan = mask_nan;%circshift(mask_nan,grain_struct.extra_shift);
grain_struct.angles = angles;
grain_struct.xsuperGrid_s = dat_supergrid.xsuperGrid_s;
grain_struct.ysuperGrid_s = dat_supergrid.ysuperGrid_s;
grain_struct.dspace = struct_centroidShift.dspace;%circshift(struct_centroidShift.dspace,grain_struct.extra_shift);
grain_struct.tilt_x = struct_centroidShift.tilt_x;%circshift(struct_centroidShift.tilt_x,grain_struct.extra_shift); 
grain_struct.tilt_y = struct_centroidShift.tilt_y;%circshift(struct_centroidShift.tilt_y,grain_struct.extra_shift); 
grain_struct.radial_proj = struct_centroidShift.radial_proj;%circshift(struct_centroidShift.tilt_y,grain_struct.extra_shift); 
grain_struct.window_for_maps = [1,50,1,446];
grain_struct.contour_values_down = [1:89]*0.01;
grain_struct.Xval = dat_supergrid.xsuperGrid_s(1,:);
grain_struct.Yval = dat_supergrid.ysuperGrid_s(:,1);
grain_struct.theta_last = dat_supergrid.scan(end).theta;

grain_struct.Yval_range = [min(grain_struct.Yval) max(grain_struct.Yval)];
grain_struct.Yval_step =  (max(grain_struct.Yval)-min(grain_struct.Yval))/numel(grain_struct.Yval);
grain_struct.Yval = ([1:numel(grain_struct.Yval)]-(numel(grain_struct.Yval)/2))*grain_struct.Yval_step;

grain_struct.ylim_for_radial_proj = [20.5 20.8];

fields_to_plot = {'radial_proj'};
figNum = 60;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,fields_to_plot,...
'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Radial_map_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');



fields_to_plot = {'map2D_SumInt'};
figNum = 6;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,fields_to_plot,...
'Yval_lim',[-2 2],'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',5,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Integrated_diff_int_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%{

ND_display_data.display2Dmap(grain_struct.map2D_SumInt_shift,'figNum',6,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Integrated diffracted intensity','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

savefig(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpng');
%}

fields_to_plot = {'Xcentroid'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [170 225];
figNum = 7;

figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Xcentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');



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
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [180 210];
figNum = 8;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Ycentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');

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
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [82.5 82.9];
figNum = 9;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'Yval_lim',[-2 2],...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Thcentroids_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%{
ND_display_data.display2Dmap(grain_struct.Thcentroid_shift.*grain_struct.mask_nan_shift,'figNum',9,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Theta Centroids ','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([82.0 83.60]);

savefig(figure(9),['results_paper/Figure3/Thcentroid_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(9),['results_paper/Figure3/Thcentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(9),['results_paper/Figure3/Thcentroids_grain ' num2str(Grain_array(1))],'-dpng');
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
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/dspace_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%{
ND_display_data.display2Dmap(grain_struct.dspace_shift.*grain_struct.mask_nan_shift,'figNum',11,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Lattice spacing','font',30,'size_figure',[-8 334 785 367]);caxis([200 215]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
%caxis([130 220]);
caxis([3.78 3.86]);

savefig(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1))],'-dpng');
%}

%figNum = figNum_final+2;



fields_to_plot = {'tilt_x'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [-0.15 .2];
figNum = 13;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Tiltx_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');

%{
ND_display_data.display2Dmap(grain_struct.tilt_x_shift.*grain_struct.mask_nan_shift,'figNum',13,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Azhimutal tilt','font',30,'size_figure',[-8 334 785 367]);
hold on;contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([-1.0 .6]);

savefig(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1))],'-dpng');

%}



fields_to_plot = {'tilt_y'};
grain_struct.(['ylim_for_' fields_to_plot{1}]) = [-0.15 .2];
figNum = 13;


figNum = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
    fields_to_plot,'extra_shift',grain_struct.extra_shift,...
    'window',grain_struct.window_for_maps,'plot_contours',0,...
    'min_contour',5,'spec_ylim',0,'size_figure',[-97 33 1378 672],...
    'figNum',figNum,'font',30);


namefig = ['results_paper/Figure3/Tilty_grain' num2str(Grain_array(1))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%{

ND_display_data.display2Dmap(grain_struct.tilt_y_shift.*grain_struct.mask_nan_shift,'figNum',14,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Out-of-plane tilt','font',30,'size_figure',[-8 334 785 367]);
hold on;contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([-1.0 .6]);

savefig(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1))],'-dpng');

%}

%{
figNum = 100;
figure(figNum);
subplot(211);
imagesc(dat_supergrid.xsuperGrid_s);
colorbar;

subplot(212);
imagesc(dat_supergrid.ysuperGrid_s);
colorbar;

figNum = 101;
figure(figNum);
subplot(211);
imagesc(dat_orig.imapx);
colorbar;

subplot(212);
imagesc(dat_supergrid.imapy);
colorbar;
%}


%%%%% integrated CCD

%  ccd or each angle:
%pixel = [25,26];%-extra_shift;
%pixel = [36,227];%-extra_shift;
pixel = [35,197];%

XRF_maps_to_show = [13 29];%[2:10:numel(dat.scan)];

counter = 1;
clear dat_supergrid_toshow;
for kk = XRF_maps_to_show
    dat_supergrid_toshow.scan(counter) = dat_supergrid.scan(kk);
    counter = counter + 1;
end
 dat_supergrid_toshow.ii = dat_supergrid.ii;

 for kkkk = 1:numel(dat_supergrid.scan)
 dat_supergrid_toshow.thetalist(kkkk) = dat_supergrid.scan(kkkk).theta;
 end

 dat_supergrid_toshow.Xcentroids = grain_struct.Xcentroid;
 dat_supergrid_toshow.Ycentroids = grain_struct.Ycentroid;

 figNum = 20;
[figNum_final] = ND_paper_figures.displayCCD(dat_supergrid_toshow,pixel,angles,ROIystart,ROIxstart,ROIysize,ROIxsize,figNum);

rock_curve_pixel = ND_paper_figures. displayCCDsingletheta(dat_supergrid_toshow,pixel,4);


namefig = ['results_paper/Figure3/sum_CCD' num2str(Grain_array(1)) '_pixel' num2str(pixel(1)) 'and' num2str(pixel(2))];
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

figure(900);
clf;
subplot(121);
imagesc(log10(test_im.*mask_centroid));
colorbar;
axis image;
title(['x\_centroid = ' num2str(tempxcen) 'y\_centroid = ' num2str(tempycen)]);

% do the same for the ful CCD (including fringes)
test_im = dat_supergrid_toshow.ii(pixel(1)).jj(pixel(2)).im;
mask_centroid_2 = test_im > 1e-4*max(max(test_im)); % corresponds to a percentage of 0.5
integrated_int_Bragg = sum(sum(test_im.*mask));

% total integrated diffracted intensity in the summed CCD
total_int = sum(sum(test_im));

% percentage:

int_Bragg_percent = 100*integrated_int_Bragg/total_int;

[tempycen_2, tempxcen_2] = ND_analysis.computeCentroids_single_CCD(dat_supergrid_toshow,pixel,'mask',mask_centroid_2);

figure(900);
subplot(122);
imagesc(log10(test_im.*mask_centroid_2));
colorbar;
axis image;
title(['x\_centroid = ' num2str(tempxcen_2) 'y\_centroid = ' num2str(tempycen_2)]);

tempycen_error = abs(tempycen_2-tempycen)/tempycen;
tempxcen_error = abs(tempxcen_2-tempxcen)/tempxcen;

figNum = 900;
namefig = ['results_paper/Figure3/error_sum_CCD' num2str(Grain_array(1)) '_pixel' num2str(pixel(1)) 'and' num2str(pixel(2))];
savefig(figure(figNum),[namefig  '.fig']);
print(figure(figNum),namefig,'-dpdf','-bestfit');
print('-r600',figure(figNum),namefig,'-dpng');


%{
counter = 1;
clear dat_supergrid_toshow;
for kk = XRF_maps_to_show
    dat_supergrid_toshow.scan(counter) = dat_supergrid.scan(kk);
    counter = counter + 1;
end
 dat_supergrid_toshow.ii = dat_supergrid.ii;

 for kkkk = 1:numel(dat_supergrid.scan)
 dat_supergrid_toshow.thetalist(kkkk) = dat_supergrid.scan(kkkk).theta;
 end

rock_curve_pixel = ND_paper_figures. displayCCDsingletheta(dat_supergrid_toshow,pixel,4);
%}