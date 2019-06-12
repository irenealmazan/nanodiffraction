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

load([filename_toload num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');    

[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('x_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('y_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);

dat_orig = dat;

load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

dat_supergrid = dat;

%%%% Display:

XRF_maps_to_show = [19 23];%[2:10:numel(dat.scan)];


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

ND_paper_figures.displayAlignedmaps(XRF_aligned_toshow,'XRF','figNum',1,'size_figure',[6 70 473 481],'numRows',numel(XRF_maps_to_show));

savefig(figure(1),['results_paper/Figure3/original_XRF_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(1),['results_paper/Figure3/original_XRF_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(1),['results_paper/Figure3/original_XRF_grain ' num2str(Grain_array(1))],'-dpng');

counter = 1;
clear XRF_aligned_toshow
for kk = XRF_maps_to_show
   
    XRF_aligned_toshow.scan(counter).XRF = circshift(dat_supergrid.scan(kk).XRF,[0 20]);
    XRF_aligned_toshow.scan(counter).s_pos = dat_supergrid.xsuperGrid_s(1,:);
    XRF_aligned_toshow.scan(counter).y_pos = dat_supergrid.ysuperGrid_s(:,1);
    XRF_aligned_toshow.scan(counter).z_pos = dat_supergrid.xsuperGrid_s(1,:)./tand(dat_supergrid.scan(kk).theta);
    XRF_aligned_toshow.scan(counter).theta = dat_supergrid.scan(kk).theta;

    counter = counter + 1;
end

ND_paper_figures.displayAlignedmaps(XRF_aligned_toshow,'XRF','figNum',2,'size_figure',[6 70 473 481],'numRows',numel(XRF_maps_to_show));
savefig(figure(2),['results_paper/Figure3/supergrid_XRF_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(2),['results_paper/Figure3/supergrid_XRF_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(2),['results_paper/Figure3/supergrid_XRF_grain ' num2str(Grain_array(1))],'-dpng');


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
savefig(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(3),['results_paper/Figure3/diffraction_maps_grain ' num2str(Grain_array(1))],'-dpng');


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

savefig(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(4),['results_paper/Figure3/ccd_image1_grain ' num2str(Grain_array(1))],'-dpng');


savefig(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(5),['results_paper/Figure3/ccd_image2_grain ' num2str(Grain_array(1))],'-dpng');

varnames = {['theta'] ['intensity']};
filename = ['results_paper/Figure3/rock_curve_single_pixel_grain' num2str(Grain_array(1)) '.xlsx'];
T = table(dat_supergrid_toshow.thetalist',rock_curve_pixel','VariableNames',varnames);
writetable(T,filename,'Sheet',1,'Range','D1');


% integrated diffracted intensity map and centroid maps, shifted in order
% to avoid the splitting in the figure

[~,dat_supergrid] = ND_data_processing.calculateMask(dat_supergrid,0.01);

[mask_nan] = ND_data_processing.turnMaskInNan(dat_supergrid.mask);

[centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat_supergrid,'mask',dat_supergrid.mask,'do_plot',0);
dat_supergrid.Xcentroids = centroid_struct.Xcentroids;
dat_supergrid.Ycentroids = centroid_struct.Ycentroids;
dat_supergrid.Thcentroids = centroid_struct.Thcentroids;
    
[struct_centroidShift,angles] = ND_analysis.computeCentroidShiftAndStrain(dat_supergrid,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],1);


grain_struct.extra_shift = [0 100];
grain_struct.map2D_SumInt_shift = circshift(dat_supergrid.map2D_SumInt,grain_struct.extra_shift);
grain_struct.Xcentroid_shift = circshift(dat_supergrid.Xcentroids,grain_struct.extra_shift);
grain_struct.Ycentroid_shift = circshift(dat_supergrid.Ycentroids,grain_struct.extra_shift);
grain_struct.Thcentroid_shift = circshift(dat_supergrid.Thcentroids,grain_struct.extra_shift);
grain_struct.mask_shift = circshift(dat_supergrid.mask,grain_struct.extra_shift);
grain_struct.mask_nan_shift = circshift(mask_nan,grain_struct.extra_shift);
grain_struct.angles = angles;
grain_struct.xsuperGrid_s = dat_supergrid.xsuperGrid_s;
grain_struct.ysuperGrid_s = dat_supergrid.ysuperGrid_s;
grain_struct.dspace_shift = circshift(struct_centroidShift.dspace,grain_struct.extra_shift);
grain_struct.tilt_x_shift = circshift(struct_centroidShift.tilt_x,grain_struct.extra_shift); 
grain_struct.tilt_y_shift = circshift(struct_centroidShift.tilt_y,grain_struct.extra_shift); 


window_for_maps = [1,50,1,446];
ND_display_data.display2Dmap(grain_struct.map2D_SumInt_shift,'figNum',6,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Integrated diffracted intensity','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

savefig(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(6),['results_paper/Figure3/Integrated_diff_int_grain ' num2str(Grain_array(1))],'-dpng');




ND_display_data.display2Dmap(grain_struct.Xcentroid_shift.*grain_struct.mask_nan_shift,'figNum',7,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Column Centroids','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([130 220]);

savefig(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(7),['results_paper/Figure3/Xcentroids_grain ' num2str(Grain_array(1))],'-dpng');


ND_display_data.display2Dmap(grain_struct.Ycentroid_shift.*grain_struct.mask_nan_shift,'figNum',8,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Rows Centroids ','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
%caxis([170 210]);
caxis([130 220]);

savefig(figure(8),['results_paper/Figure3/Ycentroid_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(8),['results_paper/Figure3/Ycentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(8),['results_paper/Figure3/Ycentroids_grain ' num2str(Grain_array(1))],'-dpng');



ND_display_data.display2Dmap(grain_struct.Thcentroid_shift.*grain_struct.mask_nan_shift,'figNum',9,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Theta Centroids ','font',30,'size_figure',[-8 334 785 367]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([82.0 83.60]);

savefig(figure(9),['results_paper/Figure3/Thcentroid_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(9),['results_paper/Figure3/Thcentroids_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(9),['results_paper/Figure3/Thcentroids_grain ' num2str(Grain_array(1))],'-dpng');





ND_paper_figures.displayCCDrowsColsDistribution(dat_supergrid,pixel,101)
savefig(figure(101),['results_paper/Figure3/CCD_rows_distrib_grain ' num2str(Grain_array(1)) '.fig']);
savefig(figure(102),['results_paper/Figure3/CCD_cols_distrib_grain ' num2str(Grain_array(1)) '.fig']);

% strain and tilt:


ND_display_data.display2Dmap(grain_struct.dspace_shift.*grain_struct.mask_nan_shift,'figNum',11,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Lattice spacing','font',30,'size_figure',[-8 334 785 367]);caxis([200 215]);
hold on;
contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
%caxis([130 220]);
caxis([3.78 3.86]);

savefig(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(11),['results_paper/Figure3/Lattice_grain ' num2str(Grain_array(1))],'-dpng');


%figNum = figNum_final+2;


ND_display_data.display2Dmap(grain_struct.tilt_x_shift.*grain_struct.mask_nan_shift,'figNum',13,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Azhimutal tilt','font',30,'size_figure',[-8 334 785 367]);
hold on;contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([-1.0 .6]);

savefig(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(13),['results_paper/Figure3/Tiltx_grain ' num2str(Grain_array(1))],'-dpng');



ND_display_data.display2Dmap(grain_struct.tilt_y_shift.*grain_struct.mask_nan_shift,'figNum',14,'Xval',grain_struct.xsuperGrid_s(1,:),'Yval',grain_struct.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Out-of-plane tilt','font',30,'size_figure',[-8 334 785 367]);
hold on;contour(grain_struct.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),grain_struct.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),grain_struct.mask_shift(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
caxis([-1.0 .6]);

savefig(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1)) '.fig']);
print(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1))],'-dpdf');
print(figure(14),['results_paper/Figure3/Tilty_grain ' num2str(Grain_array(1))],'-dpng');



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


