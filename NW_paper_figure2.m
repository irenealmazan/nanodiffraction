%Highlight one grain similar to Fig. 2 in JSR paper: 
%Te, Cd, Cu, Zn, XBIC, 
%total Bragg intensity, 
%focus on one pixel to show CCD image
% all this on grain 1

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


% what to do?
flag_read_HXN_parameters = 1;

Grain_array = [1];%[1];

flag_read_HXN_parameters = 1;

eval(['RockCurve' num2str(Grain_array(1)) '_script']);


% load fully aligned data
load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','XRF_aligned','centroid_struct','struct_centroidShift');

% add the xrf channels for:
 
[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('detsum_Cd_L'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
[dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('detsum_Te_L'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
  


% do supergrid for the new XRFs maps

[dat,xsuperGrid_s,ysuperGrid_s] = ND_data_processing.makeSuperGrid(dat,'detsum_Cd_L',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
[dat] = ND_data_processing.makeSuperGrid(dat,'detsum_Te_L',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);


% apply the alignement shift: 

    for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).detsum_Cd_L = dat.scan(ii).detsum_Cd_L_supergrid;
        XRF_struct.scan(ii).detsum_Te_L = dat.scan(ii).detsum_Te_L_supergrid;
        XRF_struct.scan(ii).x_pos = xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = ysuperGrid_s(:,1);
        xshifts_ini(ii) = XRF_aligned.scan(ii).xshifts;
        yshifts_ini(ii) = XRF_aligned.scan(ii).yshifts;
    end
    
    
    
XRF_aligned_Cd = ND_data_processing.doAlignment(XRF_struct,'detsum_Cd_L','scanid',lst,'thetalist',angs,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',0);
XRF_aligned_Te = ND_data_processing.doAlignment(XRF_struct,'detsum_Te_L','scanid',lst,'thetalist',angs,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',0);

save([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned','XRF_aligned_Cd','XRF_aligned_Te','-v7.3');

%{
ND_display_data.displayAlignedXRFmaps(XRF_aligned_Cd,'detsum_Cd_L','figNum',1);

XRFMaps_to_show = [2:10:numel(XRF_aligned_Te.scan) numel(XRF_aligned_Te.scan)+1];

for jj = XRFMaps_to_show 
    figure(jj);
    pause();
end

ND_display_data.displayAlignedXRFmaps(XRF_aligned_Te,'detsum_Te_L','figNum',1);

XRFMaps_to_show = [2:10:numel(XRF_aligned_Te.scan)];

for jj = XRFMaps_to_show 
    figure(jj);
    pause();
end
%}



load(['./results_paper/Figure5_grain' num2str(Grain_array(1)) '_struct.mat'],'grain_struct');

grain_struct.scan_num = numel(dat.scan);
grain_struct.detsum_Cd_L = circshift(XRF_aligned_Cd.scan(grain_struct.scan_num).detsum_Cd_L,grain_struct.extra_shift);
grain_struct.detsum_Te_L = circshift(XRF_aligned_Te.scan(grain_struct.scan_num).detsum_Te_L,grain_struct.extra_shift);
grain_struct.detsum_Cu_K = circshift(dat.scan(grain_struct.scan_num).XRF,grain_struct.extra_shift);
grain_struct.PC = circshift(dat.scan(grain_struct.scan_num).PC,grain_struct.extra_shift);

save(['./results_paper/Figure2_5_grain' num2str(Grain_array(1)) '_struct.mat'],'grain_struct','-v7.3');

return;

load(['./results_paper/Figure2_5_grain' num2str(Grain_array(1)) '_struct.mat'],'grain_struct');


fig_num = 100;
hfig = ND_paper_figures.display2DmapContoursfig5(grain_struct,{'detsum_Te_L','detsum_Cd_L','detsum_Cu_K','PC','map2D_SumInt','dspace','tilt_tot'},'extra_shift',grain_struct.extra_shift,'window',grain_struct.window_for_maps,'contours_to_plot',grain_struct.contours_to_plot,'plot_contours',0,'spec_ylim',5,'color_array',grain_struct.color_array,'size_figure',[6 60 472 691],'figNum',fig_num);


savefig(figure(100),['results_paper/Figure2_5/maps2D_' num2str(Grain_array(1)) '.fig']);
print(figure(100),['results_paper/Figure2_5/maps2D_' num2str(Grain_array(1))],'-dpdf');
print(figure(100),['results_paper/Figure2_5/maps2D_' num2str(Grain_array(1))],'-dpng');


