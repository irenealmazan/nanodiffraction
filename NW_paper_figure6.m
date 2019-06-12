% In this script we prepare Figure 6 for the paper in Photovoltaics
% journal, where we compare the strain and tilt of the maps of the
% different measured grains
%(a) d-spacing of all grains (b) statistics of all grains


clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

filename = 'results/data_scan';
filename_toload = 'results/data_scan_onlyread_nopad';

% what to do?
flag_read_first_time = 0;
flag_read_HXN_parameters = 1;

Grain_array = [1 7 9 11 12];

lst_struct(1).lst = [45193:1:45233]; % Rock curve 1
lst_struct(2).lst = [45477:1:45487 45449:1:45459 45460:1:45470]; % Rock curve 7
lst_struct(3).lst = [45526:1:45546];% Rock curve 9
lst_struct(4).lst = [45569:1:45629]; % Rock curve 11
lst_struct(5).lst = [45640:1:45650]; % Rock curve 12


percent = 0.008; % start with a very low mask to calculate the centroids and to apply a further mask for the figures

if flag_read_first_time == 1
    
    for gg = numel(Grain_array)
        load([filename_toload '_analysis_aligned' num2str(lst_struct(gg).lst(1)) '_' num2str(lst_struct(gg).lst(end)) '.mat']);
        
        grain_struct(gg).map2D_SumInt = dat.map2D_SumInt;
        
        mask0 = dat.map2D_SumInt > percent*max(max(dat.map2D_SumInt));
        
        grain_struct(gg).mask = mask0;
        
        eval(['RockCurve' num2str(Grain_array(gg)) '_script']);
        
        [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',mask0,'do_plot',0);
        grain_struct(gg).Xcentroids = centroid_struct.Xcentroids;
        grain_struct(gg).Ycentroids = centroid_struct.Ycentroids;
        grain_struct(gg).Thcentroids = centroid_struct.Thcentroids;
        
        
        
        [grain_struct(gg).struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(grain_struct(gg),twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);
        
        
        
        
        
        grain_struct(gg).Xval = dat.xsuperGrid_s(1,:);
        grain_struct(gg).Yval = dat.ysuperGrid_s(:,1);
        %grain_struct(gg).dat = dat;
        %grain_struct(gg).centroid_struct =  centroid_struct;
        %grain_struct(gg).struct_centroidShift =  struct_centroidShift;
        
        
    end
    
    save('results_paper/Figure5_grain_structure.mat','grain_struct','-v7.3');
else
    
    for gg = 1:numel(Grain_array)
        
       allgrains(gg).grain_struct = load(['results_paper/Figure5_grain' num2str(Grain_array(gg)) '_struct.mat']);
       
      grain_struct(gg).map2D_SumInt =  allgrains(gg).grain_struct.grain_struct.map2D_SumInt;
      grain_struct(gg).mask = allgrains(gg).grain_struct.grain_struct.mask0;
      grain_struct(gg).Xcentroids = allgrains(gg).grain_struct.grain_struct.Xcentroids;
      grain_struct(gg).Ycentroids = allgrains(gg).grain_struct.grain_struct.Ycentroids;
      grain_struct(gg).Thcentroids = allgrains(gg).grain_struct.grain_struct.Thcentroids;
      grain_struct(gg).struct_centroidShift.dspace =  allgrains(gg).grain_struct.grain_struct.dspace;
      grain_struct(gg).struct_centroidShift.tilt_x =  allgrains(gg).grain_struct.grain_struct.tilt_x;
      grain_struct(gg).struct_centroidShift.tilt_y =  allgrains(gg).grain_struct.grain_struct.tilt_y;
      grain_struct(gg).struct_centroidShift.tilt_tot =  allgrains(gg).grain_struct.grain_struct.tilt_tot;
      grain_struct(gg).Xval =  allgrains(gg).grain_struct.grain_struct.Xval;
      grain_struct(gg).Yval =  allgrains(gg).grain_struct.grain_struct.Yval;
      grain_struct(gg).extra_shift =  allgrains(gg).grain_struct.grain_struct.extra_shift;


    end
end

    


percent = 0.05;

for gg = 1:numel(Grain_array)
    
    mask = grain_struct(gg).map2D_SumInt > percent*max(max(grain_struct(gg).map2D_SumInt));

    
    [grain_struct_todisplay(gg).struct_centroidShift] = grain_struct(gg).struct_centroidShift; %ND_analysis.computeCentroidShiftAndStrain(grain_struct(gg),twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);
    grain_struct_todisplay(gg).Xval = grain_struct(gg).Xval;
    grain_struct_todisplay(gg).Yval = grain_struct(gg).Yval;
    grain_struct_todisplay(gg).mask = mask; 
    grain_struct_todisplay(gg).extra_shift = grain_struct(gg).extra_shift; 
   
    if gg == 1
        grain_struct_todisplay(gg).extra_shift = [0 100];
    end
    
end



%ND_paper_figures.display2Dmapfig6(grain_struct_todisplay,{'strain','tilt_x','tilt_y','tilt_tot'});
ND_paper_figures.display2Dmapfig6(grain_struct_todisplay,{'dspace','tilt_tot'},'size_figure',[6 60 472 691]);

savefig(figure(1),['results_paper/Figure6_7/dspace_allgrains.fig']);
print(figure(1),['results_paper/Figure6_7/dspace_allgrains'],'-dpdf');
print(figure(1),['results_paper/Figure6_7/dspace_allgrains'],'-dpng');

savefig(figure(2),['results_paper/Figure6_7/total_tilt_allgrains.fig']);
print(figure(2),['results_paper/Figure6_7/total_tilt_allgrains'],'-dpdf');
print(figure(2),['results_paper/Figure6_7/total_tilt_allgrains'],'-dpng');




