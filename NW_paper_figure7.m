% In this script we prepare Figure 7 for the paper in Photovoltaics
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
    
    percent = 0.05; % percent to show
    
    
    for gg = 1:numel(Grain_array)
        
        flag_read_HXN_parameters = 1;

        eval(['RockCurve' num2str(Grain_array(gg)) '_script']);
        
        %allgrains(gg) = load(['results_paper/Figure5_grain' num2str(Grain_array(gg)) '_struct.mat']);
        load(['results_paper/Figure5_grain' num2str(Grain_array(gg)) '_struct.mat']);
        
        
        
        mask = grain_struct.map2D_SumInt > percent*max(max(grain_struct.map2D_SumInt));
        grain_struct.mask0 = mask;
          %{
        mask = allgrains(gg).grain_struct.map2D_SumInt > percent*max(max(allgrains(gg).grain_struct.map2D_SumInt));
        
      
        grain_struct.map2D_SumInt =  allgrains(gg).grain_struct.map2D_SumInt;
        grain_struct.mask = mask;%allgrains(gg).grain_struct.grain_struct.mask0;
        grain_struct.Xcentroids = allgrains(gg).grain_struct.Xcentroids;
        grain_struct.Ycentroids = allgrains(gg).grain_struct.Ycentroids;
        grain_struct.Thcentroids = allgrains(gg).grain_struct.Thcentroids;
        grain_struct.struct_centroidShift.dspace =  allgrains(gg).grain_struct.dspace;
        grain_struct.struct_centroidShift.tilt_x =  allgrains(gg).grain_struct.tilt_x;
        grain_struct.struct_centroidShift.tilt_y =  allgrains(gg).grain_struct.tilt_y;
        grain_struct.struct_centroidShift.tilt_tot =  allgrains(gg).grain_struct.tilt_tot;
        grain_struct.Xval =  allgrains(gg).grain_struct.Xval;
        grain_struct.Yval =  allgrains(gg).grain_struct.Yval;
        grain_struct.extra_shift =  allgrains(gg).grain_struct.extra_shift;
        grain_struct.window_for_maps =  allgrains(gg).grain_struct.window_for_maps;
        %}
        
        
        switch gg
            
            case 1                
                grain_struct.extra_shift = [0 100];
                grain_struct.ylim_for_tilt_tot = [0.0 0.4];
                grain_struct.ylim_for_dspace = [3.8 3.86];
            case 2
                grain_struct.ylim_for_tilt_tot = [0.15 0.7];
                grain_struct.ylim_for_dspace = [3.73 3.77];
            case 3
                grain_struct.ylim_for_tilt_tot = [0.0 0.4];
                grain_struct.ylim_for_dspace = [3.76 3.8];
            case 4
                grain_struct.extra_shift = [0 10];
                grain_struct.ylim_for_tilt_tot = [0.0 2.0];
                grain_struct.ylim_for_dspace = [3.67 3.75];
            case 5
                grain_struct.extra_shift = [0 20];
                grain_struct.ylim_for_tilt_tot = [0.0 0.2];
                grain_struct.ylim_for_dspace = [3.74 3.775];
        end
        
        
        grain_struct.Yval_range = [min(grain_struct.Yval) max(grain_struct.Yval)];
        grain_struct.Yval_step =  (max(grain_struct.Yval)-min(grain_struct.Yval))/numel(grain_struct.Yval);
        grain_struct.Yval = ([1:numel(grain_struct.Yval)]-(numel(grain_struct.Yval)/2))*grain_struct.Yval_step;
        
        grain_struct.theta_last = th_end;

        
        
        
        fields_to_plot = {'dspace'};
        
        
        fig_num = 100 + (gg-1);
        fig_number = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
            fields_to_plot,'extra_shift',grain_struct.extra_shift,...
        'Yval_lim',[-2 2],...
            'window',grain_struct.window_for_maps,'contours_to_plot',grain_struct.contours_to_plot,'plot_contours',0,...
            'min_contour',5,'spec_ylim',0,'color_array',grain_struct.color_array,'size_figure',[-97 33 1378 672],...
            'figNum',fig_num,'font',30);
        
        
        name_fig = ['results_paper/Figure6_7/' fields_to_plot{1} '_' num2str(Grain_array(gg)) ];
        savefig(figure(fig_number),[name_fig '.fig']);
        print(figure(fig_number),name_fig,'-dpdf','-bestfit');%,'PaperOrientation','landscape');
        print('-r600',figure(fig_number),name_fig,'-dpng');
        
        fields_to_plot = {'tilt_tot'};
        
        
        fig_num = 200 + (gg-1);
        fig_number = ND_paper_figures.display2DmapContoursfig5(grain_struct,...
            fields_to_plot,'extra_shift',grain_struct.extra_shift,...
        'Yval_lim',[-2 2],...
            'window',grain_struct.window_for_maps,'contours_to_plot',grain_struct.contours_to_plot,'plot_contours',0,...
            'min_contour',5,'spec_ylim',0,'color_array',grain_struct.color_array,'size_figure',[-97 33 1378 672],...
            'figNum',fig_num,'font',30);
        
        
        name_fig = ['results_paper/Figure6_7/' fields_to_plot{1} '_' num2str(Grain_array(gg)) ];
        savefig(figure(fig_number),[name_fig '.fig']);
        print(figure(fig_number),name_fig,'-dpdf','-bestfit');%,'PaperOrientation','landscape');
        print('-r600',figure(fig_number),name_fig,'-dpng');
        
        
        
    end
end

    
%{


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
ND_paper_figures.display2Dmapfig6(grain_struct_todisplay,{'dspace','tilt_tot'},'size_figure',[-97 33 1378 672]);

savefig(figure(1),['results_paper/Figure6_7/dspace_allgrains.fig']);
print(figure(1),['results_paper/Figure6_7/dspace_allgrains'],'-dpdf');
print(figure(1),['results_paper/Figure6_7/dspace_allgrains'],'-dpng');

savefig(figure(2),['results_paper/Figure6_7/total_tilt_allgrains.fig']);
print(figure(2),['results_paper/Figure6_7/total_tilt_allgrains'],'-dpdf');
print(figure(2),['results_paper/Figure6_7/total_tilt_allgrains'],'-dpng');


%}

