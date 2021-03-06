% this script reads the data, aligns the diffraction maps and the XBIC maps
% using the fluorescence map, displays the rocking curve and calculates the
% strain map

%if ~exist('flag_read_HXN_parameters')    
%    clear all; close all;
%end


addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

filename = 'results/data_scan';
filename_toload = 'results/data_scan_onlyread_nopad';



%% what are we going to do?
flag_read_forfirsttime = 0;
flag_make_supergrid = 0;
flag_doAlignment = 0;
flag_analyze = 0;
flag_display = 0;
flag_strain_analysis = 1;
flag_manual_alignement = 0;
flag_median = 0;

RockCurve1_script;
%RockCurve7_script;
%RockCurve9_script;
%RockCurve11_script;
%RockCurve12_script;

%if exist('flag_read_HXN_parameters')    
%    return;    
%end



inneraxis = 'z';

use_fitXRF = 1;
if use_fitXRF == 0
    XRFchan = 'Det1_Cu';
else
    XRFchan = 'detsum_Cu_K_';
end
normchan = 'sclr1_ch4';
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

struct_median_I0 =  load('median_I0_grain1.mat');
median_I0 = struct_median_I0.median_I0;

eval('Init_parameters');

%return;


if flag_read_forfirsttime == 1
    [dat]=ND_read_data.ThetaScan_film_onlyread(datapath,lst,XRFchan,normchan,median_I0,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',0,'do_align',0,'do_Ref_XRF0',0,'use_fitXRF',use_fitXRF);
    save([filename '_onlyread_nopad' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','-v7.3');
else
    
    load([filename_toload num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');
end



%% Analisis section

if flag_make_supergrid == 1
%     if exist('dat','var')
%         load([filename_toload num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');
%     end
    
    delta_x_microns = 0.1; % step in the sample frame in microns
    dat.delta_x_microns = delta_x_microns;
    
    [dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('x_pos'),lst,dat,'donorm',0,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
    [dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('y_pos'),lst,dat,'donorm',0,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
   
    [dat] = ND_data_processing.convertHXNtoSampleFrame(dat);
    
    [dat] = ND_data_processing.makeSuperGrid(dat,'XRF',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    [dat,xsuperGrid_s,ysuperGrid_s] = ND_data_processing.makeSuperGrid(dat,'PC',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    %[dat,xsuperGrid_s,ysuperGrid_s] = ND_data_processing.makeSuperGrid(dat,'sclr1_ch4',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);

    for ii = 1:numel(dat.scan)
        for jj = 1:size(dat.scan(1).dataout,3)
            chan = ['dataout' num2str(jj)];
            dat.scan(ii).(chan) = dat.scan(ii).dataout(:,:,jj);
        end
    end
    
    for jj = 1:size(dat.scan(1).dataout,3)
        chan = ['dataout' num2str(jj)];
        [dat,xsuperGrid_s,ysuperGrid_s] = ND_data_processing.makeSuperGrid(dat,chan,delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    end
    
    for ii = 1:numel(dat.scan)
        for jj = 1:size(dat.scan(1).dataout,3)
            chan = ['dataout' num2str(jj) '_supergrid'];
            dat.scan(ii).dataout_supergrid(:,:,jj) = dat.scan(ii).(chan);
        end
    end
    
    
  
    
end

if flag_median == 1
    [dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('sclr1_ch4'),lst,dat,'donorm',0,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
 
     if strcmp(Grain_label,'Grain1')
       median_I0 = median(dat.scan(32).sclr1_ch4(:));
       dat.median_I0 = median_I0;
       save('median_I0_grain1.mat','median_I0');
    else
       dat.median_I0 =  load('median_I0_grain1.mat');
    end
end

if flag_doAlignment == 1
    
    for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = dat.scan(ii).XRF_supergrid;
        XRF_struct.scan(ii).x_pos = xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = ysuperGrid_s(:,1);
        xshifts_ini(ii) = 0;
        yshifts_ini(ii) = 0;
    end
    
    XRF_aligned = ND_data_processing.doAlignment(XRF_struct,'XRF','scanid',lst,'thetalist',angs,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',1,'do_Ref_XRF0',0);
    save([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned','-v7.3');

else
      load([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned');

end
    

if flag_analyze == 1
    
    
    
    for ii = 1:numel(dat.scan)
        XRF_aligned.scan(ii).x_pos = dat.xsuperGrid_s(1,:);%dat.scan(ii).x_pos_supergrid;
        XRF_aligned.scan(ii).y_pos = dat.ysuperGrid_s(:,1);%dat.scan(ii).y_pos_supergrid;
       % xshifts_ini(ii) = 0;
       % yshifts_ini(ii) = 0;
    end
    
    

    
    dat = ND_data_processing.diffrIntenMaps(dat,XRF_aligned,'do_supergrid',1,'delta_x_microns',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    
    
    [mask_nan,dat] = ND_data_processing.calculateMask(dat,0.05);
    
    %dat.mask_rock = mask0;
    
    %%{
    pixel = [19,566]; % rows and columns, Y and X
    
    [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',dat.mask,'do_plot',0);
    dat.Xcentroids = centroid_struct.Xcentroids;
    dat.Ycentroids = centroid_struct.Ycentroids;
    dat.Thcentroids = centroid_struct.Thcentroids;
    
    
    
    [struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(dat,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);
    
    save([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','XRF_aligned','centroid_struct','struct_centroidShift','-v7.3');
   
else
   load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','centroid_struct','struct_centroidShift');

end

%% Display results:



if flag_display == 1
    
    window_for_maps = [1,50,15,350];
    
   % load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat']);
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned,'figNum',1);
    
    ND_display_data.displayRockCurveShift(XRF_aligned,'list',angs,'title_spec','theta','figNum',4000);
    
   [dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurve(dat,'do_mask',1,'figNum',9000);

    
    [dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurveMaps(lst,dat,'do_mask',1,'plot_XBIC',1,'window',window_for_maps,'size_figure',[100 100 1300 200],'figNum',9000);

    [dat.mask_nan] = ND_data_processing.turnMaskInNan(dat.mask);
    
    

   

    
    figNum = 5000;
    
     ND_display_data.display2Dmap(dat.map2D_SumInt.*mask_nan,'figNum',figNum,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Integrated diffracted intensity','font',30,'size_figure',[90 151 1203 510]);
    
     window_for_maps = [1,49,14,349];
     figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),dat.mask(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
     figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),dat.mask_rock(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-k','LineWidth',3.0);

    
    ND_display_data.display2Dmap(centroid_struct.Xcentroids,'figNum',figNum+22,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Column Centroids','font',30,'size_figure',[90 151 1203 510]);caxis([200 215]);
    figure(figNum+22);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(centroid_struct.Ycentroids,'figNum',figNum+23,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Rows Centroids ','font',30,'size_figure',[90 151 1203 510]);caxis([192 207]);
    figure(figNum+23);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(centroid_struct.Thcentroids.*mask_nan,'figNum',figNum+24,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Theta Centroids ','font',30,'size_figure',[90 151 1203 510]);caxis([82.6 82.85])
    figure(figNum+24);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(struct_centroidShift.strain.*mask_nan,'figNum',figNum+25,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Strain','font',30,'size_figure',[90 151 1203 510]);caxis([-0.04 0.025]*1e-1);
    figure(figNum+25);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
    ND_display_data.display2Dmap(struct_centroidShift.tilt_x.*mask_nan,'figNum',figNum+26,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Tilt around x ','font',30,'size_figure',[90 151 1203 510]);caxis([-0.2 0.2]);
    figure(figNum+26);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
   
    ND_display_data.display2Dmap(struct_centroidShift.tilt_y.*mask_nan,'figNum',figNum+27,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Tilt around y ','font',30,'size_figure',[90 151 1203 510]);caxis([-0.1 0.1]*1);
    figure(figNum+27);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
    ND_display_data.display2Dmap(struct_centroidShift.tilt_tot.*mask_nan,'figNum',figNum+28,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Total tilt modulus ','font',30,'size_figure',[90 151 1203 510]);%caxis([84.28 84.5])
    figure(figNum+28);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
     ND_display_data.display2Dmap(struct_centroidShift.dspace.*mask_nan,'figNum',figNum+29,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','D space ','font',30,'size_figure',[90 151 1203 510]);%caxis([84.28 84.5])
    %figure(figNum+29);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),'-r','LineWidth',1.0);
    figure(figNum+29);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_nan(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
 
    
    
          
end

if flag_strain_analysis == 1
    
    [mask_nan,dat] = ND_data_processing.calculateMask(dat,0.008);
    
    dat.mask_rock = mask_nan;
    
    %%{
    pixel = [19,566]; % rows and columns, Y and X
    
    [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',mask_nan,'do_plot',0);
    dat.Xcentroids = centroid_struct.Xcentroids;
    dat.Ycentroids = centroid_struct.Ycentroids;
    dat.Thcentroids = centroid_struct.Thcentroids;
    
    
    
    [struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(dat,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);
    
    
    
    contour_values_up = [2:1:100]*0.01;%%[0.01 0.02 0.03 0.04 0.05 0.06];%
    contour_values_down = [1:1:99]*0.01;%[1:1:70]*0.01;
    [distr_struct_dspace,mask_struct] = ND_analysis.computeStrainOrTiltContours(dat,struct_centroidShift,'dspace',contour_values_up,contour_values_down);
    distr_struct.dspace = distr_struct_dspace.distr;
    distr_struct.dspace_sigma = distr_struct_dspace.sigma;
     
    [distr_struct_tilt_tot,mask_struct] = ND_analysis.computeStrainOrTiltContours(dat,struct_centroidShift,'tilt_tot',contour_values_up,contour_values_down);
     distr_struct.tilt_tot = distr_struct_tilt_tot.distr;
    distr_struct.tilt_tot_sigma = distr_struct_tilt_tot.sigma;
    
    [distr_struct_grad,mask_struct] = ND_analysis.computeGradTiltContours(dat,struct_centroidShift,contour_values_up,contour_values_down);
    distr_struct.grad_tilt_tot = distr_struct_grad.distr_grad_tilt_tot;
    distr_struct.grad_tilt_tot_sigma = distr_struct_grad.sigma_grad_tilt_tot;
    
    window_for_maps = [1,size(dat.xsuperGrid_s,1),1,size(dat.xsuperGrid_s,2)];

    
    color_array = ['r','k','g','b','m','r','k','g','b','m','r','k','g','b','m','r','k','g','b','m'];

    
    contour_min_to_plot = 5;
    [mask_toshow,dat] = ND_data_processing.calculateMask(dat,0.08);%contour_values_down(contour_min_to_plot)

    [mask_toshow_nan] = ND_data_processing.turnMaskInNan(double(mask_toshow));

    
    dat.mask_toshow_nan = mask_toshow_nan;
    
    
   figNum = 5000;
   ND_display_data.display2Dmap(dat.map2D_SumInt.*mask_toshow_nan,'Xval',dat.xsuperGrid_s(1,:),'YVal',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',[' Integrated diffraction intensity'],'font',30,'size_figure',[90 151 1203 510],'figNum',figNum);
   caxis([min(min(dat.map2D_SumInt.*mask_toshow_nan)) max(max(dat.map2D_SumInt.*mask_toshow_nan))]);

    
    index_color = 1;

    %{
    for kk = contour_min_to_plot:15:numel(mask_struct)
         

        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_down(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));
        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_up(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));

        %ND_display_data.display2Dmap(struct_centroidShift.strain.*mask_struct(kk).mask,'figNum',figNum+25+kk,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['Strain at ' num2str(contour_values(kk)) '%'],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);
        %figure(figNum+25+kk);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',1.0);
        %caxis([-.02 .015])
        
        index_color = index_color + 1;
    end
    %}
    
    figNum = 6000;
    ND_display_data.display2Dmap(struct_centroidShift.dspace.*mask_toshow_nan,'figNum',figNum,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['D space '],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);
   caxis([min(min(struct_centroidShift.dspace.*mask_toshow_nan)) max(max(struct_centroidShift.dspace.*mask_toshow_nan))]);


    index_color = 1;
  
    %{
    for kk = contour_min_to_plot:15:numel(mask_struct)
        
       %ND_display_data.display2Dmap(struct_centroidShift.dspace.*mask_struct(kk).mask,'figNum',figNum+kk,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['D space at ' num2str(contour_values(kk)) '%'],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);

        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_down(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));
        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_up(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));

        %caxis([3.78 3.86]);
        
        index_color = index_color + 1;
    end  
    %}
    
    figNum = 7000;
    ND_display_data.display2Dmap(gradient(struct_centroidShift.tilt_tot).*mask_toshow_nan,'figNum',figNum,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['Total tilt'],'font',30,'size_figure',[90 151 1203 510]);
    caxis([min(min(gradient(struct_centroidShift.tilt_tot).*mask_toshow_nan)) max(max(gradient(struct_centroidShift.tilt_tot).*mask_toshow_nan))]);


    index_color = 1;
   
    %{
    for kk = contour_min_to_plot:15:numel(mask_struct)
        
       %ND_display_data.display2Dmap(struct_centroidShift.dspace.*mask_struct(kk).mask,'figNum',figNum+kk,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['D space at ' num2str(contour_values(kk)) '%'],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);

        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_down(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));
        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_up(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));

        %caxis([3.78 3.86]);
        
        index_color = index_color + 1;
    end
    %}
    figNum = 8000;
    
    figure(figNum);
    clf;
    subplot(131);
    errorbar(contour_values_down,distr_struct.dspace,distr_struct.dspace_sigma)
    ylim([3.78 3.9]);
    xlabel(' bins of % of diffr. inten.');
    ylabel('D [Angstroms]');
    set(gca,'FontSize',20);
    title('dspace vs int XRD');
    
    subplot(132);
    errorbar(contour_values_down,distr_struct.tilt_tot,distr_struct.tilt_tot_sigma)
    ylim([0.0 1.0]);
    xlabel(' bins of % of diffr. inten.');
    ylabel('total tilt [degrees]');
    set(gca,'FontSize',20);
    title('total tilt vs int XRD');
    
    subplot(133);
    errorbar(contour_values_down,distr_struct.grad_tilt_tot,distr_struct.grad_tilt_tot_sigma)
     ylim([0.0 1.0]);
    xlabel(' bins of % of diffr. inten.');
    ylabel('grad total tilt');
    set(gca,'FontSize',20);
    title('grad. of total tilt vs int XRD');
    
    
    
end




if flag_manual_alignement == 1
    
initial_index = numel(dat.scan);
final_index = numel(dat.scan);

     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = dat.scan(ii).XRF_supergrid;
        XRF_struct.scan(ii).x_pos = dat.xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = dat.ysuperGrid_s(:,1);
        if ii > initial_index
            xshifts_ini(ii) = XRF_aligned.scan(ii).xshifts;%(ii-initial_index-1)*5;
            yshifts_ini(ii) = XRF_aligned.scan(ii).yshifts;
        else            
            xshifts_ini(ii) = XRF_aligned.scan(ii).xshifts+150;
            yshifts_ini(ii) = XRF_aligned.scan(ii).yshifts;
        end
    end
    
    XRF_aligned_2 = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'initial_index',initial_index,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',0);
    
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned,'figNum',1);

end


%return;

%{
pixel = [48,48]; % rows and columns, Y and X


ROICenter = [ ...
                ROIxstart + ROIxsize / 2 ...
                Ndet - ( ROIystart + ROIysize / 2 ) ...
                ];
            
ROIxdisplay = ROIxsize / 6;  
ROIydisplay = ROIysize / 6;  
            
figure(102);
clf;
imagesc(log10(dat.ii(pixel(1)).jj(pixel(2)).im));% (ROICenter(2)-ROIydisplay:ROICenter(2)+ROIydisplay-1,ROICenter(1)-ROIxdisplay:ROICenter(1)+ROIxdisplay-1));%(ROIxstart:ROIxsize+ROIxstart,ROIystart:ROIysize+ROIystart));
hold on;
plot(ROICenter(1),ROICenter(2),'or')
axis image;
colormap hot;
colorbar;
title(['summed ccd']);





%}