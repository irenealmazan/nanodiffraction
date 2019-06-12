% this script reads the data, aligns the diffraction maps and the XBIC maps
% using the fluorescence map, displays the rocking curve and calculates the
% strain map

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

filename = 'results/data_scan';
filename_toload = 'results/data_scan_onlyread_nopad';

%% Rocking curve 9, HXN beamtimes line 461 and Rocking curve 10, HXN beamtimes line 465
%lst1_theo = [45346:1:45366];
%%{
lst1 = [45526:1:45546]; 
lst2 = [45547:1:45567];

lst = [];
for ii = 1:numel(lst1)
    lst = [lst lst1(ii) lst2(ii)];
end

th_middle = 84.426;
th_start = th_middle - 0.5;
th_end = th_middle + 0.5;
delta_th = 1.0/(numel(lst)-1);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns

th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];


del = -16.1;
gam = 13.1;%10.65;
twoTheta = 21.2;
detdist = 0.35; % in meters


ROIxstart = 89;%50;
ROIystart = 6;%81;
ROIxsize = 262;%200;
ROIysize = 506;%200;

innerpts = 60; % columns
outerpts = 20; % rows
innerpts_zeropad = 10;
outerpts_zeropad = 10;

inneraxis = 'z';


inneraxis = 'z';

use_fitXRF = 1;
if use_fitXRF == 0
    XRFchan = 'Det1_Cu';
else
    XRFchan = 'detsum_Cu_K_';
end
XBICchan = 'sclr1_ch3';
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};



eval('Init_parameters');

% What are we going to do?
flag_read_forfirsttime = 0;
flag_make_supergrid = 1;
flag_doAlignment = 0;
flag_analyze = 1;
flag_display = 1;
flag_manual_alignement = 0;

if flag_read_forfirsttime == 1

 [dat]=ND_read_data.ThetaScan_film_onlyread(datapath,lst,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',0,'do_align',0,'do_Ref_XRF0',0,'use_fitXRF',use_fitXRF);
 save([filename '_onlyread_nopad' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','-v7.3');
 
else

    load([filename_toload num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');    
 
end


%% Analisis section

if flag_make_supergrid == 1
    delta_x_microns = 0.1; % step in the sample frame in microns
    
    [dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('x_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
    [dat] = ND_read_data.read_tiff_and_getLinearDatain2DMap(datapath,char('y_pos'),lst,dat,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'do_padding',0);
    [dat] = ND_data_processing.convertHXNtoSampleFrame(dat);
    
    
    [dat,xsuperGrid_s,ysuperGrid_s] = ND_data_processing.makeSuperGrid(dat,'XRF',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    [dat] = ND_data_processing.makeSuperGrid(dat,'PC',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
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
    
    
    for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = dat.scan(ii).XRF_supergrid;
        XRF_struct.scan(ii).x_pos = xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = ysuperGrid_s(:,1);
        xshifts_ini(ii) = 0;
        yshifts_ini(ii) = 0;
    end
    
end

if flag_doAlignment == 1
    
    XRF_aligned = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',1,'do_Ref_XRF0',0);
    save([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned','-v7.3');

else
      load([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat']);

end
    

if flag_analyze == 1
    
    
    
    for ii = 1:numel(dat.scan)
        XRF_aligned.scan(ii).x_pos = xsuperGrid_s(1,:);%dat.scan(ii).x_pos_supergrid;
        XRF_aligned.scan(ii).y_pos = ysuperGrid_s(:,1);%dat.scan(ii).y_pos_supergrid;
       % xshifts_ini(ii) = 0;
       % yshifts_ini(ii) = 0;
    end
    
    

    
    dat = ND_data_processing.diffrIntenMaps(dat,XRF_aligned,'do_supergrid',1,'delta_x_microns',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    
    
    [mask,dat] = ND_data_processing.calculateMask(dat,0.1);
    
    dat.mask_rock = mask;
    
    %%{
    pixel = [19,566]; % rows and columns, Y and X
    
    [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',mask,'do_plot',0);
    dat.Xcentroids = centroid_struct.Xcentroids;
    dat.Ycentroids = centroid_struct.Ycentroids;
    dat.Thcentroids = centroid_struct.Thcentroids;
    
    
    
    [struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(dat,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],1);
    
    save([filename_toload '_analysis' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','XRF_aligned','centroid_struct','struct_centroidShift','-v7.3');
   
%else
 %  load([filename_toload '_analysis' num2str(lst(1)) '_' num2str(lst(end)) '.mat']);

end

%% Display results:


if flag_display == 1
    
    load([filename_toload '_analysis' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned');
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned,'figNum',1);
    
    ND_display_data.displayRockCurveShift(XRF_aligned,'list',angs,'title_spec','theta','figNum',4000);
    
   [dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurve(dat,'do_mask',1,'figNum',9000);

    
    %[dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurveMaps(lst,dat,'do_mask',1,'plot_XBIC',0,'size_figure',[100 100 1300 200],'figNum',9000);

    
    window_for_maps = [1,34,173,690];
    
      
    [mask,dat] = ND_data_processing.calculateMask(dat,0.1);
    
    dat.mask_rock = mask;
    
    ND_display_data.display2Dmap(dat.map2D_SumInt,'figTitle','Integrated diffraction intensity','figNum',5000);
    figure(5000);  hold on;contour(mask,'-r','LineWidth',3.0);
    
    
    ND_display_data.display2Dmap(centroid_struct.Xcentroids,'figNum',22,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'figTitle','Row Centroids','size_figure',[22 422 395 322]);caxis([220 240]);
    ND_display_data.display2Dmap(centroid_struct.Ycentroids,'figNum',23,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'figTitle','Column Centroids ','size_figure',[22 422 395 322]);caxis([110 135]);
    ND_display_data.display2Dmap(centroid_struct.Thcentroids,'figNum',24,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'figTitle','Theta Centroids ','size_figure',[22 422 395 322]);caxis([84 84.6]);
    
end




if flag_manual_alignement == 1
    
initial_index = 34;

     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = dat.scan(ii).XRF_supergrid;
        XRF_struct.scan(ii).x_pos = dat.xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = dat.ysuperGrid_s(:,1);
        if ii > initial_index
            xshifts_ini(ii) = XRF_aligned_2.scan(ii).xshifts-2;%(ii-initial_index-1)*5;
            yshifts_ini(ii) = XRF_aligned_2.scan(ii).yshifts;
        else            
            xshifts_ini(ii) = XRF_aligned.scan(ii).xshifts;
            yshifts_ini(ii) = XRF_aligned.scan(ii).yshifts;
        end
    end
    
    XRF_aligned_3 = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'initial_index',initial_index,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',1,'do_Ref_XRF0',0);
    
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned_3,'figNum',1);
    ND_display_data.displayRockCurveShift(XRF_aligned_3,'list',angs,'title_spec','theta','figNum',4000);


end
%{

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
plot(dat.Ycentroids(pixel(1),pixel(2)),dat.Xcentroids(pixel(1),pixel(2)),'or')
axis image;
colormap hot;
colorbar;
title(['summed ccd at pixel col = ' num2str(pixel(2)) ' row = ' num2str(pixel(1))]);

figure(103);
clf;
plot(dat.thetalist,dat.ii(pixel(1)).jj(pixel(2)).intensity,'LineWidth',3.0);
hold on;
plot(dat.Thcentroids(pixel(1),pixel(2)),dat.Thcentroids(pixel(1),pixel(2)),'or')
title(['Rock. curve at pixel col = ' num2str(pixel(2)) ' row = ' num2str(pixel(1))]);


do_correlation = 1;

if do_correlation 
    %lst9 = [45526:1:45546];
    %lst10 = [45547:1:45567]; 
    filenames = {['results/dat_Rock9.mat']...
        ,['results/dat_Rock10.mat'],...
        ['results/data_scan_alignXRF0' num2str(lst(1)) '_' num2str(lst(end)) '.mat']};
   Rtest = ND_analysis.calculateCorrCoef(dat,mask,filenames,1);
end


%}