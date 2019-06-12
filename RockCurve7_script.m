if ~exist('flag_read_HXN_parameters')    
    clear all; close all;
end

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

filename = 'results/data_scan';


filename_toload = 'results/data_scan_onlyread_nopad';


%% Rocking curve 7, HXN beamtimes line 385 - 390
%%{
lst = [45477:1:45487 45449:1:45459 45460:1:45470]; 
th_middle = 84.576;
th_start = th_middle - 0.75;
th_end = th_middle + 0.75;
delta_th = 1.5/(numel(lst)-1);
th_drift_perang = -8*delta_th; % correction for drift of maps for each angle, in microns


th_step = (th_end-th_start)/(numel(lst)-1);
angs = [th_start:th_step:th_end];

del = -18.1;
gam = 10.80;%10.65;
twoTheta = 21.2;
detdist = 0.500; % in meters



ROIxstart = 1;%50;
ROIystart = 132;%81;
ROIxsize = 404;%200;
ROIysize = 380;%200;

innerpts = 41;
outerpts = 21;
innerpts_zeropad = 40;
outerpts_zeropad = 30;

if exist('flag_read_HXN_parameters')    
    return;    
end

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
% what are we going to do?
flag_read_forfirsttime = 0;
flag_make_supergrid = 1;
flag_doAlignment = 0;
flag_analyze = 0;
flag_display = 0;
flag_manual_alignement = 1;
flag_strain_analysis = 0;


if flag_read_forfirsttime == 1
    
    [dat]=ND_read_data.ThetaScan_film_onlyread(datapath,lst,XRFchan,'XBICchan',XBICchan,'thetalist',angs,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'showmerlin',0,'inneraxis',inneraxis,'do_padding',0,'do_align',0,'do_Ref_XRF0',0,'use_fitXRF',use_fitXRF);
    save([filename '_onlyread_nopad' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat','-v7.3');
    
end


%% Analisis section

if flag_make_supergrid == 1
     load([filename '_onlyread_nopad' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat'); 
    
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
            %scan = dat.scan(ii);
            %scan = rmfield(scan,chan);
            %dat.scan(ii) = scan;
        end
    end
    
    
   
    
end

if flag_doAlignment == 1
    
     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = dat.scan(ii).XRF_supergrid;
        XRF_struct.scan(ii).x_pos = xsuperGrid_s(1,:);
        XRF_struct.scan(ii).y_pos = ysuperGrid_s(:,1);
        xshifts_ini(ii) = xvalues(ii);
        yshifts_ini(ii) = 0;
     end
    
    initial_index = 1;
    final_index = numel(dat.scan);
    XRF_aligned = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'initial_index',initial_index,'final_index',final_index,'do_align',0,'do_Ref_XRF0',0);
    
    counter = 1;
    for jj = [1:numel(dat.scan)]%[initial_index:final_index]
       XRF_to_show.scan(counter) = XRF_aligned.scan(jj);
       counter = counter + 1;
    end
     ND_display_data.displayAlignedXRFmaps(XRF_to_show,'figNum',1);

    
    save([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned','-v7.3');

else
      load([filename_toload '_alignement_structure' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'XRF_aligned');

end
    

if flag_analyze == 1
    
   
    dat = ND_data_processing.diffrIntenMaps(dat,XRF_aligned,'do_supergrid',1,'delta_x_microns',delta_x_microns,'do_zeropadding',1,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad);
    
    
    [mask0,dat] = ND_data_processing.calculateMask(dat,0.1);
    
    dat.mask_rock = mask0;
    
    %%{
    pixel = [19,566]; % rows and columns, Y and X
    
    [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',mask0,'do_plot',0);
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
    
   % load([filename_toload '_analysis_aligned' num2str(lst(1)) '_' num2str(lst(end)) '.mat']);
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned,'figNum',1);
    
    ND_display_data.displayRockCurveShift(XRF_aligned,'list',angs,'title_spec','theta','figNum',4000);
    
   [dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurve(dat,'do_mask',1,'figNum',9000);

    
    %[dat.rock_curve,dat.thetalist] = ND_display_data.displayRockCurveMaps(lst,dat,'do_mask',1,'plot_XBIC',0,'size_figure',[100 100 1300 200],'figNum',9000);

     [mask0,dat] = ND_data_processing.calculateMask(dat,0.1);
    
    dat.mask_rock = mask0;
    
    mask_nan = mask0;
    
    for ii = 1:size(mask,1)
        for jj = 1:size(mask,2)
            if mask(ii,jj) == 0
                mask_nan(ii,jj) = NaN;
            end
        end
    end
    
    window_for_maps = [1,size(xsuperGrid_s,1),1,size(xsuperGrid_s,2)];

    
    figNum = 5000;
    
    
    ND_display_data.display2Dmap(dat.map2D_SumInt,'figNum',figNum,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Integrated diffr. Intensity','font',30,'size_figure',[90 151 1203 510]);
    figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
    ND_display_data.display2Dmap(centroid_struct.Xcentroids,'figNum',figNum+22,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Column Centroids','font',30,'size_figure',[90 151 1203 510]);caxis([190 210]);
    figure(figNum+22);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(centroid_struct.Ycentroids,'figNum',figNum+23,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Rows Centroids ','font',30,'size_figure',[90 151 1203 510]);caxis([167 195]);
    figure(figNum+23);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(centroid_struct.Thcentroids,'figNum',figNum+24,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Theta Centroids ','font',30,'size_figure',[90 151 1203 510]);caxis([83.7 85.0])
    figure(figNum+24);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    ND_display_data.display2Dmap(struct_centroidShift.strain.*mask0,'figNum',figNum+25,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Strain','font',30,'size_figure',[90 151 1203 510]);caxis([-2 1]*5e-3);
    figure(figNum+25);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',1.0);
    %figure(figNum+25); hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),struct_centroidShift.strain(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color','k');
   
    
    ND_display_data.display2Dmap(struct_centroidShift.tilt_x.*mask0,'figNum',figNum+26,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Tilt around x ','font',30,'size_figure',[90 151 1203 510]);caxis([-1 1]*5e-1);
    figure(figNum+26);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);
   
    ND_display_data.display2Dmap(struct_centroidShift.tilt_y.*mask0,'figNum',figNum+27,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Tilt around y ','font',30,'size_figure',[90 151 1203 510]);caxis([-1 1]*5e-1);
    figure(figNum+27);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
    ND_display_data.display2Dmap(struct_centroidShift.tilt_tot.*mask0,'figNum',figNum+28,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','Total tilt modulus ','font',30,'size_figure',[90 151 1203 510]);caxis([0 1]*1);
    figure(figNum+28);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask0(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',3.0);

    
     ND_display_data.display2Dmap(struct_centroidShift.dspace,'figNum',figNum+29,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle','D space ','font',30,'size_figure',[90 151 1203 510]);%caxis([84.28 84.5])
    %figure(figNum+29);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),'-r','LineWidth',1.0);
    figure(figNum+29); hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),struct_centroidShift.dspace(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),fliplr(distr_struct.dspace),'LineWidth',3.0,'Color','k');
 
    
    
          
end

if flag_strain_analysis == 1
    
     [mask,dat] = ND_data_processing.calculateMask(dat,0.008);
    
    dat.mask_rock = mask;
    
    %%{
    pixel = [19,566]; % rows and columns, Y and X
    
    [centroid_struct] = ND_analysis.computeCentroidsRockCurve(dat,'mask',mask,'do_plot',0);
    dat.Xcentroids = centroid_struct.Xcentroids;
    dat.Ycentroids = centroid_struct.Ycentroids;
    dat.Thcentroids = centroid_struct.Thcentroids;
    
    
    
    [struct_centroidShift] = ND_analysis.computeCentroidShiftAndStrain(dat,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,[22 543 548 201],0);
    
    
    
    contour_values_up = [2:1:11 20:10:80]*0.01;%[0.01 0.02 0.03 0.04 0.05 0.06];%
    contour_values_down = [1:1:10 10:10:70]*0.01;%[1:1:70]*0.01;
    [distr_struct,mask_struct] = ND_analysis.computeStrainContours(dat,struct_centroidShift,contour_values_up,contour_values_down);
    
    
    
    window_for_maps = [1,size(dat.xsuperGrid_s,1),1,size(dat.xsuperGrid_s,2)];

    
    color_array = ['r','k','g','b'];

    
   figNum = 5000;
    ND_display_data.display2Dmap(dat.map2D_SumInt,'Xval',dat.xsuperGrid_s(1,:),'YVal',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',[' Integrated diffraction intensity'],'font',30,'size_figure',[90 151 1203 510],'figNum',figNum);
    index_color = 1;

    for kk = 1:5:numel(mask_struct)
         

        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_down(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));
        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_up(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));

        %ND_display_data.display2Dmap(struct_centroidShift.strain.*mask_struct(kk).mask,'figNum',figNum+25+kk,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['Strain at ' num2str(contour_values(kk)) '%'],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);
        %figure(figNum+25+kk);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'-r','LineWidth',1.0);
        %caxis([-.02 .015])
        
        index_color = index_color + 1;
    end
    
    
    figNum = 6000;
    ND_display_data.display2Dmap(struct_centroidShift.dspace,'figNum',figNum,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['D space '],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);


    index_color = 1;
    
    for kk = 1:5:numel(mask_struct)
        
       %ND_display_data.display2Dmap(struct_centroidShift.dspace.*mask_struct(kk).mask,'figNum',figNum+kk,'Xval',dat.xsuperGrid_s(1,:),'Yval',dat.ysuperGrid_s(:,1),'window',window_for_maps,'figTitle',['D space at ' num2str(contour_values(kk)) '%'],'font',30,'size_figure',[90 151 1203 510]);%caxis([220 250]);

        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_down(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));
        figure(figNum);  hold on;contour(dat.xsuperGrid_s(1,window_for_maps(3):window_for_maps(4)),dat.ysuperGrid_s(window_for_maps(1):window_for_maps(2),1),mask_struct(kk).mask_up(window_for_maps(1):window_for_maps(2),window_for_maps(3):window_for_maps(4)),1,'LineWidth',3.0,'Color',color_array(index_color));

        caxis([3.7 3.77]);
        
        index_color = index_color + 1;
    end
    
    
    figure(figNum+25+kk+1);
    clf;
%     subplot(121);
%     plot(contour_values,distr_struct.strain);
%     xlabel('% of diffracted intensity map');
%     ylabel('strain \Delta d/d');
%     set(gca,'FontSize',30)
    
    %subplot(122);
    %plot(contour_values_down,distr_struct.dspace);
    
    
    x = contour_values_down; % start of bar
    y = zeros(length(x),1);
    dx = diff([x 1]); % width of bar
    dy = distr_struct.dspace;
    
    
     hold on;
    for ii=1:length(x)
        rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)])
    end
    
   index_color = 1;

    for ii=1:5:length(x)
        rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)],'FaceColor',color_array(index_color));
        index_color = index_color + 1;
    end
    
    %axis([0.5 2 0 4.1])
    
    %ylabel('Prob density')
    %xlabel('Time')
    ylim([3.7 3.77]);
    xlabel(' bins of % of diffr. inten.');
    ylabel('D [Angstroms]');
    set(gca,'FontSize',30);
    
end



if flag_manual_alignement == 1
 
index = 11;    
initial_index = 2;%index-1;
final_index = 6;%;numel(dat.scan)%index+1;

     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = XRF_aligned.scan(ii).XRF;
        XRF_struct.scan(ii).x_pos = XRF_aligned.scan(ii).x_pos;
        XRF_struct.scan(ii).y_pos = XRF_aligned.scan(ii).y_pos;
        if     ismember(ii, [initial_index:1:final_index])
            xshifts_ini(ii) = 9*(ii-initial_index);%XRF_aligned.scan(ii).xshifts-10;%(ii-initial_index-1)*5;
            yshifts_ini(ii) = 0;%XRF_aligned.scan(ii).yshifts;
        else            
            xshifts_ini(ii) = XRF_aligned.scan(ii).xshifts;
            yshifts_ini(ii) = XRF_aligned.scan(ii).yshifts;
        end
     end
    
  XRF_aligned_2 = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'initial_index',initial_index,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',1);
  
  
  initial_index = 6;%index-1;
final_index = 13;%;numel(dat.scan)%index+1;

     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = XRF_aligned.scan(ii).XRF;
        XRF_struct.scan(ii).x_pos = XRF_aligned.scan(ii).x_pos;
        XRF_struct.scan(ii).y_pos = XRF_aligned.scan(ii).y_pos;
        if     ismember(ii, [initial_index:1:final_index])
            xshifts_ini(ii) = 15*(ii-initial_index);%XRF_aligned.scan(ii).xshifts-10;%(ii-initial_index-1)*5;
            yshifts_ini(ii) = 0;%XRF_aligned.scan(ii).yshifts;
        else            
            xshifts_ini(ii) = XRF_aligned_2.scan(ii).xshifts;
            yshifts_ini(ii) = XRF_aligned_2.scan(ii).yshifts;
        end
     end
    
  XRF_aligned_2 = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'initial_index',initial_index,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',1);
    %XRF_aligned_toshow.scan(1) = XRF_aligned_2.scan(index);%XRF_aligned.scan(15);%
   XRF_aligned_toshow = XRF_aligned_2;%XRF_aligned.scan(15);%
    
   ND_display_data.displayAlignedXRFmaps(XRF_aligned_toshow,'figNum',1);    
     
     
initial_index = 13;%index-1;
final_index = 20;%;numel(dat.scan)%index+1;

     for ii = 1:numel(dat.scan)
        XRF_struct.scan(ii).XRF = XRF_aligned.scan(ii).XRF;
        XRF_struct.scan(ii).x_pos = XRF_aligned.scan(ii).x_pos;
        XRF_struct.scan(ii).y_pos = XRF_aligned.scan(ii).y_pos;
        if     ismember(ii, [initial_index:1:final_index])
            xshifts_ini(ii) = 40 + 10*(ii-initial_index);%1*(ii-initial_index);%XRF_aligned.scan(ii).xshifts-10;%(ii-initial_index-1)*5;
            yshifts_ini(ii) = 0;%XRF_aligned.scan(ii).yshifts;
        else            
            xshifts_ini(ii) = XRF_aligned_2.scan(ii).xshifts;
            yshifts_ini(ii) = XRF_aligned_2.scan(ii).yshifts;
        end
    end
    
    XRF_aligned_2 = ND_data_processing.doAlignment(XRF_struct,'scanid',lst,'thetalist',angs,'initial_index',initial_index,'xshifts_ini',xshifts_ini,'yshifts_ini',yshifts_ini,'do_align',0,'do_Ref_XRF0',1);
    %XRF_aligned_toshow.scan(1) = XRF_aligned_2.scan(index);%XRF_aligned.scan(15);%
    XRF_aligned_toshow = XRF_aligned_2;%XRF_aligned.scan(15);%
    
    ND_display_data.displayAlignedXRFmaps(XRF_aligned_toshow,'figNum',1);

    
   XRF_aligned.scan(index) =  XRF_aligned_2.scan(index);%;%
    
    

    for ii = 1:numel(XRF_aligned_2.scan)
        XRF_aligned.scan(ii).XRF = XRF_aligned_2.scan(ii).XRF;
        XRF_aligned.scan(ii).indexMap = XRF_aligned_2.scan(ii).indexMap;
        XRF_aligned.scan(ii).yshifts = XRF_aligned.scan(ii).yshifts+XRF_aligned_2.scan(ii).yshifts;
        XRF_aligned.scan(ii).xshifts = XRF_aligned.scan(ii).xshifts+XRF_aligned_2.scan(ii).xshifts;
        XRF_aligned.scan(ii).scanid = XRF_aligned_2.scan(ii).scanid;
        XRF_aligned.scan(ii).theta = XRF_aligned_2.scan(ii).theta;
        XRF_aligned.scan(ii).x_pos = XRF_aligned_2.scan(ii).x_pos;
        XRF_aligned.scan(ii).y_pos = XRF_aligned_2.scan(ii).y_pos;
    end
    ND_display_data.displayAlignedXRFmaps(XRF_aligned,'figNum',1);
    
end
%{

pixel = [23,234]; % rows and columns, Y and X


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

