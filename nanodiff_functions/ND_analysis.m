classdef ND_analysis
    % This library contains all the functions to analyse the
    % nanodiffraction pattern in order to extract the strain and the tilt
    % maps
    properties(Constant)
    end
    
    
    methods(Static)
   
        function [centroid_struct,pixel_show] = computeCentroidsRockCurve(fly2Dmaps,varargin)
            % this funtion computes the centroid for each pixel of a rocking curve
            % (over many angles)
            
            
            p = inputParser;
            
            addRequired(p,'fly2Dmaps',@isstruct)
            addParameter(p,'mask',[],@ismatrix)
            addParameter(p,'thetalist',[],@isnumeric)
            addParameter(p,'do_plot',0,@isnumeric)
            addParameter(p,'show_pixel_index',[],@isnumeric)
            
            parse(p,fly2Dmaps,varargin{:})
            
            fly2Dmaps = p.Results.fly2Dmaps;
            mask = p.Results.mask;
            
            if isempty(mask)
                mask = ones(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            end
            
            if isempty(p.Results.thetalist)
               for ii = 1:numel(fly2Dmaps.scan)
                  thetalist(ii) = fly2Dmaps.scan(ii).theta; 
               end
            else
                thetalist = p.Results.thetalist;
            end
            
            tempim = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            tempxcen = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            tempycen = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            tempthcen = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            
            for kk = 1:numel(fly2Dmaps.ii)
                for ll = 1:numel(fly2Dmaps.ii(1).jj)
                    if mask(kk,ll) == 1
                        tempim(kk,ll) = fly2Dmaps.ii(kk).jj(ll).SumInt; %sum_angle sum(sum(ccd))
                        imgin =fly2Dmaps.ii(kk).jj(ll).im; %sum_angle ccd_angle -> 512 x 512 image
                        imgin_theta = squeeze(fly2Dmaps.ii(kk).jj(ll).intensity); % sum(sum(ccd))_angle -> array with number of angles in rock curve entries
                        
                        % vertical
                        line1=sum(imgin,1);  
                        sumt1 = sum(line1);
                        
                        for mm=1:size(line1,2)
                            tempycen(kk,ll)=tempycen(kk,ll)+mm*line1(mm)/ sumt1 ;
                        end
                        
                        % horizontal
                        line2=sum(imgin,2);  
                        sumt2 = sum(line2);
                        for mm=1:size(line2,1)
                            tempxcen(kk,ll)=tempxcen(kk,ll)+mm*line2(mm)/sumt2;
                        end
                        
                        sumt = sum(imgin_theta);
                        for tt = 1:size(imgin_theta,2)
                            tempthcen(kk,ll) = tempthcen(kk,ll) + thetalist(tt)*imgin_theta(tt)/sumt;
                        end
                        
                        if ~isempty(p.Results.show_pixel_index)
                            if kk == p.Results.show_pixel_index(1)
                                if ll ==  p.Results.show_pixel_index(2)
                                    pixel_show.line1 = line1/sumt1;
                                    pixel_show.line2 = line2/sumt2;
                                end
                            end
                        else
                            pixel_show.line1 = 0;
                            pixel_show.line2 = 0;
                        end
                    else
                       pixel_show.line1 = 0;
                       pixel_show.line2 = 0; 
                    end
                end
            end
            
          
            
            centroid_struct.Xcentroids = tempxcen; % pixel units
            centroid_struct.Ycentroids = tempycen;
            centroid_struct.Thcentroids = tempthcen;
            
            if p.Results.do_plot
                ND_display_data.display2Dmap(centroid_struct.Xcentroids,'figNum',22,'Xval',fly2Dmaps.xsuperGrid_s(1,:),'Yval',fly2Dmaps.ysuperGrid_s(:,1),'figTitle','Row Centroids');
                ND_display_data.display2Dmap(centroid_struct.Ycentroids,'figNum',23,'Xval',fly2Dmaps.imapx,'Yval',fly2Dmaps.imapy,'figTitle',['Column Centroids ']);
                ND_display_data.display2Dmap(centroid_struct.Thcentroids,'figNum',24,'Xval',fly2Dmaps.imapx,'Yval',fly2Dmaps.imapy,'figTitle',['Theta Centroids ']);
            end
            
        end
        
        function [struct_centroidShift,angles] = computeCentroidShiftAndStrain(fly2Dmaps,twoTheta,del,gam,detdist,ROIxstart,ROIxsize,ROIystart,ROIysize,fig_position,plotflag)
            
            eval('Init_parameters');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Transformation for XCentroid and YCentroid from matlab convention to
            %%% detector space convention as seen on Merlin.
            
            ROICenter = [ ...
                ROIxstart + ROIxsize / 2 ...
                Ndet - ( ROIystart + ROIysize / 2 ) ...
                ];
            
            detectorCenter = [ 1 1 ] * ( Ndet/2 );
            ROIOffset = ROICenter - detectorCenter;
            
            kb = 2*pi*Ekev/12.39842; % beam momentum 1/A
            lam_A = 2*pi/kb;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%
            %%{
            %CIGS_diffraction_detector;
            angles = ND_analysis.qmatrix_MAR(Ekev,del,gam,0,0,detdist,Ndet,pixsize,0);
            
            radial_D = angles.anglemat(:,:,1);
            azi_D = angles.anglemat(:,:,2);
            
            
            ROI_radial = radial_D(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
            ROI_azi = azi_D(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
           
            % in 1/A
            ROI_qmat.qmat_1 = angles.qmat(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1,1); % recip_x component
            ROI_qmat.qmat_2 = angles.qmat(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1,2); % recip_y component
            ROI_qmat.qmat_3 = angles.qmat(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1,3); % recip_z component
            ROI_qmat.qmat_4 = angles.qmat(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1,4); % magnitude
            
            xcen = round(fly2Dmaps.Xcentroids(:));
            ycen = round(fly2Dmaps.Ycentroids(:));
            thcen = twoTheta - fly2Dmaps.Thcentroids(:);
            
            radial_proj = zeros(size(xcen));
            azi_proj = zeros(size(xcen));
            theta_proj = zeros(size(thcen));
            mask = fly2Dmaps.mask(:);%zeros(size(xcen));
            
            meanval_x = [];  %calculate the mean centroid value, excluding fluo
            meanval_y = [];
            meanval_theta_proj = [];
            
            for ii = 1:numel(xcen)
                if  mask (ii) == 1;%xcen(ii) ~=0 && isnan(xcen(ii)) == 0 && ycen(ii) ~=0 && isnan(ycen(ii)) == 0
                    %mask(ii) = 1;
                    
                    radial_proj(ii) = ROI_radial(ycen(ii),xcen(ii)); % radial component
                    azi_proj(ii) = ROI_azi(ycen(ii),xcen(ii));% azimutal component           
                    theta_proj(ii) = thcen(ii);
                    
                    %if fly2Dmaps.mask(ii)
                        %if XRFtemp(ii)>minXRF
                        meanval_x = [meanval_x, radial_proj(ii)];
                        meanval_y = [meanval_y, azi_proj(ii)];
                        meanval_theta_proj = [meanval_theta_proj, theta_proj(ii)];
                   % end
                end
                
            end
            
            
            meanval_radial = mean(meanval_x);
            meanval_azi = mean(meanval_y);
            meanval_theta = mean(meanval_theta_proj);
                        
            radial_proj = reshape(radial_proj,size(fly2Dmaps.Xcentroids) );
            azi_proj = reshape(azi_proj, size(fly2Dmaps.Xcentroids));
            theta_proj = reshape(theta_proj, size(fly2Dmaps.Xcentroids));
            mask = reshape(mask, size(fly2Dmaps.Xcentroids));
            xcen = reshape(xcen, size(fly2Dmaps.Xcentroids));%reshape(xcen, size(fly2Dmaps.Xcentroids));
            ycen = reshape(ycen, size(fly2Dmaps.Xcentroids));
            thcen = reshape(thcen, size(fly2Dmaps.Xcentroids));
            meanval_radial_matrix = meanval_radial .*ones(size(fly2Dmaps.Xcentroids) );
            
            
            
   
            dspace_mean = lam_A./(2*sind( meanval_radial_matrix/2));
            dspace_strain = lam_A./(2*sind((radial_proj)/2));
            
            strain = (dspace_strain-dspace_mean)./dspace_mean;
            
            tilt = azi_proj - meanval_azi;
            delta_theta_proj = theta_proj -   meanval_theta ;
            tilt_tot = sqrt( tilt.^2 + delta_theta_proj.^2);
            
            
            struct_centroidShift.angles = angles;
            struct_centroidShift.dspace = dspace_strain;
            struct_centroidShift.tilt_x = tilt;
            struct_centroidShift.tilt_y = delta_theta_proj;
            struct_centroidShift.tilt_tot = tilt_tot;
            struct_centroidShift.strain = strain;
            struct_centroidShift.radial_proj = radial_proj;
            struct_centroidShift.azi_proj = azi_proj;
            struct_centroidShift.theta_proj = theta_proj;
            struct_centroidShift.xcen = xcen;
            struct_centroidShift.ycen = ycen;
            struct_centroidShift.thcen = thcen;
            struct_centroidShift.meanval_radial = meanval_radial;
            
            %%{
            if plotflag
                
             
                figure(7); clf; imagesc(fly2Dmaps.imapx(1,:),fly2Dmaps.imapy(:,1),strain.*mask); colorbar; title('strain');axis image;
                set(gcf,'Position',fig_position);
                figure(8); clf; imagesc(fly2Dmaps.imapx(1,:),fly2Dmaps.imapy(:,1),tilt.*mask); colorbar; title('tilt around x');axis image;
                set(gcf,'Position',fig_position);
                figure(9); clf; imagesc(fly2Dmaps.imapx(1,:),fly2Dmaps.imapy(:,1),delta_theta_proj.*mask); colorbar; title('tilt around y');axis image;
                set(gcf,'Position',fig_position);
                figure(10); clf; imagesc(fly2Dmaps.imapx(1,:),fly2Dmaps.imapy(:,1),tilt_tot.*mask); colorbar; title('tilt total');axis image;
                set(gcf,'Position',fig_position);
                
                figure(11);clf;imagesc(ROI_radial);axis image;title('Radial component in ccd');colorbar;colormap jet;
                figure(12);clf;imagesc(ROI_azi);axis image;title('Azhimutal component in ccd');colorbar;colormap jet;
            end
            %}
            
        end
        
        function [distr_struct,mask_struct] = computeStrainOrTiltContours(dat,struct_centroidShift,field,contour_values_up,contour_values_down)
            
            map2D_SumInt = dat.map2D_SumInt.*dat.mask;
            dat.map2D_SumInt = map2D_SumInt;
            for kk = 1:numel(contour_values_up)
                [mask_struct(kk).mask_up] = ND_data_processing.calculateMask(dat,contour_values_up(kk));
                [mask_struct(kk).mask_down] = ND_data_processing.calculateMask(dat,contour_values_down(kk));

                mask_ring = zeros(size(mask_struct(kk).mask_up));
                mask_ring (mask_struct(kk).mask_down>0) = 1;
                mask_ring(mask_struct(kk).mask_up>0) = 0; 
                
                mask_struct(kk).mask = mask_ring;
                
                distrib = struct_centroidShift.(field).*mask_ring;
               
                weight = mask_ring.*map2D_SumInt./sum(sum(map2D_SumInt.*mask_ring));
                
                distrib_weight = struct_centroidShift.(field).*mask_ring.*map2D_SumInt./sum(sum(map2D_SumInt.*mask_ring));
                
                sigma(kk) = std(distrib(:),weight(:),'omitnan');
                
                median_field(kk) = sum(sum(distrib_weight,'omitnan'),'omitnan');%/sum(sum(dat.map2D_SumInt.*mask_ring));
               
                distr_struct.hist(kk).distrib = distrib;
                distr_struct.hist(kk).distrib_weight = distrib_weight;
                distr_struct.hist(kk).weight = weight;
                %strain_weight = struct_centroidShift.strain.*mask_struct(kk).mask.*dat.map2D_SumInt;
                %distr_struct.strain(kk) = sum(sum(strain_weight,'omitnan'),'omitnan')/sum(sum(dat.map2D_SumInt.*mask_struct(kk).mask));
            end
            
             distr_struct.distr = median_field;
                
             distr_struct.sigma = sigma;
            
            
        end
       
        
        
        function f_struct = qmatrix_MAR(Ekev,Tth,Gam,xdet,ydet,Rdet, Ndet,pixsize,plotflag)
            % function to create momentum transfer matrix for MAR detector
            % Rdet cylindrical, ydet vertically up, xdet horizontal in detector plane - scale in mm
            % used in finding pixel track out as a function of energy, use syntax
            % tempchi = (qq2(:,:,1)-qq1(400,400,1)).^2+(qq2(:,:,2)-qq1(400,400,2)).^2+(qq2(:,:,3)-qq1(400,400,3)).^2;
            % [temprow tempcol]=find(tempchi==min(min(tempchi)));
            
            %Ndet = round(2048/rebinrate); %MAR
            %Ndet = round(1690/rebin_rate);  %CSPAD
            %Ndet = 516;
            
            Ncen = round(Ndet/2);
            imaxis = (1-Ncen):(Ndet-Ncen);
            
            %pxsz = 165/Ndet; %Pixel size MAR165 in mm
            %pxsz = 186/Ndet; %Pixel size CSPAD in mm
            %pxsz = 0.055;
             kb = 2*pi*Ekev/12.39842; % beam momentum 1/A
          
            
            %Detector center and pixel plane vectors
            Detcen = [Rdet*sind(Tth)*cosd(Gam) Rdet*sind(Gam) Rdet*cosd(Tth)*cosd(Gam)];
            Detunit = Detcen./(sqrt(dot(Detcen,Detcen)));
            jvec = cross([0 1 0],Detunit); % horizontal pixels unit vector NSLSII (high pixel number, low tth)
            jvec = jvec./(sqrt(dot(jvec,jvec))); %make sure it is a unit vector for nonzero gamma
            ivec = cross(jvec,Detunit); % vertical pixels unit vector for NSLSII (low pixel number, more up)
            Detcen = Detcen - ydet*ivec + xdet*jvec; %ydet moves it up, xdet moves it outboard
            
            kfmat=zeros(Ndet,Ndet,4);
            qmat=zeros(Ndet,Ndet,4);
            
            %qscale = 16;
            %smallq=zeros(Ndet/qscale,Ndet/qscale,4);smallx=1:Ndet/qscale;
            
            % Assign magnitude and vector q to active detector pixels within q and gamma range
            %display('assigning q vectors...');
            %Generate pixel arrays (um)
            pix_x = repmat((imaxis.*(pixsize)),Ndet,1);
            pix_y = repmat((imaxis.*(pixsize))',1,Ndet);
            %Calculate real space pixel positions in mm as a vector sum
            kfmat(:,:,1) = pix_x(:,:).*jvec(1)+pix_y(:,:).*ivec(1)+Detcen(1);
            kfmat(:,:,2) = pix_x(:,:).*jvec(2)+pix_y(:,:).*ivec(2)+Detcen(2);
            kfmat(:,:,3) = pix_x(:,:).*jvec(3)+pix_y(:,:).*ivec(3)+Detcen(3);
            kfmat(:,:,4) = sqrt(kfmat(:,:,1).^2+kfmat(:,:,2).^2+kfmat(:,:,3).^2);
            %Assign kfinal unit vector to each pixel based on position
            kfmat(:,:,1) = kfmat(:,:,1)./kfmat(:,:,4);
            kfmat(:,:,2) = kfmat(:,:,2)./kfmat(:,:,4);
            kfmat(:,:,3) = kfmat(:,:,3)./kfmat(:,:,4);
            %Calculate q vectors
            qmat(:,:,1) = kb*(kfmat(:,:,1)-0);
            qmat(:,:,2) = kb*(kfmat(:,:,2)-0);
            qmat(:,:,3) = kb*(kfmat(:,:,3)-1);
            qmat(:,:,4) = sqrt(qmat(:,:,1).^2+qmat(:,:,2).^2+qmat(:,:,3).^2);
            
            anglemat(:,:,1) = acosd(kfmat(:,:,3)); % radial
            anglemat(:,:,2) = atand(-kfmat(:,:,2)./kfmat(:,:,1)); %azimutal
            %anglemat(:,:,3) = 180 - acosd(qmat(:,:,4)/(2*kb)); %delta two theta
            
            %smallq(smallx,smallx,:)=qmat(smallx*qscale,smallx*qscale,:);
            
            if (plotflag==1)
                scrsz = get(0,'ScreenSize');  % Plot size parameters
                imd = scrsz(3)/6;
                imb = scrsz(4)-2*imd;
                %%{
                figure(11);set(11,'Position',[imd imb imd imd]);clf;imagesc(qmat(:,:,4));colormap hot;axis square;colorbar;
                title('Magnitude q per pixel');
                figure(12);set(12,'Position',[2*imd imb imd imd]);clf;imagesc(qmat(:,:,1));colormap hot;axis square;colorbar;
                title('Magnitude qx per pixel');
                figure(13);set(13,'Position',[3*imd imb imd imd]);clf;imagesc(qmat(:,:,2));colormap hot;axis square;colorbar;
                title('Magnitude qy per pixel');
                figure(14);set(14,'Position',[4*imd imb imd imd]);clf;imagesc(qmat(:,:,3));colormap hot;axis square;colorbar;
                title('Magnitude qz per pixel');
                figure(15);set(15,'Position',[5*imd imb imd imd]);clf;imagesc(anglemat(:,:,1));colormap hot;axis square;colorbar;
                title('Two Theta polar angle (deg) per pixel');
                figure(16);set(16,'Position',[6*imd imb imd imd]);clf;imagesc(anglemat(:,:,2));colormap hot;axis square;colorbar;
                title('Tilt (Eta) from horizontal plane per pixel');
                %}
                %    figure(18);clf;imagesc(qmat(:,:,4));colormap hot; axis square tight;
                %    title('Magnitude q per pixel');
                %   figure(40);clf;surf(qmat(:,:,1),qmat(:,:,2),qmat(:,:,3));
                
                %    figure(40);clf;surf(smallq(:,:,1),smallq(:,:,2),smallq(:,:,3));
                %    title('Lab frame momentum transfer (1/A)');
                %    hold on;scatter3(0,0,0,400,'blue','LineWidth',5);
            end
            
            f_struct.anglemat = anglemat;
            f_struct.qmat = qmat;
            f_struct.jvec = jvec;
            f_struct.ivec = ivec;
            f_struct.Detcen = Detcen;
            
        end
        
        function Rtest = calculateCorrCoef(dat,mask,filenames,plotflag)
            
            dat_align = load(filenames{1},'dat');%load(['results/data_scan_align' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');
            dat_alignXRF0 = load(filenames{2},'dat');%load(['results/data_scan_alignXRF0' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');
            dat_nonalign = load(filenames{3},'dat');;%load(['results/data_scan_nonalign' num2str(lst(1)) '_' num2str(lst(end)) '.mat'],'dat');
            
            dat_align.dat.thetalist = dat.thetalist;
            dat_alignXRF0.dat.thetalist = dat.thetalist;
            dat_nonalign.dat.thetalist = dat.thetalist;
            
            [centroid_struct_align] = ND_analysis.computeCentroidsRockCurve(dat_align.dat,'mask',mask,'do_plot',0);
            [centroid_struct_alignXRF0] = ND_analysis.computeCentroidsRockCurve(dat_alignXRF0.dat,'mask',mask,'do_plot',0);
            [centroid_struct_nonalign] = ND_analysis.computeCentroidsRockCurve(dat_nonalign.dat,'mask',mask,'do_plot',0);
            
            Xcentroids_corr(:,1) = centroid_struct_align.Xcentroids(:);
            Xcentroids_corr(:,2) = centroid_struct_alignXRF0.Xcentroids(:);
            Xcentroids_corr(:,3) = centroid_struct_nonalign.Xcentroids(:);
            
            Rtest.Xcentroids_corrcoef = corrcoef(Xcentroids_corr,'rows','pairwise');
            
            Ycentroids_corr(:,1) = centroid_struct_align.Ycentroids(:);
            Ycentroids_corr(:,2) = centroid_struct_alignXRF0.Ycentroids(:);
            Ycentroids_corr(:,3) = centroid_struct_nonalign.Ycentroids(:);
            
            Rtest.Ycentroids_corrcoef = corrcoef(Ycentroids_corr,'rows','pairwise');
            
            Thcentroids_corr(:,1) = centroid_struct_align.Thcentroids(:);
            Thcentroids_corr(:,2) = centroid_struct_alignXRF0.Thcentroids(:);
            Thcentroids_corr(:,3) = centroid_struct_nonalign.Thcentroids(:);
            
            Rtest.Thcentroids_corrcoef = corrcoef(Thcentroids_corr,'rows','pairwise');
            
            if plotflag
                
                figure(1);
                imagesc(Rtest.Xcentroids_corrcoef);
                colorbar;
                axis square;
                title('Correlation coefficients for Ycentroids')
                
                figure(2);
                imagesc(Rtest.Ycentroids_corrcoef);
                colorbar;
                axis square;
                title('Correlation coefficients for Ycentroids');
               
                figure(3);
                imagesc(Rtest.Thcentroids_corrcoef);
                colorbar;
                axis square;
                title('Correlation coefficients for Thcentroids')
                
            end
            
        end
        
    end
    
end
