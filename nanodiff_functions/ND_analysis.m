classdef ND_analysis
    % This library contains all the functions to analyse the
    % nanodiffraction pattern
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function diff_data = computeCentroidsfromDet(merlimgs,ROIinteg)
           
            numimgs = size(merlimgs,3);
            diff_data = zeros(numimgs,4);
            
            for ii = 1:size(merlimgs,3)
            
                ccd = merlimgs(:,:,ii);
                
                diff_data(ii,1) = sum(sum(ccd));
                
                if(~isempty(ROIinteg))
                    diff_data(ii,4) = sum(sum(ccd(ROIinteg(1):ROIinteg(2), ROIinteg(3):ROIinteg(4))));
                end
                
               
                    line1=sum(ccd,1);  % vertical
                    line2=sum(ccd,2);  % horizontal
                    for kk=1:size(line1,2)
                        diff_data(ii,2)=diff_data(ii,2)+kk*line1(kk)/diff_data(ii,1);
                    end
                    for kk=1:size(line2,1)
                        diff_data(ii,3)=diff_data(ii,3)+kk*line2(kk)/diff_data(ii,1);
                    end
                
            end
        end
        
        
        function fly2Dmaps = computeCentroids_rockCurve(fly2Dmaps)
           
            tempim = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            tempxcen = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            tempycen = zeros(numel(fly2Dmaps.ii),numel(fly2Dmaps.ii(1).jj));
            
            for kk = 1:numel(fly2Dmaps.ii)
                for ll = 1:numel(fly2Dmaps.ii(1).jj)
                    tempim(kk,ll) = fly2Dmaps.ii(kk).jj(ll).SumInt;
                    imgin =fly2Dmaps.ii(kk).jj(ll).im;
                    
                    line1=sum(imgin,1);  % vertical
                    line2=sum(imgin,2);  % horizontal
                    sumt = sum(sum(imgin));
                    
                    if fly2Dmaps.mask_rock(kk,ll) == 1
                        if fly2Dmaps.mask(kk,ll) == 1
                            if sumt==0
                                tempycen(kk,ll)= 0;
                                tempxcen(kk,ll)= 0;
                            else
                                
                                for mm=1:size(line1,2)
                                    tempycen(kk,ll)=tempycen(kk,ll)+mm*line1(mm)/sumt;
                                end
                                for mm=1:size(line2,1)
                                    tempxcen(kk,ll)=tempxcen(kk,ll)+mm*line2(mm)/sumt;
                                end
                            end
                        end
                    end
                    
                end
            end
            
            fly2Dmaps.Xcentroids = tempxcen;
            fly2Dmaps.Ycentroids = tempycen;

        end
        
        
        
        function [struct_centroidShift] = computeCentroidShift(dat1,ROIxstart,ROIxsize,ROIystart,ROIysize,plotflag)
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%
            %%{
            %CIGS_diffraction_detector;
            angles = ND_analysis.qmatrix_MAR(Ekev,del,gam,0,0,detdist/100,Ndet,pixsize,0);
            
            strain_D = angles(:,:,1);
            tilt_D = angles(:,:,2);
            ROI_radial = strain_D(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
            ROI_azi = tilt_D(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
            
            
            xcen = round(dat1.Xcentroids(:));
            ycen = round(dat1.Ycentroids(:));
            radial_proj = zeros(size(xcen));
            azi_proj = zeros(size(xcen));
            
            meanval_x = [];  %calculate the mean centroid value, excluding fluo
            meanval_y = [];
            XRFtemp = dat1.scan(1).XRF(:).*dat1.mask(:).*dat1.mask_rock(:);
            mask_array = dat1.mask(:).*dat1.mask_rock(:);
            minXRF = max(XRFtemp)*0.1;
            for ii = 1:numel(xcen)
                if xcen(ii) ~=0
                        radial_proj(ii) = ROI_radial(ycen(ii),xcen(ii));
                        azi_proj(ii) = ROI_azi(ycen(ii),xcen(ii));
                        
                        if dat1.mask(ii)
                        %if XRFtemp(ii)>minXRF
                            meanval_x = [meanval_x, radial_proj(ii)];
                            meanval_y = [meanval_y, azi_proj(ii)];
                        end
                end
            end
            
            meanval_radial = mean(meanval_x);
            meanval_azi = mean(meanval_y);
            
        
            
            radial_proj = reshape(radial_proj, size(dat1.Xcentroids));
            azi_proj = reshape(azi_proj, size(dat1.Xcentroids));
            xcen = reshape(xcen, size(dat1.Xcentroids));
            ycen = reshape(ycen, size(dat1.Xcentroids));
            
            
          
            
            %{
figure(1); clf; imagesc(radial_proj); axis image; colorbar; title('radial projected centroid');
figure(2); clf; imagesc(azi_proj); axis image; colorbar; title('tangential projected centroid');
            %}
            dspace = lam/(2*sind(meanval_radial/2));
            
            tilt = azi_proj - meanval_azi;
            %strain = -cotd(meanval_radial/2).*(meanval_radial- radial_proj)/2;
            strain = -cotd(meanval_radial/2).*(-meanval_radial+radial_proj)/2;
            
            struct_centroidShift.angles = angles;
            struct_centroidShift.dspace = dspace;
            struct_centroidShift.tilt = tilt;
            struct_centroidShift.strain = strain;
            struct_centroidShift.radial_proj = radial_proj;
            struct_centroidShift.azi_proj = azi_proj;
            struct_centroidShift.xcen = xcen;
            struct_centroidShift.ycen = ycen;
            
            %%{
            if plotflag
                
                ND_display_data.display2Dmap(radial_proj,'figNum',1,'figTitle','radial projected centroid');
                ND_display_data.display2Dmap(azi_proj,'figNum',1,'figTitle','tangential projected centroid');

                %{
                figure(3); clf;
                subplot(3,1,1); s=imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),strain.*msk); axis image; colorbar; title('strain');
                subplot(3,1,2); imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),tilt.*msk); axis image; colorbar; title('tilt');
                tmp  = dat1.scan(1).PC;
                subplot(3,1,3); p3=imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),dat1.scan(2).PC.*msk); axis image; colorbar; title('photocurrent'); caxis([.8*median(tmp(:)) 1.2*median(tmp(:))])
                
                alpha(p3,0.5)
                %}
                figure(1); clf; imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),xcen); colorbar; title('xcen');axis image;
                figure(2); clf; imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),ycen); colorbar; title('ycen');axis image;
                figure(3); clf; imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),radial_proj); colorbar; title('radial cen');axis image;
                figure(4); clf; imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),azi_proj); colorbar; title('tangential cen');axis image;
                figure(5); clf; imagesc(strain); colorbar; title('strain-v2');axis image;
                figure(6); clf; imagesc(tilt); colorbar; title('tilt-v2');axis image;
         
            end
            %}
            
        end
        
        function [strain_struct] = calculateStrain(dat1,ROIxstart,ROIxsize,ROIystart,ROIysize,plotflag)
            % Convert to strain
            eval('Init_parameters');
            
            
            ROICenter = [ ...
                ROIxstart + ROIxsize / 2 ...
                Ndet - ( ROIystart + ROIysize / 2 ) ...
                ];
            
            detectorCenter = [ 1 1 ] * ( Ndet/2 );
            ROIOffset = ROICenter - detectorCenter;
            
            a = detdist * kf;
            b = detdist * cosd( twoTheta ) * ki;
            x_strain = b - a;
            x_strain = x_strain / norm( x_strain );
            
            R = [ ...
                0 -1 0 ; ...
                1 0 0 ; ...
                0 0 1 ];
            x_tilt = R * x_strain;
            trans = [ x_strain x_tilt ];
            
            XcentroidsCorrected = ROIOffset(1) + dat1.Xcentroids;
            YcentroidsCorrected = ROIOffset(2) + dat1.Ycentroids;
            
            
            centroids = ( ROIOffset(1) + dat1.Xcentroids(:) ) * [ 1 0 0 ] + ( ROIOffset(2) + dat1.Ycentroids(:) ) * [ 0 1 0 ];
            centroidsTransformed = centroids * trans;
            
            centroids_strain = reshape(centroidsTransformed(:,1),size(dat1.Xcentroids));
            centroids_tilt = reshape(centroidsTransformed(:,2),size(dat1.Ycentroids));
            
            tilt = centroids_tilt.*degperpix;
            
            strain = centroids_strain.*degperpix;
            rel_th = mean(strain(:));
            strain = -cotd(twoTheta/2).*(rel_th - strain);
            
            strain_struct.strain = strain;
            strain_struct.rel_th = rel_th;
            strain_struct.tilt = tilt;
            strain_struct.XcentroidsCorrected = XcentroidsCorrected;
            strain_struct.YcentroidsCorrected = YcentroidsCorrected;
            %}
            %%{
            %ROI_radial = Angs_radial(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
            %ROI_azi = Angs_azi(ROIystart:ROIystart+ROIysize-1,ROIxstart:ROIxstart+ROIxsize-1);
            %%{
            if plotflag
                figure(9); s=imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),strain); colorbar; title('strain-v1');
                figure(10); imagesc(dat1.xaxis(1,:),dat1.yaxis(:,1),tilt); colorbar; title('tilt-v1');
                %{
                figure(1); clf; imagesc(Angs_radial); colorbar; ca = caxis; axis image;
                figure(2); clf; imagesc(ROI_radial); colorbar; caxis(ca); axis image;
                figure(3); clf; imagesc(Angs_azi); colorbar; ca = caxis; axis image;
                figure(4); clf; imagesc(ROI_azi); colorbar; caxis(ca); axis image; colormap hot;
                %}
            end
            %}
            
            
        end
        
        function f = qmatrix_MAR(Ekev,Tth,Gam,xdet,ydet,Rdet, Ndet,pixsize,plotflag)
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
            
            anglemat(:,:,1) = acosd(kfmat(:,:,3)); % two theta
            anglemat(:,:,2) = atand(-kfmat(:,:,2)./kfmat(:,:,1)); %eta
            
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
            
            f = anglemat;
            
        end
       
        
        function [mask,dat1] = calculateMask(dat1,percent)
            
            mask_rock = ones(numel(dat1.ii),numel(dat1.ii(1).jj));
            for kk = 1:numel(dat1.ii)
                for ll = 1:numel(dat1.ii(kk).jj)
                    map2D_SumInt(kk,ll) = dat1.ii(kk).jj(ll).SumInt;
                    [max_val,max_rock_curve(kk,ll)] = max(dat1.ii(kk).jj(ll).intensity);
                    
                    if max_rock_curve(kk,ll) == numel(dat1.ii(kk).jj(ll).intensity)
                        mask_rock(kk,ll) = 0 ;
                    end
                    
                end
            end
            
            % mask for intensities
            mask = map2D_SumInt > percent*max(max(map2D_SumInt));
            
            % mask for rocking curve
            
            dat1.mask = mask;
            dat1.mask_rock = mask_rock;
        end
        
    end
    
end