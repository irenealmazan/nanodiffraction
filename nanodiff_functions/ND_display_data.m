
classdef ND_display_data
    % This library contains all the functions to display the data read by
    % loadscan
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function hfig = showmerlin_function(scandata,pResults,pixx,pixy,fignum)
            
            figure(fignum);
            clf reset;
            if(pResults.outerpts>1)
                hfig = imagesc(scandata(1,:,3),scandata(:,1,2),scandata(:,:,1));%./scandata(:,:,4)*mean(mean(scandata(:,:,4))));
                colormap jet; shading interp;colorbar; set(gca, 'YDir', 'reverse');axis image
            else
                %hfig=plot(scandata(1,:,3),scandata(1,:,1)./scandata(1,:,4)*mean(scandata(1,:,4)));
                hfig = plot(scandata(1,:,3),scandata(1,:,1));
            end
            
            %set(gcf, 'Position', [800 1100 600 400]);
            %set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
            %set(gcf, 'Color' ,'w');
            %set(gca, 'FontSize',8);
            title(['Scan: ' num2str(pResults.scanid) '  Detector: ' num2str(pResults.detchan)], 'Interpreter', 'none', 'FontSize', 10);
            
            
            figure(fignum+10);
            clf reset;
            if(pResults.outerpts>1)
                %hfig = imagesc(scandata(1,:,3),scandata(:,1,2),scandata(:,:,5)./scandata(:,:,4)*mean(mean(scandata(:,:,4))));
                hfig = imagesc(scandata(1,:,3),scandata(:,1,2),scandata(:,:,5));
                colormap jet; shading interp;colorbar;set(gca, 'YDir', 'reverse');axis image;
            else
                %hfig=plot(scandata(1,:,3),scandata(1,:,1)./scandata(1,:,5)*mean(scandata(1,:,4)));
                hfig=plot(scandata(1,:,3),scandata(1,:,5));
            end
            
            %set(gcf, 'Position', [1400 1100 600 400]);
            %set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
            %set(gcf, 'Color' ,'w');
            %set(gca, 'FontSize',8);
            title(['Scan: ' num2str(pResults.scanid) '  Detector: Integrated diffraction'], 'Interpreter', 'none', 'FontSize', 10);
            
            
            if(~isempty(pResults.ROIinteg))
                figure(fignum+11);
                imagesc(scandata(:,:,9));
                axis image;
            end
            
            
            
            
            
        end
        
        
        function hfig = show_centroid(scandata,pResults,pixx,pixy,fignum)
            
            
            figure(fignum+3);
            clf reset;
            if(pResults.outerpts>1)
                hfig = imagesc(scandata(1,:,3),scandata(:,1,2),scandata(:,:,6));
                colormap hot; shading interp;colorbar;set(gca, 'YDir', 'reverse'); axis image;
            else
                hfig=plot(scandata(1,:,3),scandata(1,:,6));
            end
            %set(gcf, 'Position', [2000 1100 600 400]);
            %set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
            %set(gcf, 'Color' ,'w');
            %set(gca, 'FontSize',8);
            title(['Scan: ' num2str(pResults.scanid) '  Detector: Y Centroid'], 'Interpreter', 'none', 'FontSize', 10);
            
            pass2click.imrows = pResults.outerpts;
            pass2click.imcols = pResults.innerpts;
            pass2click.xaxis = [min(scandata(1,:,3)) max(scandata(1,:,3)) size(scandata(1,:,3),2)];
            pass2click.yaxis = [min(scandata(:,1,2)) max(scandata(:,1,2)) size(scandata(:,1,2),1)];
            pass2click.merlfilepath = [pResults.datapath '/scan_' num2str(pResults.scanid) '.h5'];
            pass2click.merlhdfpath = pResults.datapath;
            pass2click.pixx = pixx;
            pass2click.pixy = pixy;
            pass2click.ccdnums = scandata(:,:,8);
            
            
            
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
            set(dcm_obj, 'UpdateFcn', @click4ccd_nsls);
            
            figure(fignum+4);
            clf reset;
            if(pResults.outerpts>1)
                hfig = imagesc(scandata(1,:,3),scandata(:,1,2),scandata(:,:,7));
                colormap hot; shading interp;colorbar;set(gca, 'YDir', 'reverse'); axis image;
            else
                hfig=plot(scandata(1,:,3),scandata(1,:,7));
            end
            %set(gcf, 'Position', [2600 1100 600 400]);
            %set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
            %set(gcf, 'Color' ,'w');
            %set(gca, 'FontSize',8);
            title(['Scan: ' num2str(pResults.scanid) '  Detector: X Centroid'], 'Interpreter', 'none', 'FontSize', 10);
            
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
            set(dcm_obj, 'UpdateFcn', @click4ccd_nsls);
            
            
        end
        
        function hfig = displayRockCurveLine(scanid,fly2Dmaps,pix_list,varargin)
            p = inputParser;
            
            addRequired(p,'scanid');
            addRequired(p,'fly2DMaps');
            addRequired(p,'pix_list');
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'titleXray',1,@isnumeric);
            addParameter(p,'titleFluo',1,@isnumeric);
            addParameter(p,'titleXBIC',1,@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,scanid,fly2Dmaps,pix_list,varargin{:});
            
            numRows_unique = unique(p.Results.pix_list(:,1));
            numCols_unique = unique(p.Results.pix_list(:,2));
            
            num_supblots = numel(numRows_unique);
            
            
            
            for kk = 1:numel(p.Results.pix_list(:,1))
                for ll = 1:numel(p.Results.pix_list(:,2))
                    for jj = 1:numel(scanid)
                        
                        diff_map(kk,jj,ll) = fly2Dmaps.ii(p.Results.pix_list(kk,1)).jj(p.Results.pix_list(ll,2)).intensity(jj);
                        thetalist(jj) = fly2Dmaps.scan(jj).theta;
                    end
                end
                
                
            end
            
            figure(p.Results.figNum);
            clf;
            
            if numel(numRows_unique) == 1
                imagesc(p.Results.pix_list(:,2),thetalist,squeeze(diff_map(1,:,:)));
                title(['Rocking curve at row = ' num2str(p.Results.pix_list(1,1))]);
                colormap jet;
            elseif numel(numCols_unique) == 1
                imagesc(p.Results.pix_list(:,1),thetalist,squeeze(diff_map(:,1,:)));
                title(['Rocking curve at column = ' num2str(p.Results.pix_list(1,2))]);
                colormap jet;
                %axis image;
            end
            
        end
        
        function [thetalist,rock_curve] = displayRockCurve(fly2Dmaps,varargin)
             p = inputParser;
            
            
            addRequired(p,'fly2Dmaps');
            %addParameter(p,'list',[1:numel(fly2Dmaps.scan)],@isnumeric);
            addParameter(p,'title_spec','scanid');
            addParameter(p,'do_mask',0,@isnumeric);
            addParameter(p,'figNum',1,@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,fly2Dmaps,varargin{:});
            
            %fly2Dmaps = p.Results.fly2Dmaps;
            
            for jj = 1:numel(fly2Dmaps.scan)
                for kk = 1:numel(fly2Dmaps.ii)
                    for ll = 1:numel(fly2Dmaps.ii(kk).jj)
                        if p.Results.do_mask == 0
                            diff_map(kk,ll) = fly2Dmaps.ii(kk).jj(ll).intensity(jj);
                        else
                            diff_map(kk,ll) = fly2Dmaps.mask(kk,ll)*fly2Dmaps.mask_rock(kk,ll)*fly2Dmaps.ii(kk).jj(ll).intensity(jj);
                        end
                    end
                end
                
                rock_curve(jj) =  sum(sum(diff_map));
                thetalist(jj) = fly2Dmaps.scan(jj).theta;
                
            end
            
            figure(p.Results.figNum);
            plot(thetalist,rock_curve,'LineWidth',3.0);
            xlabel('Angle (deg)');
            title('Rocking curve');
            
        end
        
        function [rock_curve,thetalist] = displayRockCurveMaps(scanid,fly2Dmaps,varargin)
            p = inputParser;
            
            addRequired(p,'scanid');
            addRequired(p,'fly2DMaps');
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'do_mask',0,@isnumeric);
            addParameter(p,'size_figure',[100 100 1000 800],@isnumeric);
            addParameter(p,'window',[1,numel(fly2Dmaps.ii),1,numel(fly2Dmaps.ii(1).jj)],@isnumeric);
            addParameter(p,'plot_XBIC',1,@isnumeric);
            addParameter(p,'titleXray','Diffr. Int');
            addParameter(p,'titleFluo','XRF');
            addParameter(p,'titleXBIC','Photo current');
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,scanid,fly2Dmaps,varargin{:});
            
            if p.Results.plot_XBIC
                numRows = 3;
            else
                numRows = 2;
            end
            
            numCols = 6;
            
            numFigs = ceil(numel(scanid)/numCols);
            
            counter = 1;
            
            for jj = 1:numel(scanid)
                for kk = 1:numel(fly2Dmaps.ii)
                    for ll = 1:numel(fly2Dmaps.ii(kk).jj)
                        if p.Results.do_mask == 0
                            diff_map(kk,ll) = fly2Dmaps.ii(kk).jj(ll).intensity(jj);
                            PC_map(kk,ll) = fly2Dmaps.scan(jj).PC(kk,ll);
                        else
                            diff_map(kk,ll) = fly2Dmaps.mask(kk,ll)*fly2Dmaps.mask_rock(kk,ll)*fly2Dmaps.ii(kk).jj(ll).intensity(jj);
                            PC_map(kk,ll) = fly2Dmaps.scan(jj).PC(kk,ll)*fly2Dmaps.mask_rock(kk,ll)*fly2Dmaps.mask(kk,ll);
                        end
                    end
                end
                
                rock_curve(jj) =  sum(sum(diff_map));
                thetalist(jj) = fly2Dmaps.scan(jj).theta;
                
                if mod(jj-1,numCols) == 0
                    figure(p.Results.figNum+counter);
                    clf;
                    set(gcf,'Name',['Rock curve' num2str(counter) '/' num2str(numFigs)]);
                    set(gcf, 'Position', p.Results.size_figure);
                    set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
                    set(gcf, 'Color' ,'w');
                    set(gca, 'FontSize',8);
                    
                    counter = counter + 1;
                    
                end
                
                %subplot(numel(scanid),4,jj);
                subplot(numRows,numCols,mod(jj-1,numCols)+1);
                h_struct(jj).hfig = imagesc('XData',fly2Dmaps.xsuperGrid_s(1,p.Results.window(3):p.Results.window(4)),'YData',fly2Dmaps.ysuperGrid_s(1,p.Results.window(1):p.Results.window(2)),'CData',diff_map(p.Results.window(1):p.Results.window(2),p.Results.window(3):p.Results.window(4)));
                title({[p.Results.titleXray ' # ' num2str(scanid(jj))],[ ' th = ' num2str(thetalist(jj))]});
                colormap jet;
                %axis image;
                xlim([-5 50]);
                colorbar;
                %{
                if jj == 1
                    colorbar;
                end
                %}
                
                %subplot(numel(scanid),4,2*4+jj);
                subplot(numRows,numCols,(numCols)+mod(jj-1,numCols)+1);
                h_struct(jj).hfig = imagesc('XData',fly2Dmaps.xsuperGrid_s(1,p.Results.window(3):p.Results.window(4)),'YData',fly2Dmaps.ysuperGrid_s(1,p.Results.window(1):p.Results.window(2)),'CData',fly2Dmaps.scan(jj).XRF);
                title({[p.Results.titleFluo  ' # ' num2str(scanid(jj))],[ ' th = ' num2str(thetalist(jj))]});
                colormap jet;
                xlim([-5 50]);
                %axis image;
                colorbar;
                %{
                if jj == 1
                    colorbar;
                end
                %}
                
                if p.Results.plot_XBIC
                    subplot(numRows,numCols,(numCols)*2+mod(jj-1,numCols)+1);
                    h_struct(jj).hfig = imagesc('XData',fly2Dmaps.xsuperGrid_s(1,p.Results.window(3):p.Results.window(4)),'YData',fly2Dmaps.ysuperGrid_s(1,p.Results.window(1):p.Results.window(2)),'CData',PC_map);
                    title({[p.Results.titleXBIC  ' # ' num2str(scanid(jj))] , [' th = ' num2str(thetalist(jj))]});
                    colormap jet;
                    xlim([-5 50]);
                    %axis image;
                    colorbar;
                end
                %{
                if jj == 1
                    colorbar;
                end
                %}
                
            end
            
            
            
            figure(numFigs+500);
            plot(thetalist,rock_curve,'LineWidth',3.0);
            xlabel('Angle (deg)');
            title('Rocking curve');
            
        end
        
        function displayRockCurveShift(fly2Dmaps,varargin)
            p = inputParser;
            
            
            addRequired(p,'fly2Dmaps');
            addParameter(p,'list',[1:numel(fly2Dmaps.scan)],@isnumeric);
            addParameter(p,'title_spec','scanid');
            addParameter(p,'figNum',1,@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,fly2Dmaps,varargin{:});
            
            %fly2Dmaps = p.Results.fly2Dmaps;
            
            for jj = 1:numel(fly2Dmaps.scan)
                xshiftlist(jj) = fly2Dmaps.scan(jj).xshifts;
                yshiftlist(jj) = fly2Dmaps.scan(jj).yshifts;
                thetalist(jj) = fly2Dmaps.scan(jj).theta;
            end
            
            figure(p.Results.figNum);
            clf;
            
            subplot(2,1,1);
            %plot(thetalist,xshiftlist);
            %title('Xshift vs theta');
            plot(p.Results.list,xshiftlist);
            title(['Xshift vs ' p.Results.title_spec] );
            
            subplot(2,1,2);
            %plot(thetalist,yshiftlist);
            %title('Yshift vs theta');
            plot(p.Results.list,yshiftlist);
            title(['Yshift vs ' p.Results.title_spec]);
            
        end
        
        function displayXRFmapsFrames(XRF_align_struct,varargin)
            p = inputParser;
            
            addRequired(p,'XRF_align_struct');
            addParameter(p,'show_xmaps',0,@isnumeric);
            addParameter(p,'theta_shift_deg',0,@isnumeric);
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'size_figure',[6 437 1203 314],@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,XRF_align_struct,varargin{:});
            
            %XRF_align_struct = p.Results.XRF_align_struct;
            
            numCols = 1;
            numRows = 1;
            
            numFigs = ceil(numel(XRF_align_struct.scan)/numCols);
            
            counter = 1;
            
            for jj = 1:numel(XRF_align_struct.scan)
                
                
                
                if mod(jj-1,numCols) == 0
                    figure(p.Results.figNum+counter);
                    clf;
                    set(gcf,'Name',['Rock curve' num2str(counter) '/' num2str(numFigs)]);
                    set(gcf, 'Position', p.Results.size_figure);
                    set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
                    set(gcf, 'Color' ,'w');
                    set(gca, 'FontSize',8);
                    
                    counter = counter + 1;
                    
                end
                
                %subplot(numel(scanid),4,jj);
                %subplot(numRows,numCols,mod(jj-1,numCols)+1);
                subplot(211);
                h_struct(jj).hfig = imagesc('XData',XRF_align_struct.scan(jj).x_pos,'YData',XRF_align_struct.scan(jj).y_pos,'CData',XRF_align_struct.scan(jj).XRF);
                title({['Map. int. # ' num2str(XRF_align_struct.scan(jj).scanid)],...
                    [ ' th = ' num2str(XRF_align_struct.scan(jj).theta)]});
                colormap jet;
                axis image;
                colorbar;
                %caxis([-7 4]);
                clim_1 = get(gca,'CLim');
                
                subplot(212);
                %h_struct(jj).hfig = imagesc(XRF_align_struct.scan(jj).XRF);
               h_struct(jj).hfig = imagesc('XData',XRF_align_struct.scan(jj).x_pos_supergrid,'YData',XRF_align_struct.scan(jj).y_pos_supergrid,'CData',XRF_align_struct.scan(jj).XRF_supergrid);
                title({['Map. int. # ' num2str(XRF_align_struct.scan(jj).scanid)],...
                    [ ' th = ' num2str(XRF_align_struct.scan(jj).theta)]});
                colormap jet;
                %axis image;
                colorbar;
                caxis([clim_1(1) clim_1(2)]);
                
                
            end
            
            
        end
        
        function displayAlignedXRFmaps(XRF_align_struct,field_to_show,varargin)
            p = inputParser;
            
            addRequired(p,'XRF_align_struct');
            addRequired(p,'field_to_show');
            addParameter(p,'show_xmaps',0,@isnumeric);
            addParameter(p,'numCols',1,@isnumeric);
            addParameter(p,'numRows',1,@isnumeric);
            addParameter(p,'theta_shift_deg',0,@isnumeric);
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'size_figure',[6 437 1203 314],@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,XRF_align_struct,field_to_show,varargin{:});
            
            %XRF_align_struct = p.Results.XRF_align_struct;
            
            numCols = p.Results.numCols;
            numRows = p.Results.numRows;
            
            numFigs = ceil(numel(XRF_align_struct.scan)/numCols);
            
            counter = 1;
            
            for jj = 1:numel(XRF_align_struct.scan)
                
                
                
                if mod(jj-1,numCols) == 0
                    figure(p.Results.figNum+counter);
                    clf;
                    set(gcf,'Name',['Rock curve' num2str(counter) '/' num2str(numFigs)]);
                    set(gcf, 'Position', p.Results.size_figure);
                    set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
                    set(gcf, 'Color' ,'w');
                    set(gca, 'FontSize',8);
                    
                    counter = counter + 1;
                    
                end
                
                %subplot(numel(scanid),4,jj);
                subplot(numRows,numCols,mod(jj-1,numCols)+1);
                if p.Results.show_xmaps
                    h_struct(jj).hfig = imagesc(XRF_align_struct.scan(jj).x_pos+p.Results.theta_shift_deg*cosd(XRF_align_struct.scan(jj).theta)*(jj-1));
                    title({['XYMap. int. # ' num2str(XRF_align_struct.scan(jj).scanid)],[ ' th = ' num2str(XRF_align_struct.scan(jj).theta)],...
                        ['xshift = ' num2str(XRF_align_struct.scan(jj).xshifts) 'pixels' ], ['yshift = ' num2str(XRF_align_struct.scan(jj).yshifts) 'pixels']});
                else
                    h_struct(jj).hfig = imagesc('XData',XRF_align_struct.scan(jj).x_pos,'YData',XRF_align_struct.scan(jj).y_pos,'CData',XRF_align_struct.scan(jj).(field_to_show));
                    title({['Map. int. # ' num2str(XRF_align_struct.scan(jj).scanid)],[ ' th = ' num2str(XRF_align_struct.scan(jj).theta)],...
                        ['xshift = ' num2str(XRF_align_struct.scan(jj).xshifts) 'pixels' ], ['yshift = ' num2str(XRF_align_struct.scan(jj).yshifts) 'pixels']});
                end
                colormap jet;
                %axis image;
                colorbar;
                %caxis([-7 4]);
            end
            
            
        end
        
        function hfig = display2Dmap(map2D,varargin)
            p = inputParser;
            
            addRequired(p,'map2D');
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'Xval',[1:size(map2D,2)],@isnumeric);
            addParameter(p,'Yval',[1:size(map2D,1)],@isnumeric);
            addParameter(p,'window',[1 size(map2D,1) 1 size(map2D,2)],@isnumeric);
            addParameter(p,'figTitle','',@ischar);
            addParameter(p,'size_figure',[6 437 1203 314],@isnumeric);
            addParameter(p,'font',10,@isnumeric);
            %addParameter(p,'figAxis','image',@ischar);
            %addParameter(p,'figColormap','hot',@ischar);
            
            parse(p,map2D,varargin{:});
            
            Yval = p.Results.Yval(p.Results.window(1):p.Results.window(2));
            Xval = p.Results.Xval(p.Results.window(3):p.Results.window(4));
            
            figure(p.Results.figNum);
            clf;
            hfig=pcolor(Xval,Yval,map2D(p.Results.window(1):p.Results.window(2),p.Results.window(3):p.Results.window(4)));
            %hfig=imagesc('YData',p.Results.Yval(p.Results.window(1):p.Results.window(2)),'XData',p.Results.Xval(p.Results.window(3):p.Results.window(4)),'CData',map2D(p.Results.window(1):p.Results.window(2),p.Results.window(3):p.Results.window(4)));
            %hfig=imagesc('YData',p.Results.Yval,'XData',p.Results.Xval,'CData',map2D);
           % axis image;
           shading interp;%
           colormap jet;
            colorbar
            title(p.Results.figTitle);
            set(gcf, 'Position', p.Results.size_figure);
            set(gca, 'FontSize', p.Results.font);
        end
        
        function hfig = display2Dmap_toclick(map2D,varargin)
            p = inputParser;
            
            addRequired(p,'map2D');
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            addParameter(p,'figTitle','',@ischar);
            %addParameter(p,'figAxis','image',@ischar);
            %addParameter(p,'figColormap','hot',@ischar);
            
            parse(p,map2D,varargin{:});
            
            figure(p.Results.figNum);
            hfig=imagesc(p.Results.Xval,p.Results.Yval,map2D);
            axis image;
            colormap hot;
            
            pass2click.xaxis = [min(p.Results.Xval) max(p.Results.Xval) size(p.Results.Xval,2)];
            pass2click.yaxis = [min(p.Results.Yval) max(p.Results.Yval) size(p.Results.Yval,1)];
            datacursormode on;
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
            set(dcm_obj, 'UpdateFcn', @click4rock_film);
            
            
        end
        
    end
    
    
end

