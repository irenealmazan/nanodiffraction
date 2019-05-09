
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
       
        
        function [rock_curve,thetalist] = displayRockCurveMaps(scanid,fly2Dmaps,varargin)
            p = inputParser;
            
            addRequired(p,'scanid');
            addRequired(p,'fly2DMaps');
            addParameter(p,'figNum',1,@isnumeric);
            addParameter(p,'do_mask',0,@isnumeric);
            addParameter(p,'titleXray',1,@isnumeric);
            addParameter(p,'titleFluo',1,@isnumeric);
            addParameter(p,'titleXBIC',1,@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,scanid,fly2Dmaps,varargin{:});
            
            numRows = 3;
            numCols = 5;
            
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
                    set(gcf, 'Position', [100 100 1300 800]);
                    set(gcf, 'PaperPosition', [.25 6.75 4.5 4]);
                    set(gcf, 'Color' ,'w');
                    set(gca, 'FontSize',8);
                    
                    counter = counter + 1;
                    
                end
                
                %subplot(numel(scanid),4,jj);
                subplot(numRows,numCols,mod(jj-1,numCols)+1);
                h_struct(jj).hfig = imagesc(diff_map);
                title({['Diffr. int. # ' num2str(scanid(jj))],[ ' th = ' num2str(thetalist(jj))]});
                colormap jet;
                axis image;
                colorbar;
                %{
                if jj == 1
                    colorbar;
                end
                %}
                
                %subplot(numel(scanid),4,2*4+jj);
                subplot(numRows,numCols,(numCols)+mod(jj-1,numCols)+1);
                h_struct(jj).hfig = imagesc(fly2Dmaps.scan(jj).XRF./fly2Dmaps.scan(jj).PC);
                title({['X-ray Fluo.  # ' num2str(scanid(jj))],[ ' th = ' num2str(thetalist(jj)) 'norm by PC']});
                colormap jet;
                axis image;
                colorbar;
                %{
                if jj == 1
                    colorbar;
                end
                %}
                
                subplot(numRows,numCols,(numCols)*2+mod(jj-1,numCols)+1);
                h_struct(jj).hfig = imagesc(PC_map);
                title({['Photo. current # ' num2str(scanid(jj))] , [' th = ' num2str(thetalist(jj))]});
                colormap jet;
                axis image;
                colorbar;
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
        
         function displayRockCurveShift(scanid,fly2Dmaps,varargin)
            p = inputParser;
            
            addRequired(p,'scanid');
            addRequired(p,'fly2DMaps');
            addParameter(p,'figNum',1,@isnumeric);
            %addParameter(p,'Xval',[1:size(map2D,1)],@isnumeric);
            %addParameter(p,'Yval',[1:size(map2D,2)],@isnumeric);
            %addParameter(p,'figTitle','',@ischar);
            
            parse(p,scanid,fly2Dmaps,varargin{:});
            
           
            
            
            for jj = 1:numel(scanid)
                xshiftlist(jj) = fly2Dmaps.scan(jj).xshift;
                yshiftlist(jj) = fly2Dmaps.scan(jj).yshift;
                thetalist(jj) = fly2Dmaps.scan(jj).theta;
               
            end
            
            figure(p.Results.figNum);           
            clf;
            
            subplot(2,1,1);
            %plot(thetalist,xshiftlist);
            %title('Xshift vs theta');
            plot(scanid,xshiftlist);
            title('Xshift vs scanid');
                
            subplot(2,1,2);
            %plot(thetalist,yshiftlist);
            %title('Yshift vs theta');
            plot(scanid,yshiftlist);
            title('Yshift vs scanid');    
            
        end
        
        function hfig = display2Dmap(map2D,varargin)
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
            clf;
            hfig=imagesc(p.Results.Xval,p.Results.Yval,map2D);
            axis image;
            colormap jet;
            colorbar
            title(p.Results.figTitle);

            
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

