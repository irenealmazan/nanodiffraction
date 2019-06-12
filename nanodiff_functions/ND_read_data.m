classdef ND_read_data
    % This library contains all the functions to read the data from HXN and
    % sector 26
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function [scandata,ccdimgs] = loadscan_HXN(datapath,scanid,detchan,varargin)
            % loadscan_HXN reads the ccd stored in .h5 for a single scan
            % (2D map in real space at 1 angle) and creates 
            
            % Create instance of inputParser class.
            p = inputParser;
           
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid', @isnumeric);
            addRequired(p,'detchan', @ischar);
            addParameter(p,'use_fitXRF',0,@isnumeric);
            addParameter(p,'MonChan','sclr1_ch3' ,@ischar);
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            addParameter(p,'showmerlin', 0, @isnumeric);
            addParameter(p,'inneraxis', 'z', @ischar);
            addParameter(p,'flyscan', 1, @isnumeric);
            addParameter(p,'domedian', 1, @isnumeric);
            addParameter(p,'ROIinteg', [], @isnumeric);
            addParameter(p,'hotpixels', [], @isnumeric);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);          
            addParameter(p,'do_padding', 0, @isnumeric);
            addParameter(p,'outerpts_zeropad',0,@isnumeric);
            addParameter(p,'innerpts_zeropad',0,@isnumeric);
            parse(p,datapath,scanid,detchan,varargin{:});
            
            %read the inputs:
            flag_struct.showmerlin   = p.Results.showmerlin;
            flag_struct.flyscan      = p.Results.flyscan;
            flag_struct.domedian     = p.Results.domedian;
         
            
            pResults = p.Results;
            
            if p.Results.inneraxis=='y'
                innerchan='zpssy';%2
                outerchan='zpssx';%1;
            elseif p.Results.inneraxis=='z'
                innerchan='zpssz';%3;
                outerchan='zpssy';%2;
            else
                innerchan='zpssx';%1;
                outerchan='zpssy';%2;
            end
            
            [ccdimgs,data_sum_ccd] = ND_read_data.readCCDimages(p.Results.datapath,p.Results.scanid,p.Results.hotpixels,p.Results.domedian);
                 
            numimgs = size(ccdimgs,3);
            
             pResults.innerchan = innerchan;
            pResults.outerchan = outerchan;
                        
             if not(p.Results.innerpts)==0
                innerpts  = p.Results.innerpts;                
                if not(p.Results.outerpts)==0
                    outerpts=round(numimgs/innerpts);
                else
                    outerpts = p.Results.outerpts;
                end
            else
                disp('did not specify inner points')
                return
             end
            
              %outerpts=round(numimgs/innerpts);
            pResults.outerpts = outerpts;
            
            [scandata] = ND_read_data.getLinearDatain2DMap(pResults,data_sum_ccd);
            
            if p.Results.do_padding
                scandata_pad =  ND_data_processing.padData(scandata,p.Results.outerpts_zeropad,p.Results.innerpts_zeropad);
                scandata = scandata_pad;
            end

            
        end
        
        function [merlimgs,diff_data]= readCCDimages(datapath,scanid,hotpixels,do_median)
            
            
            merlimgs = h5read([datapath '/scan_' num2str(scanid) '.h5'],'/entry/instrument/detector/data');
            merlimgs = permute(merlimgs,[2 1 3]);
            numimgs = size(merlimgs,3);
            
            %pixx is horizontal in hutch
            %pixy is vertical in hutch
            pixx = size(merlimgs,2);
            pixy = size(merlimgs,1);
                      
                     
           
            % store ccd images
            
            for ii = 1:numimgs
                %    ccd = flipud(rot90(double(merlimgs(:,:,ii)))); %old version of hd5
                %    reader
                ccd = double(merlimgs(:,:,ii));
                
                if(~isempty(hotpixels))
                    for ll = size(hotpixels,1)
                        ccd(hotpixels(ll,2),hotpixels(ll,1)) = 0;
                    end
                end
                
                if do_median
                    ccd1 = zeros(size(ccd,1),size(ccd,2),5);
                    ccd1(:,:,1) = ccd;
                    ccd1(:,:,2) = circshift(ccd,[0,1]);
                    ccd1(:,:,3) = circshift(ccd,[1,0]);
                    ccd1(:,:,4) = circshift(ccd,[0,-1]);
                    ccd1(:,:,5) = circshift(ccd,[-1,0]);
                    ccd2 = median(ccd1,3);
                    ccdmask = ccd>ccd2+50;
                    ccd = ccd.*(1-ccdmask)+ccd2.*ccdmask;
                end
                
                merlimgs(:,:,ii) = ccd;
                diff_data(ii,1) = sum(sum(ccd));
                
                
            end
            
        end
        
        function [scandata] = getLinearDatain2DMap(pResults,diff_data)

           
            if(~isempty(pResults.ROIinteg))
                scandata = zeros(pResults.outerpts,pResults.innerpts,6);
            else
                scandata = zeros(pResults.outerpts, pResults.innerpts, 7);
            end
            
            % get linear data
            filename = [pResults.datapath '/scan_' num2str(pResults.scanid) '.txt'];
            
            [scan,variable_names_cell] = ND_read_data.importfile(filename, 1, pResults.innerpts*pResults.outerpts,'prefix',pResults.prefix);

            if pResults.use_fitXRF == 0
                temp1 = scan{:,pResults.detchan}; % fluorescence map
            else
                % read fluorescence maps from tiff files:
               t = Tiff([pResults.datapath '/xrf_fit/output_tiff_scan2D_' num2str(pResults.scanid) '/' pResults.detchan num2str(pResults.scanid) '_norm.tiff'],'r');
               temp = read(t);
               temp1 = double(temp);
            end
            
            temp11 = scan{:,pResults.MonChan}; % monitor map (I don't know what chanel it is)
            temp2 = scan{:,pResults.outerchan}; % motor piezo values for the inner loop of the raster scan
            temp3 = scan{:,pResults.innerchan}; % motor piezo values for the outer loop of the raster scan. I don't know the convention at HXN row/columns - inner/outer
            
           
            % construct the 2D maps in real space 
            if(pResults.flyscan)
                for ii=1:pResults.outerpts
                    for jj=1:pResults.innerpts
                        if pResults.use_fitXRF == 0
                            scandata(ii,jj,1) = temp1((ii-1)*pResults.innerpts+jj);% Fluorescence map
                        else
                            scandata(ii,jj,1) = temp1(ii,jj);% Fluorescence map
                        end
                        scandata(ii,jj,4) = temp11((ii-1)*pResults.innerpts+jj);% monitor map??
                        scandata(ii,jj,2) = temp2((ii-1)*pResults.innerpts+jj); % map of the motor piezo positions
                        scandata(ii,jj,3) = temp3((ii-1)*pResults.innerpts+jj);% map of the motor piezo positions
                        scandata(ii,jj,5) = diff_data((ii-1)*pResults.innerpts+jj,1); % diffraction map for each angle = sum(sum(ccd))
                        %scandata(ii,jj,6) = diff_data((ii-1)*pResults.innerpts+jj,2); % vertical centroid
                        %scandata(ii,jj,7) = diff_data((ii-1)*pResults.innerpts+jj,3); % horizontal centroid
                        scandata(ii,jj,6) = (ii-1)*pResults.innerpts+jj;
                        if(~isempty(pResults.ROIinteg))
                            scandata(ii,jj,7) = diff_data((ii-1)*pResults.innerpts+jj,4);
                        end
                    end
                end
            else
                for ii=1:pResults.outerpts
                    for jj=1:pResults.innerpts
                        scandata(ii,jj,2) = ii*0.1;
                        scandata(ii,jj,3) = jj*0.1;
                        scandata(ii,jj,4) = 1;
                        scandata(ii,jj,5) = diff_data((ii-1)*pResults.innerpts+jj,1);
                        %scandata(ii,jj,6) = diff_data((ii-1)*pResults.innerpts+jj,2);
                        %scandata(ii,jj,7) = diff_data((ii-1)*pResults.innerpts+jj,3);
                        scandata(ii,jj,6) = (ii-1)*pResults.innerpts+jj;
                    end
                end
            end
            
         
        end
        
        function [fly2Dmap_return] = read_tiff_and_getLinearDatain2DMap(datapath,detchan,scanid,fly2Dmap,varargin)
            
            p = inputParser;
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'detchan', @ischar);
            addRequired(p,'scanid', @isnumeric);
            addRequired(p,'fly2Dmap', @isstruct);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);
            addParameter(p,'do_padding', 0, @isnumeric);
            addParameter(p,'outerpts_zeropad',0,@isnumeric);
            addParameter(p,'innerpts_zeropad',0,@isnumeric);
            
            parse(p,datapath,detchan,scanid,fly2Dmap,varargin{:});
            
            fly2Dmap_return = p.Results.fly2Dmap;
            
            for ii = 1:numel(p.Results.scanid)
                t = Tiff([p.Results.datapath '/xrf_fit/output_tiff_scan2D_' num2str(p.Results.scanid(ii)) '/' p.Results.detchan '_' num2str(p.Results.scanid(ii)) '.tiff'],'r');
                temp = read(t);
                
                if p.Results.do_padding
                    if isempty(strfind(p.Results.detchan,'pos'))
                        scandata_pad = ND_data_processing.padData(double(temp),p.Results.outerpts_zeropad,p.Results.innerpts_zeropad);
                    else
                        if ~isempty(strfind(p.Results.detchan,'x_pos'))
                            scandata_pad = ND_data_processing.padPosData(double(temp),p.Results.outerpts+p.Results.outerpts_zeropad,p.Results.innerpts+p.Results.innerpts_zeropad,1) ;
                        elseif ~isempty(strfind(p.Results.detchan,'y_pos'))
                            scandata_pad = ND_data_processing.padPosData(double(temp),p.Results.innerpts+p.Results.innerpts_zeropad,p.Results.outerpts+p.Results.outerpts_zeropad,2) ;
                        end
                    end
                    
                    fly2Dmap_return.scan(ii).(p.Results.detchan) =  scandata_pad;
                else
                    fly2Dmap_return.scan(ii).(p.Results.detchan) = double(temp);
                end
                
                
            end
            
            
            
           
            
            
        end
        
        
        
        function [ fly2Dmaps] = ThetaScan_film(datapath, scanid, detchan,varargin)
            
              p = inputParser; 
            
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid');
            addRequired(p,'detchan', @ischar);
            addParameter(p,'use_fitXRF',0,@isnumeric);
            addParameter(p,'thetalist',[]);
            addParameter(p,'MonChan','sclr1_ch3',@ischar);
            addParameter(p,'XBICchan','sclr1_ch3',@ischar);
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            addParameter(p,'showmerlin', 1, @isnumeric);
            addParameter(p,'plotflag', 0, @isnumeric);
            addParameter(p,'inneraxis', 'z', @ischar);
            addParameter(p,'flyscan', 1, @isnumeric);
            addParameter(p,'domedian', 1, @isnumeric);
            addParameter(p,'ROIinteg', [], @isnumeric);
            addParameter(p,'hotpixels', [], @isnumeric);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);
            addParameter(p,'do_align', 0, @isnumeric);
            addParameter(p,'do_Ref_XRF0', 0, @isnumeric);
            addParameter(p,'do_padding', 0, @isnumeric);
            addParameter(p,'outerpts_zeropad',0,@isnumeric);
            addParameter(p,'innerpts_zeropad',0,@isnumeric);
            parse(p,datapath,scanid,detchan,varargin{:});
             
            %read the inputs:
%             flag_struct.showmerlin   = p.Results.showmerlin;
%             flag_struct.flyscan      = p.Results.flyscan;
%             flag_struct.domedian     = p.Results.domedian;
%             flag_struct.do_centroids = p.Results.do_centroids;
            
            thetalist = p.Results.thetalist;
            inneraxis = p.Results.inneraxis;
            ROIinteg  = p.Results.ROIinteg;
            hotpixels = p.Results.hotpixels;
            detchan  = p.Results.detchan;
            MonChan = p.Results.MonChan;
            prefix = p.Results.prefix;
            XBICchan = p.Results.XBICchan;
            
            pResults = p.Results;
            
            
            
            
            if isempty(thetalist)
                thetalist = 76.350+0.05*((1:numel(scanid))-1);
                %thetalist = 76.5+0.1*((1:numel(fly2Dscanlist))-1);
            end

           
            [dataout,imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(1),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',p.Results.use_fitXRF);
            [datatrash] = ND_read_data.loadscan_HXN(datapath,scanid(1), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',0); %% reads out the photo current;
                                   
            XRF0 = dataout(:,:,1);%./dataout(:,:,4);
            PC0 = datatrash(:,:,1);%./datatrash(:,:,4);
            
            XRF_struct.scan(1).XRF = XRF0;
            
          
            
            xshifts = zeros(numel(scanid),1);
            yshifts = zeros(numel(scanid),1);
            %chi2 = zeros(10,10);% even numbers here please
            
            fly2Dmaps.scan(1).XRF = XRF0;
            fly2Dmaps.scan(1).PC = PC0;
            fly2Dmaps.scan(1).xshift = xshifts(1);
            fly2Dmaps.scan(1).yshift = yshifts(1);
            fly2Dmaps.scan(1).theta = thetalist(1);
            fly2Dmaps.scan(1).imgsout = imgsout;
            fly2Dmaps.xaxis = dataout(:,:,3);
            fly2Dmaps.yaxis = dataout(:,:,2);
            
              for kk = 1:size(dataout,1)
                for ll = 1:size(dataout,2)
                    %if p.Results.do_padding == 0
                        if dataout(kk,ll,6) == 0
                            fly2Dmaps.ii(kk).jj(ll).im = 0.0;
                        else
                            fly2Dmaps.ii(kk).jj(ll).im = double(imgsout(:,:,dataout(kk,ll,6)))./(numel(scanid)); % sum_angle ccd_angle -> 512x512
                        end
                    %end
                    fly2Dmaps.ii(kk).jj(ll).intensity(1) = dataout(kk,ll,5);%sum(sum(ccd))_angle ->scalar
                    fly2Dmaps.ii(kk).jj(ll).SumInt = dataout(kk,ll,5); % sum_angle(sum(sum(ccd))) -> scalar
                end
            end
            
            for ii=2:numel(scanid)
       
                
                [dataout,imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(ii),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',p.Results.use_fitXRF);
                [datatrash] = ND_read_data.loadscan_HXN(datapath,scanid(ii), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',0); %% reads out the photo current;
               
                
                XRF1= dataout(:,:,1);%./dataout(:,:,4); %normalize by monitor
                PC1 = datatrash(:,:,1);%./datatrash(:,:,4); %photo current normalized by monitor;
                
                
                
                if p.Results.do_align == 0
                    XRF_struct.scan(2).XRF = XRF0;                 
                else
                    XRF_struct.scan(2).XRF = XRF1;  
                    if p.Results.do_Ref_XRF0 == 0
                        XRF_struct.scan(1).XRF = fly2Dmaps.scan(ii-1).XRF;
                    else
                        XRF_struct.scan(1).XRF = XRF0;
                    end
                end
                
                
                
                XRF_align_struct = ND_analysis.doAlignment(XRF_struct,'scanid',p.Results.scanid,'thetalist',p.Results.thetalist,'do_Ref_XRF0',1);    
     %{           
               
                for jj = 1:size(chi2,1)
                    for kk = 1:size(chi2,2)
                        chi2(jj,kk) = sum(sum((XRF0 - circshift(XRF1,[jj-size(chi2,1)/2,kk-size(chi2,2)/2])).^2,'omitnan'),'omitnan');                       
                    end
                end
                [xcen,ycen,intemp] = find(chi2==min(min(chi2)));
                %figure(4);imagesc(chi2);axis image tight;colormap hot
                xshifts(ii) = xshifts(ii-1)+xcen-size(chi2,1)/2;
                yshifts(ii) = yshifts(ii-1)+ycen-size(chi2,2)/2;
           
     %}   
                
                fly2Dmaps.scan(ii).XRF = XRF_align_struct.scan(2).XRF;%circshift(XRF1,[yshifts(ii),xshifts(ii)]);
                fly2Dmaps.scan(ii).PC =  circshift(PC1, [XRF_align_struct.scan(2).yshifts,XRF_align_struct.scan(2).xshifts]);
                fly2Dmaps.scan(ii).indexMap = XRF_align_struct.scan(2).indexMap;%xshifts(ii);
                fly2Dmaps.scan(ii).xshift = XRF_align_struct.scan(2).xshifts;%xshifts(ii);
                fly2Dmaps.scan(ii).yshift = XRF_align_struct.scan(2).yshifts;%yshifts(ii);
                fly2Dmaps.scan(ii).theta = thetalist(ii);
                fly2Dmaps.scan(ii).imgsout = imgsout;
               
                
                ccdnums = circshift(dataout(:,:,6),[XRF_align_struct.scan(2).yshifts,XRF_align_struct.scan(2).xshifts]);
                tempints = circshift(dataout(:,:,5),[XRF_align_struct.scan(2).yshifts,XRF_align_struct.scan(2).xshifts]);
                
                for kk = 1:size(dataout,1)
                    for ll = 1:size(dataout,2)
                        %if p.Results.do_padding == 0
                            if ccdnums(kk,ll) == 0
                                fly2Dmaps.ii(kk).jj(ll).im = 0.0;
                            else
                                fly2Dmaps.ii(kk).jj(ll).im = fly2Dmaps.ii(kk).jj(ll).im + double(imgsout(:,:,ccdnums(kk,ll)))./(numel(scanid));
                            end
                       % end
                        fly2Dmaps.ii(kk).jj(ll).intensity(ii) = tempints(kk,ll);
                        fly2Dmaps.ii(kk).jj(ll).SumInt = fly2Dmaps.ii(kk).jj(ll).SumInt+tempints(kk,ll);
                    end
                end
                
            
                if not(p.Results.plotflag)
                    ND_display_data.display2Dmap(XRF0,'figNum',2,'figTitle',['XRF0: ' num2str(thetalist(ii))]);
                    ND_display_data.display2Dmap(circshift(XRF1,[fly2Dmaps.scan(ii).yshift,fly2Dmaps.scan(ii).xshift]),'figNum',3,'figTitle',['Theta: ' num2str(thetalist(ii)) ' row shift: ' num2str(fly2Dmaps.scan(ii).yshift) ' col shift: ' num2str(fly2Dmaps.scan(ii).xshift)]);
                end
                
            end
           
         
            fly2Dmaps.imapx = dataout(1,:,3);
            fly2Dmaps.imapy = dataout(:,1,2);
            %imapx = dataout(1,:,3);
            %imapy = dataout(:,1,2);
            %sumim = tempim;
            
           
        end
        
        function [ fly2Dmaps] = ThetaScan_film_onlyread(datapath, scanid, detchan,varargin)
            
            p = inputParser;
            
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid');
            addRequired(p,'detchan', @ischar);
            addParameter(p,'use_fitXRF',0,@isnumeric);
            addParameter(p,'thetalist',[]);
            addParameter(p,'MonChan','sclr1_ch3',@ischar);
            addParameter(p,'XBICchan','sclr1_ch3',@ischar);
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            addParameter(p,'showmerlin', 1, @isnumeric);
            addParameter(p,'plotflag', 0, @isnumeric);
            addParameter(p,'inneraxis', 'z', @ischar);
            addParameter(p,'flyscan', 1, @isnumeric);
            addParameter(p,'domedian', 1, @isnumeric);
            addParameter(p,'ROIinteg', [], @isnumeric);
            addParameter(p,'hotpixels', [], @isnumeric);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);
            addParameter(p,'do_align', 0, @isnumeric);
            addParameter(p,'do_Ref_XRF0', 0, @isnumeric);
            addParameter(p,'do_padding', 0, @isnumeric);
            addParameter(p,'outerpts_zeropad',0,@isnumeric);
            addParameter(p,'innerpts_zeropad',0,@isnumeric);
            parse(p,datapath,scanid,detchan,varargin{:});
            
            %read the inputs:
            %             flag_struct.showmerlin   = p.Results.showmerlin;
            %             flag_struct.flyscan      = p.Results.flyscan;
            %             flag_struct.domedian     = p.Results.domedian;
            %             flag_struct.do_centroids = p.Results.do_centroids;
            
            thetalist = p.Results.thetalist;
            inneraxis = p.Results.inneraxis;
            ROIinteg  = p.Results.ROIinteg;
            hotpixels = p.Results.hotpixels;
            detchan  = p.Results.detchan;
            MonChan = p.Results.MonChan;
            prefix = p.Results.prefix;
            XBICchan = p.Results.XBICchan;
            
            pResults = p.Results;
            
            
            
            
            if isempty(thetalist)
                thetalist = 76.350+0.05*((1:numel(scanid))-1);
                %thetalist = 76.5+0.1*((1:numel(fly2Dscanlist))-1);
            end
            
            
            for ii=1:numel(scanid)
                
                
                [dataout,imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(ii),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',p.Results.use_fitXRF);
                [datatrash] = ND_read_data.loadscan_HXN(datapath,scanid(ii), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',p.Results.do_padding,'use_fitXRF',0); %% reads out the photo current;
                
                
                XRF1= dataout(:,:,1);%./dataout(:,:,4); %normalize by monitor
                PC1 = datatrash(:,:,1);%./datatrash(:,:,4); %photo current normalized by monitor;
                
                
                fly2Dmaps.scan(ii).XRF = dataout(:,:,1);%./dataout(:,:,4); %normalize by monitor
                fly2Dmaps.scan(ii).PC  = datatrash(:,:,1);%./datatrash(:,:,4); %photo current normalized by monitor;
                
                fly2Dmaps.scan(ii).theta = thetalist(ii);
                fly2Dmaps.scan(ii).imgsout = imgsout;
                fly2Dmaps.scan(ii).dataout = dataout;
                
                
            end
            
            
            fly2Dmaps.imapx = dataout(1,:,3);
            fly2Dmaps.imapy = dataout(:,1,2);
            %imapx = dataout(1,:,3);
            %imapy = dataout(:,1,2);
            %sumim = tempim;
            
            
        end
   
        
        function [ccd_to_plot,sum_ccd] = readCCD_oneScan_onepixel(datapath,scanid,innerpts,pixel)
            
            
            for kk = 1:numel(scanid)
                merlimgs = h5read([datapath '/scan_' num2str(scanid(kk)) '.h5'],'/entry/instrument/detector/data');
                merlimgs = permute(merlimgs,[2 1 3]);
                num_ccd = (pixel(1)-1)*innerpts+pixel(2);
                
                ccd_to_plot(kk).ccd = merlimgs(:,:, num_ccd );
            end
            
            sum_ccd = zeros(size(merlimgs(:,:,num_ccd)));
            
            for kk = 1:numel(scanid)
               sum_ccd = sum_ccd +  double(ccd_to_plot(kk).ccd);
            end
            
           
            
        end
     
       
        function [scan,variable_names_cell] = importfile(filename, startRow, endRow,varargin)
            
            %IMPORTFILE Import numeric data from a text file as a matrix.
            %   SCAN45417 = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
            %   the default selection.
            %
            %   SCAN45417 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
            %   STARTROW through ENDROW of text file FILENAME.
            %
            % Example:
            %   scan45417 = importfile('scan_45417.txt', 2, 862);
            %
            %    See also TEXTSCAN.
            
            % Auto-generated by MATLAB on 2019/05/07 11:22:14
            
            
            p = inputParser;
            
            addRequired(p,'filename');
            addRequired(p,'startRow');
            addRequired(p,'endRow');
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            
            parse(p,filename,startRow,endRow,varargin{:});

            
            filename = p.Results.filename;
            startRow = p.Results.startRow;
            endRow = p.Results.endRow;
            prefix = p.Results.prefix;
            
             fileID = fopen(filename,'r');
            
            % Initialize variables.
            delimiter = '\t';
            if nargin<=2
                startRow = 2;
                endRow = inf;
            end
             
            % Count number of columns
            
           [var_struct, variable_names_cell] =  ND_read_data.prepare_header(fileID,'prefix',prefix);
           
            formatSpec = '';
            for kk = 1:numel(var_struct)
                if strcmp(var_struct(kk).varname,'time') == 1
                    formatSpec = [formatSpec '%s'];
                elseif strcmp(var_struct(kk).varname,'sclr1_calculations_calc5_equation') == 1
                    formatSpec = [formatSpec '%C'];
                else
                formatSpec = [formatSpec '%f'];%%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
                end
            end
             formatSpec = [formatSpec '%[^\n\r]'];
            %% Open the text file.
            
            
            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID,formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            
            % Close the text file.
            fclose(fileID);
          
    
            % Create output variable
            scan = table(dataArray{1:end-1},'VariableNames',variable_names_cell);
        end
        
        function [var_struct, variable_names_cell] =  prepare_header(fileID,varargin)
            
            p = inputParser;
            
            addRequired(p,'fileID');
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            
            parse(p,fileID,varargin{:});
            frewind(p.Results.fileID)
            
            tline = fgetl(p.Results.fileID);
            
            % count number of columns:
            %prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};
            prefix = p.Results.prefix;
            
            count_col = 0;
            count_var = 1;
            for jj = 1:numel(prefix)-1
                var_index = strfind(tline,prefix{jj});
                if ~isempty(var_index)
                    if ismember(prefix{jj},{'alive','dead','elapsed_time','scaler_alive','time'})
                        var_index = var_index(1);
                    end
                end
                count_col = count_col + numel(var_index);%count(tline,prefix{jj});
                
                for kk = 1:numel(var_index)
                    space_array = isspace(tline(var_index(kk):var_index(kk)+15));
                    index_nonzero = find(space_array == 1);
                    
                    if sum(index_nonzero) == 0
                        space_array = isspace(tline(var_index(kk):var_index(kk)+45));
                        index_nonzero = find(space_array == 1);
                    end
                    
                    
                    var_struct(count_var).varname = sscanf(tline(var_index(kk):var_index(kk)+index_nonzero(1)-1),'%s');
                    
                    count_var = count_var + 1;
                    
                end
                
                
                
            end
            
            jj = numel(prefix);
            
            var_index = strfind(tline,prefix{jj});
            if ismember(prefix{jj},{'alive','dead','elapsed_time','scaler_alive','time'})
                var_index = var_index(1);
            end
            
            count_col = count_col + numel(var_index);%count(tline,prefix{jj});
            
            for kk = 1:numel(var_index)-1
                space_array = isspace(tline(var_index(kk):end));
                index_nonzero = find(space_array == 1);
                
                
                var_struct(count_var).varname = sscanf(tline(var_index(kk):var_index(kk)+index_nonzero(1)-1),'%s');
                
                count_var = count_var + 1;
            end
            
            kk = numel(var_index);
            var_struct(count_var).varname = sscanf(tline(var_index(kk):end),'%s');
            
            
            for jj =1:numel(var_struct)
                variable_names_cell{jj} = var_struct(jj).varname;
            end
            
        end
        
    end
end
