classdef ND_read_data
    % This library contains all the functions to read the data from HXN and
    % sector 26
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function [scandata,merlimgs,scandata_pad,hfig1,hfig2] = loadscan_HXN(datapath,scanid,detchan,varargin)
            
            % Create instance of inputParser class.
            p = inputParser;
           
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid', @isnumeric);
            addRequired(p,'detchan', @ischar);
            addParameter(p,'MonChan','sclr1_ch3' ,@ischar);
            addParameter(p,'prefix',{'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'},@iscell);
            addParameter(p,'showmerlin', 1, @isnumeric);
            addParameter(p,'inneraxis', 'z', @ischar);
            addParameter(p,'flyscan', 1, @isnumeric);
            addParameter(p,'domedian', 1, @isnumeric);
            addParameter(p,'ROIinteg', [], @isnumeric);
            addParameter(p,'hotpixels', [], @isnumeric);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);
            addParameter(p,'do_centroids', 1, @isnumeric);
            addParameter(p,'do_padding', 0, @isnumeric);
            addParameter(p,'outerpts_zeropad',0,@isnumeric);
            addParameter(p,'innerpts_zeropad',0,@isnumeric);
            parse(p,datapath,scanid,detchan,varargin{:});
            
            %read the inputs:
            flag_struct.showmerlin   = p.Results.showmerlin;
            flag_struct.flyscan      = p.Results.flyscan;
            flag_struct.domedian     = p.Results.domedian;
            flag_struct.do_centroids = p.Results.do_centroids;
            
            inneraxis = p.Results.inneraxis;
            ROIinteg  = p.Results.ROIinteg;
            hotpixels = p.Results.hotpixels;
            detchan  = p.Results.detchan;
            
            pResults = p.Results;
            
            if p.Results.inneraxis=='y'
                innerchan='zpssy'; %2
                outerchan='zpssx';%1;
            elseif p.Results.inneraxis=='z'
                innerchan='zpssz';%3;
                outerchan='zpssy';%2;
            else
                innerchan='zpssx';%1;
                outerchan='zpssy';%2;
            end
            
            pResults.innerchan = innerchan;
            pResults.outerchan = outerchan;
            
            %innerchan = 2; %zpsy
            %innerchan = 1; %zpsx
            %outerchan=2;
            
            merlimgs = h5read([datapath '/scan_' num2str(scanid) '.h5'],'/entry/instrument/detector/data');
            merlimgs = permute(merlimgs,[2 1 3]);
            numimgs = size(merlimgs,3);
            
            %pixx is horizontal in hutch
            %pixy is vertical in hutch
            pixx = size(merlimgs,2);
            pixy = size(merlimgs,1);
                      
                     
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
                
                if flag_struct.domedian
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
                %diff_data(ii,1) = sum(sum(double(merlimgs(40:149,75:142,ii)).*hotmask(40:149,75:142)));
                
                
            end
            
             if(flag_struct.do_centroids)
                    diff_data = ND_analysis.computeCentroids(merlimgs,ROIinteg);
                    [scandata,scandata_pad] = ND_read_data.getLinearData(pResults,diff_data);
              end
            
           
            %%{
            if flag_struct.showmerlin==1
                if p.Results.do_padding == 1
                    hfig1 =  ND_display_data.showmerlin_function(scandata_pad,pResults,pixx,pixy,2000);
                else
                    hfig1 =  ND_display_data.showmerlin_function(scandata,pResults,pixx,pixy,2000);
                end
            end
            
            if flag_struct.do_centroids==1
                 if p.Results.do_padding == 1
                    hfig2 =  ND_display_data.show_centroid(scandata_pad,pResults,pixx,pixy,2000);
                 else
                    hfig2 =  ND_display_data.show_centroid(scandata,pResults,pixx,pixy,2000);
                 end
            end
            
            %}
            
            
        end
        
        function [scandata,scandata_pad] = getLinearData(pResults,diff_data)

            prefix = pResults.prefix;
           
            if(~isempty(pResults.ROIinteg))
                scandata = zeros(pResults.outerpts,pResults.innerpts,8);
            else
                scandata = zeros(pResults.outerpts, pResults.innerpts, 9);
            end
            
            % get linear data
            filename = [pResults.datapath '/scan_' num2str(pResults.scanid) '.txt'];
            
            [scan,variable_names_cell] = ND_read_data.importfile(filename, 1, pResults.innerpts*pResults.outerpts,'prefix',pResults.prefix);
            temp1 = scan{:,pResults.detchan};
            temp11 = scan{:,pResults.MonChan};
            temp2 = scan{:,pResults.outerchan};
            temp3 = scan{:,pResults.innerchan};
            
           
            
            if(pResults.flyscan)
                for ii=1:pResults.outerpts
                    for jj=1:pResults.innerpts                       
                        scandata(ii,jj,1) = temp1((ii-1)*pResults.innerpts+jj);%temp1(pResults.detchan);
                        scandata(ii,jj,4) = temp11((ii-1)*pResults.innerpts+jj);%temp1(53); %this is correct - IC3
                        scandata(ii,jj,2) = temp2((ii-1)*pResults.innerpts+jj);
                        scandata(ii,jj,3) = temp3((ii-1)*pResults.innerpts+jj);
                        scandata(ii,jj,5) = diff_data((ii-1)*pResults.innerpts+jj,1);
                        scandata(ii,jj,6) = diff_data((ii-1)*pResults.innerpts+jj,2);
                        scandata(ii,jj,7) = diff_data((ii-1)*pResults.innerpts+jj,3);
                        scandata(ii,jj,8) = (ii-1)*pResults.innerpts+jj;
                        if(~isempty(pResults.ROIinteg))
                            scandata(ii,jj,9) = diff_data((ii-1)*pResults.innerpts+jj,4);
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
                        scandata(ii,jj,6) = diff_data((ii-1)*pResults.innerpts+jj,2);
                        scandata(ii,jj,7) = diff_data((ii-1)*pResults.innerpts+jj,3);
                        scandata(ii,jj,8) = (ii-1)*pResults.innerpts+jj;
                    end
                end
            end
            
            if pResults.do_padding == 1
                scandata_pad = zeros(pResults.outerpts+pResults.outerpts_zeropad,pResults.innerpts+pResults.innerpts_zeropad,size(scandata,3));
                
                center_pad = round(size(scandata_pad(:,:,1))./2);
                center_orig = round(size(scandata(:,:,1))./2);
                
                if mod(size(scandata(:,:,1),1),2) == 0
                scandata_pad(center_pad(1)-center_orig(1):center_pad(1)+center_orig(1)-1,...
                    center_pad(2)-center_orig(2):center_pad(2)+center_orig(2)-1,:)= scandata;
                else
                   scandata_pad(center_pad(1)-center_orig(1)+1:center_pad(1)+center_orig(1)-1,...
                    center_pad(2)-center_orig(2)+1:center_pad(2)+center_orig(2)-1,:)= scandata; 
                end
            else
                scandata_pad = scandata;
            end
        end
        
        function [ fly2Dmaps,imapx,imapy,sumim ] = ThetaScan_film(datapath, scanid, detchan,varargin)
            
              p = inputParser; 
            
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid');
            addRequired(p,'detchan', @ischar);
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
            addParameter(p,'do_centroids', 1, @isnumeric);
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

           
            if p.Results.do_padding
                [dataout_orig,imgsout,dataout] = ND_read_data.loadscan_HXN(datapath,scanid(1),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',1);
                [datatrash_orig,imgstrash,datatrash] = ND_read_data.loadscan_HXN(datapath,scanid(1), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',1); %% reads out the photo current;
            else
                [dataout,imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(1),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'showmerlin',p.Results.showmerlin,'do_padding',0);
                [datatrash,imgstrash] = ND_read_data.loadscan_HXN(datapath,scanid(1), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'showmerlin',p.Results.showmerlin,'do_padding',0); %% reads out the photo current;
            end
          
            
            XRF0 = dataout(:,:,1);%./dataout(:,:,4);
            PC0 = datatrash(:,:,1);%./datatrash(:,:,4);
            

            %{
            figure(1);imagesc(XRF0);
            axis image tight off;
            colormap hot;
            title(['Theta: ' num2str(thetalist(1))]);
            %}
            
            xshifts = zeros(numel(scanid),1);
            yshifts = zeros(numel(scanid),1);
            chi2 = zeros(10,10);% even numbers here please
            
            fly2Dmaps.scan(1).XRF = XRF0;
            fly2Dmaps.scan(1).PC = PC0;
            fly2Dmaps.scan(1).xshift = xshifts(1);
            fly2Dmaps.scan(1).yshift = yshifts(1);
            fly2Dmaps.scan(1).theta = thetalist(1);
            fly2Dmaps.xaxis = dataout(:,:,3);
            fly2Dmaps.yaxis = dataout(:,:,2);
            
            if p.Results.do_padding == 0
                for kk = 1:size(dataout,1)
                    for ll = 1:size(dataout,2)
                        fly2Dmaps.ii(kk).jj(ll).im = double(imgsout(:,:,dataout(kk,ll,8)))./(numel(scanid));
                        fly2Dmaps.ii(kk).jj(ll).intensity(1) = dataout(kk,ll,5);
                        fly2Dmaps.ii(kk).jj(ll).SumInt = dataout(kk,ll,5);
                    end
                end
            else
                 for kk = 1:size(dataout,1)
                    for ll = 1:size(dataout,2)
                        
                        if dataout(kk,ll,8) == 0
                            fly2Dmaps.ii(kk).jj(ll).im = 0.0;
                        else
                            fly2Dmaps.ii(kk).jj(ll).im = double(imgsout(:,:,dataout(kk,ll,8)))./(numel(scanid));
                        end
                        fly2Dmaps.ii(kk).jj(ll).intensity(1) = dataout(kk,ll,5);
                        fly2Dmaps.ii(kk).jj(ll).SumInt = dataout(kk,ll,5);
                    end
                end
                
            end
            
            for ii=2:numel(scanid)
                %for ii=2:2
                
                if p.Results.do_padding
                    [dataout_orig,imgsout,dataout] = ND_read_data.loadscan_HXN(datapath,scanid(ii),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',1);
                    [datatrash_orig,imgstrash,datatrash] = ND_read_data.loadscan_HXN(datapath,scanid(ii), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'innerpts_zeropad',p.Results.innerpts_zeropad,'outerpts_zeropad',p.Results.outerpts_zeropad,'showmerlin',p.Results.showmerlin,'do_padding',1); %% reads out the photo current;
                else
                    [dataout,imgsoutt] = ND_read_data.loadscan_HXN(datapath,scanid(ii),detchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'showmerlin',p.Results.showmerlin,'do_padding',0);
                    [datatrash,imgstrash] = ND_read_data.loadscan_HXN(datapath,scanid(ii), XBICchan,'prefix',prefix,'innerpts',p.Results.innerpts,'outerpts',p.Results.outerpts,'showmerlin',p.Results.showmerlin,'do_padding',0); %% reads out the photo current;
                end
                
               
                XRF1= dataout(:,:,1);%./dataout(:,:,4); %normalize by monitor
                PC1 = datatrash(:,:,1);%./datatrash(:,:,4); %photo current normalized by monitor;
               
                for jj = 1:size(chi2,1)
                    for kk = 1:size(chi2,2)
                        chi2(jj,kk) = sum(sum((XRF0 - circshift(XRF1,[jj-size(chi2,1)/2,kk-size(chi2,2)/2])).^2,'omitnan'),'omitnan');
                       
                    end
                end
                [xcen,ycen,intemp] = find(chi2==min(min(chi2)));
                %figure(4);imagesc(chi2);axis image tight;colormap hot
                xshifts(ii) = xshifts(ii-1)+xcen-size(chi2,1)/2;
                yshifts(ii) = yshifts(ii-1)+ycen-size(chi2,2)/2;
           
        
                
                fly2Dmaps.scan(ii).XRF = circshift(XRF1,[xshifts(ii),yshifts(ii)]);
                fly2Dmaps.scan(ii).PC = circshift(PC1, [xshifts(ii),yshifts(ii)]);
                fly2Dmaps.scan(ii).xshift = xshifts(ii);
                fly2Dmaps.scan(ii).yshift = yshifts(ii);
                fly2Dmaps.scan(ii).theta = thetalist(ii);
                ccdnums = circshift(dataout(:,:,8),[xshifts(ii),yshifts(ii)]);
                tempints = circshift(dataout(:,:,5),[xshifts(ii),yshifts(ii)]);
                
                if p.Results.do_padding == 0
                    for kk = 1:size(dataout,1)
                        for ll = 1:size(dataout,2)
                            fly2Dmaps.ii(kk).jj(ll).im = fly2Dmaps.ii(kk).jj(ll).im + double(imgsout(:,:,ccdnums(kk,ll)))./(numel(scanid));
                            fly2Dmaps.ii(kk).jj(ll).intensity(ii) = tempints(kk,ll);
                            fly2Dmaps.ii(kk).jj(ll).SumInt = fly2Dmaps.ii(kk).jj(ll).SumInt+tempints(kk,ll);
                        end
                    end
                else
                    for kk = 1:size(dataout,1)
                        for ll = 1:size(dataout,2)
                            
                            if ccdnums(kk,ll) == 0
                                fly2Dmaps.ii(kk).jj(ll).im = 0.0;
                            else
                                fly2Dmaps.ii(kk).jj(ll).im = fly2Dmaps.ii(kk).jj(ll).im + double(imgsout(:,:,ccdnums(kk,ll)))./(numel(scanid));
                            end
                            fly2Dmaps.ii(kk).jj(ll).intensity(ii) = tempints(kk,ll);
                            fly2Dmaps.ii(kk).jj(ll).SumInt = fly2Dmaps.ii(kk).jj(ll).SumInt+tempints(kk,ll);
                        end
                    end
                    
                end
                
                
               
                XRF0=XRF1;
                %pause(1);
                
                if not(p.Results.plotflag)
                    ND_display_data.display2Dmap(XRF1,'figNum',3,'figTitle',['Unshifted Theta: ' num2str(thetalist(ii))]);
                    ND_display_data.display2Dmap(circshift(XRF1,[xshifts(ii),yshifts(ii)]),'figNum',2,'figTitle',['Theta: ' num2str(thetalist(ii)) ' x shift: ' num2str(yshifts(ii)) ' y shift: ' num2str(xshifts(ii))]);
                end
                
            end
            %{
            %      Creating package to pass to click...
            %     pass2click.fly2dscanlist = fly2Dscanlist;
            %     pass2click.txtfilepath = '/GPFS/XF03ID1/users/2017Q2/Hruszkewycz_2017Q2/Data';
            %    clear pass2click;
            pass2click.fly2Dmaps = fly2Dmaps;
            %     pass2click.imrows = outerpts;
            %     pass2click.imcols = innerpts;
            %}
            
            tempim = zeros(size(dataout,1),size(dataout,2));
            tempxcen = zeros(size(dataout,1),size(dataout,2));
            tempycen = zeros(size(dataout,1),size(dataout,2));
            
            for kk = 1:size(dataout,1)
                for ll = 1:size(dataout,2)
                    tempim(kk,ll) = fly2Dmaps.ii(kk).jj(ll).SumInt;
                    imgin =fly2Dmaps.ii(kk).jj(ll).im;
                    
                    line1=sum(imgin,1);  % vertical
                    line2=sum(imgin,2);  % horizontal
                    sumt = sum(sum(imgin));
                    for mm=1:size(line1,2)
                        tempycen(kk,ll)=tempycen(kk,ll)+mm*line1(mm)/sumt;
                    end
                    for mm=1:size(line2,1)
                        tempxcen(kk,ll)=tempxcen(kk,ll)+mm*line2(mm)/sumt;
                    end
                    
                end
            end
            
            fly2Dmaps.Xcentroids = tempxcen;
            fly2Dmaps.Ycentroids = tempycen;
            fly2Dmaps.imapx = dataout(1,:,3);
            fly2Dmaps.imapy = dataout(:,1,2);
            imapx = dataout(1,:,3);
            imapy = dataout(:,1,2);
            sumim = tempim;
            
            if not(p.Results.plotflag)
                ND_display_data.display2Dmap(XRF0,'figNum',1,'figTitle',['Theta: ' num2str(thetalist(1))]);
                hfig =  ND_display_data.display2Dmap_toclick(tempim,'figNum',20,'Xval',imapx,'Yval',imapy,'figTitle','Sum Int.');
                ND_display_data.display2Dmap(tempxcen,'figNum',22,'Xval',imapx,'Yval',imapy,'figTitle','Row Centroids');
                ND_display_data.display2Dmap(tempycen,'figNum',23,'Xval',imapx,'Yval',imapy,'figTitle',['Column Centroids ']);
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
