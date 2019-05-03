classdef ND_read_data
    % This library contains all the functions to read the data from HXN and
    % sector 26
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function [scandata,merlimgs,hfig1,hfig2] = loadscan_HXN(datapath,scanid,detchan,varargin)
            
            % Create instance of inputParser class.
            p = inputParser; 
            
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid', @isnumeric);
            addRequired(p,'detchan', @isnumeric);
            addParameter(p,'showmerlin', 1, @isnumeric);
            addParameter(p,'inneraxis', 'z', @ischar);
            addParameter(p,'flyscan', 1, @isnumeric);
            addParameter(p,'domedian', 1, @isnumeric);
            addParameter(p,'ROIinteg', [], @isnumeric);
            addParameter(p,'hotpixels', [], @isnumeric);
            addParameter(p,'innerpts', 0, @isnumeric);
            addParameter(p,'outerpts', 0, @isnumeric);
            addParameter(p,'do_centroids', 1, @isnumeric);
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
         
            if inneraxis=='y'
                innerchan=2;
                outerchan=1;
            elseif inneraxis=='z'
                innerchan=3;
                outerchan=2;
            else
                innerchan=1;
                outerchan=2;
            end
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
            
            diff_data = zeros(numimgs,3);
           
            
            % scan dimensions
            if(flag_struct.flyscan)
                txtfid = fopen([datapath '/scan_' num2str(scanid) '.txt']);
                
                tline = fgetl(txtfid);
                
                innerpts = 0;
                innerdone=0;
                
               
                for ii=1:numimgs
                    %if(mod(ii,10)==0) waitbar(ii/numimgs); end
                    temp1=fscanf(txtfid,'%f ',51);
                    tempx = fscanf(txtfid,'%s',2);
                    temp2 = fscanf(txtfid,'%f ',3);
                    if(not(innerdone))
                        if ii==1;
                            tempin = temp2(innerchan);
                        else
                            innerpts=innerpts+1;
                            if((abs(temp2(innerchan)-tempin)<0.05)||ii==numimgs)
                                innerdone=1;
                            end
                        end
                    end
                    
                end
                
                if not(p.Results.innerpts)==0 
                    innerpts  = p.Results.innerpts;                   
                end
                    
                 outerpts=round(numimgs/innerpts);
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
                    
                    diff_data(ii,1) = sum(sum(ccd));
                    
                    if(~isempty(ROIinteg))
                        diff_data(ii,4) = sum(sum(ccd(ROIinteg(1):ROIinteg(2), ROIinteg(3):ROIinteg(4))));
                    end
                    
                    if(flag_struct.do_centroids)
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
                %clear innerpts
                
                %% EDIT INNER SCANS
                %{
                if(ismember(scanid,[45202])) %need to change when adding scans
                    innerpts = 20;
                elseif(ismember(scanid,[45728 45732:45740 45751 45752:45776])) %need to change when adding scans
                    innerpts = 60;
                elseif(ismember(scanid,[45745 45750 45781:45830 45915:45917, 45969:45970])) %need to change when adding scans
                    innerpts = 20;
                elseif(ismember(scanid,[45853:45902,45907:45910, 45971:45975])) %need to change when adding scans
                    innerpts = 32;
                elseif(ismember(scanid,[45913:45914]))
                    innerpts = 11;
                elseif(ismember(scanid,[45918:45957, 45976:46015]))
                    innerpts = 40;
                end%close(h);
                %}
                
                if exist('innerpts')==0
                    disp('did not specify inner points')
                    return
                end
                
                
               
            else %only used if doing dscan
                innerpts=100; %SET MANUALLY
                outerpts=100; %SET MANUALLY
                for ii=1:numimgs
                    imgin =double(merlimgs(:,:,ii));
                    hotmask=double(merlimgs(:,:,ii))<2000;
                    line1=sum(imgin.*hotmask,1);  % vertical
                    line2=sum(imgin.*hotmask,2);  % horizontal
                    diff_data(ii,1)=sum(sum(imgin.*hotmask));
                    for kk=1:size(line1,2)
                        diff_data(ii,2)=diff_data(ii,2)+kk*line1(kk)/diff_data(ii,1);
                    end
                    for kk=1:size(line2,1)
                        diff_data(ii,3)=diff_data(ii,3)+kk*line2(kk)/diff_data(ii,1);
                    end
                end
            end
            
            if(~isempty(ROIinteg))
                scandata = zeros(outerpts,innerpts,8);
            else
                scandata = zeros(outerpts, innerpts, 9);
            end
            
            % get linear data
            if(flag_struct.flyscan)
                frewind(txtfid);
                tline = fgetl(txtfid);
                for ii=1:outerpts
                    for jj=1:innerpts
                        temp1=fscanf(txtfid,'%f ',51+6); %used to be 60, used to be 51 March 6, 2018
                        tempx = fscanf(txtfid,'%s',2);
                        temp2 = fscanf(txtfid,'%f ',3);
                        scandata(ii,jj,1) = temp1(detchan);
                        scandata(ii,jj,4) = temp1(53); %this is correct - IC3
                        scandata(ii,jj,2) = temp2(outerchan);
                        scandata(ii,jj,3) = temp2(innerchan);
                        scandata(ii,jj,5) = diff_data((ii-1)*innerpts+jj,1);
                        scandata(ii,jj,6) = diff_data((ii-1)*innerpts+jj,2);
                        scandata(ii,jj,7) = diff_data((ii-1)*innerpts+jj,3);
                        scandata(ii,jj,8) = (ii-1)*innerpts+jj;
                        if(~isempty(ROIinteg))
                            scandata(ii,jj,9) = diff_data((ii-1)*innerpts+jj,4);
                        end
                    end
                end
            else
                for ii=1:outerpts
                    for jj=1:innerpts
                        scandata(ii,jj,2) = ii*0.1;
                        scandata(ii,jj,3) = jj*0.1;
                        scandata(ii,jj,4) = 1;
                        scandata(ii,jj,5) = diff_data((ii-1)*innerpts+jj,1);
                        scandata(ii,jj,6) = diff_data((ii-1)*innerpts+jj,2);
                        scandata(ii,jj,7) = diff_data((ii-1)*innerpts+jj,3);
                        scandata(ii,jj,8) = (ii-1)*innerpts+jj;
                    end
                end
            end
            
            %{
merlfilepath = [datapath '/scan2Dexport_' num2str(scanid) '.h5'];
merlhdfpath = ['/HXN_' num2str(scanid) '/primary/data/merlin1'];

%merl = h5read([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid) '/primary/data/merlin1']);
x = h5read([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid) '/primary/data/zpssx']);
y= h5read([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid) '/primary/data/zpssy']);
I0= h5read([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid) '/primary/data/sclr1_ch4']);
fluo = h5read([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid) '/primary/data/Det2_' detchan]);
%merl = squeeze(merl);
header = h5readatt([datapath '/scan2Dexport_' num2str(scanid) '.h5'], ['/HXN_' num2str(scanid)], 'start');
%pixx = size(merl,2); pixy = size(merl,1);
pixx = 515; pixy = 515;

k = strfind(header, '"shape":');
st = header{1}(k{1}:(k{1}+30));
st = strsplit(st);

xpts = str2num(st{2}(2:end-1));
ypts = str2num(st{3}(1:end-2));

scandata = zeros(ypts,xpts,4);

%ccdnames = cell(ypts,xpts);
%merl2 = zeros([xpts ypts size(merl(:,:,1))]);

cnt = 1;
for ii=1:ypts
    for jj=1:xpts
        
        merlind = jj+(ii-1)*xpts;
        scandata(ii,jj,3)=x(cnt);
        scandata(ii,jj,2)=y(cnt);
        scandata(ii,jj,1)=fluo(cnt);
        scandata(ii,jj,4)=I0(cnt);
        
        merlstruct(ii,jj).A = ['h5read(''' merlfilepath ''', ''' merlhdfpath ''', [1 1 1 ' num2str(merlind) '], [' num2str([pixx pixy]) ' 1 1],[1 1 1 1])'];
        
        %merl2(ii,jj,:,:) = merl(:,:,cnt);
        %ccdnames{ii}{jj}=[datapath '/Images/' num2str(scanid) '/' AA{index1}(end-47:end)];
        cnt = cnt + 1;
    end
end
            %}
            %%{
            if flag_struct.showmerlin==1
              hfig1 =  ND_display_data.showmerlin_function(scandata,pResults,pixx,pixy,2000)
            end
            
            if flag_struct.do_centroids==1
              hfig2 =  ND_display_data.show_centroid(scandata,pResults,pixx,pixy,2000)
            end
            
            %}
        
        
        end
        
         function [ fly2Dmaps,imapx,imapy,sumim ] = ThetaScan_film(datapath, scanid, detchan,varargin)
            
              p = inputParser; 
            
            
            addRequired(p,'datapath', @ischar);
            addRequired(p,'scanid');
            addRequired(p,'detchan', @isnumeric);
            addParameter(p,'thetalist',[]);
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
            
            pResults = p.Results;
            
            
            
            
            if isempty(thetalist)
                thetalist = 76.350+0.05*((1:numel(scanid))-1);
                %thetalist = 76.5+0.1*((1:numel(fly2Dscanlist))-1);
            end

           
            
            [dataout,imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(1),detchan,'innerpts',p.Results.innerpts,'showmerlin',p.Results.showmerlin);
            [datatrash,imgstrash] = ND_read_data.loadscan_HXN(datapath,scanid(1), 46,'innerpts',p.Results.innerpts,'showmerlin',p.Results.showmerlin); %% reads out the photo current;
            
            XRF0 = dataout(:,:,1)./dataout(:,:,4);
            PC0 = datatrash(:,:,1)./datatrash(:,:,4);
            

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
            for kk = 1:size(dataout,1)
                for ll = 1:size(dataout,2)
                    fly2Dmaps.ii(kk).jj(ll).im = double(imgsout(:,:,dataout(kk,ll,8)))./(numel(scanid));
                    fly2Dmaps.ii(kk).jj(ll).intensity(1) = dataout(kk,ll,5);
                    fly2Dmaps.ii(kk).jj(ll).SumInt = dataout(kk,ll,5);
                end
            end
            
            for ii=2:numel(scanid)
                %for ii=2:2
                [dataout imgsout] = ND_read_data.loadscan_HXN(datapath,scanid(ii),detchan,'innerpts',p.Results.innerpts,'showmerlin',p.Results.showmerlin);
                [datatrash imgstrash] =  ND_read_data.loadscan_HXN(datapath,scanid(ii),46,'innerpts',p.Results.innerpts,'showmerlin',p.Results.showmerlin); %% reads out the photo current;
                
                XRF1= dataout(:,:,1)./dataout(:,:,4); %normalize by monitor
                PC1 = datatrash(:,:,1)./datatrash(:,:,4); %photo current normalized by monitor;
               
                for jj = 1:size(chi2,1)
                    for kk = 1:size(chi2,2)
                        chi2(jj,kk) = sum(sum((XRF0 - circshift(XRF1,[jj-size(chi2,1)/2,kk-size(chi2,2)/2])).^2));
                        %{
                        if ismember(scanid(ii),[40376:40407])  % The Fluo sum is only recorded up to a certain number of points,
                            %these scans had too many points, so taking a
                            %subsection of the fluo data to do the chi fitting
                            tm = min([30,size(XRF1,1)]);
                            tmat = circshift(XRF1,[jj-size(chi2,1)/2,kk-size(chi2,2)/2]);
                            chi2(jj,kk) = sum(sum((XRF0(1:tm,:) - tmat(1:tm,:)).^2));
                        else
                            chi2(jj,kk) = sum(sum((XRF0 - circshift(XRF1,[jj-size(chi2,1)/2,kk-size(chi2,2)/2])).^2));
                        end
                        %}
                    end
                end
                [xcen,ycen,intemp] = find(chi2==min(min(chi2)));
                %figure(4);imagesc(chi2);axis image tight;colormap hot
                xshifts(ii) = xshifts(ii-1)+xcen-size(chi2,1)/2;
                yshifts(ii) = yshifts(ii-1)+ycen-size(chi2,2)/2;
           
                %{
                figure(3);imagesc(XRF1);axis image tight off;colormap hot;title(['Unshifted Theta: ' num2str(thetalist(ii))]);
                figure(2);imagesc(circshift(XRF1,[xshifts(ii),yshifts(ii)]));axis image tight off;colormap hot;title(['Theta: ' num2str(thetalist(ii)) ' x shift: ' num2str(yshifts(ii)) ' y shift: ' num2str(xshifts(ii))]);
                %}
                
                %adding in to see the shifts and confirm
                %saveas(gcf,['~/Desktop/test' num2str(thetalist(ii)) '.pdf'],'pdf')
                %
                
                fly2Dmaps.scan(ii).XRF = circshift(XRF1,[xshifts(ii),yshifts(ii)]);
                fly2Dmaps.scan(ii).PC = circshift(PC1, [xshifts(ii),yshifts(ii)]);
                fly2Dmaps.scan(ii).xshift = xshifts(ii);
                fly2Dmaps.scan(ii).yshift = yshifts(ii);
                fly2Dmaps.scan(ii).theta = thetalist(ii);
                ccdnums = circshift(dataout(:,:,8),[xshifts(ii),yshifts(ii)]);
                tempints = circshift(dataout(:,:,5),[xshifts(ii),yshifts(ii)]);
                
                for kk = 1:size(dataout,1)
                    for ll = 1:size(dataout,2)
                        fly2Dmaps.ii(kk).jj(ll).im = fly2Dmaps.ii(kk).jj(ll).im + double(imgsout(:,:,ccdnums(kk,ll)))./(numel(scanid));
                        fly2Dmaps.ii(kk).jj(ll).intensity(ii) = tempints(kk,ll);
                        fly2Dmaps.ii(kk).jj(ll).SumInt = fly2Dmaps.ii(kk).jj(ll).SumInt+tempints(kk,ll);
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
            imapx = dataout(1,:,3);
            imapy = dataout(:,1,2);
            sumim = tempim;
            
            if not(p.Results.plotflag)
                ND_display_data.display2Dmap(XRF0,'figNum',1,'figTitle',['Theta: ' num2str(thetalist(1))]);
                hfig =  ND_display_data.display2Dmap_toclick(tempim,'figNum',20,'Xval',imapx,'Yval',imapy,'figTitle','Sum Int.');
                ND_display_data.display2Dmap(tempxcen,'figNum',22,'Xval',imapx,'Yval',imapy,'figTitle','Row Centroids');
                ND_display_data.display2Dmap(tempycen,'figNum',23,'Xval',imapx,'Yval',imapy,'figTitle',['Column Centroids ']);
            end
            %{
            figure(20);hfig=imagesc(dataout(1,:,3),dataout(:,1,2),tempim);axis image tight;colormap hot
            pass2click.xaxis = [min(dataout(1,:,3)) max(dataout(1,:,3)) size(dataout(1,:,3),2)];
            pass2click.yaxis = [min(dataout(:,1,2)) max(dataout(:,1,2)) size(dataout(:,1,2),1)];
            datacursormode on;
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
            set(dcm_obj, 'UpdateFcn', @click4rock_film);
            %}
            
            
            
           
            %{
            figure(22);clf
            imagesc(dataout(1,:,3),dataout(:,1,2),tempxcen);axis image tight;colormap hot;title('Row Centroids');
            figure(23);clf
            imagesc(dataout(1,:,3),dataout(:,1,2),tempycen);axis image tight;colormap hot;title('Col Centroids');
%}            
            
            %}
            
            
            
            
        end
        
   
        
    end
end
