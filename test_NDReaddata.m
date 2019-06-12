
% test of the ND_read_data.importafile to read any txt file with headers,
% appropriate number of columns and automatic identification of the Cu
% chanel

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));


clear all;close all;

datapath = './Hruszkewycz_2018Q1/Data';

%scanid = 45417;
%scanid = 45193;
%scanid = 45285;
scanid = 45418;

filename = [datapath '/scan_' num2str(scanid) '.txt'];

innerpts = 41;
outerpts = 21;
innerpts_zeropad = 0;
outerpts_zeropad = 0;
startRow = 1;
endRow = innerpts*outerpts;
XRFchan = 'Det1_Cu'; % fluorescence channel
XBICchan = 'sclr1_ch3'; % photoluminescence channel
prefix = {'seq','Det','alive','dead','elapsed_time','scaler_alive','sclr','time','xspress','zpss'};

[scandata,merlimgs,hfig1,hfig2] = ND_read_data.loadscan_HXN(datapath,scanid,XRFchan,'showmerlin',1,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'prefix',prefix,'do_padding',0);
[scandata_XBIC] = ND_read_data.loadscan_HXN(datapath,scanid,XBICchan,'showmerlin',1,'innerpts',innerpts,'outerpts',outerpts,'innerpts_zeropad',innerpts_zeropad,'outerpts_zeropad',outerpts_zeropad,'prefix',prefix,'do_padding',0);


%{
%scan_test = ND_read_data.importfile(filename,startRow , endRow, 'prefix',prefix );

% test reading the headers:
 
%fileID = fopen(filename,'r');
            
%[var_struct, variable_names_cell] =  ND_read_data.prepare_header(fileID,'prefix',prefix);

% test getLinearData

%}

%{
txtfid = fopen(filename);
tline = fgetl(txtfid);
frewind(txtfid);
tline = fgetl(txtfid);
for ii=1:outerpts
    for jj=1:innerpts
        temp1=fscanf(txtfid,'%f ',51+6); %used to be 60, used to be 51 March 6, 2018
        tempx = fscanf(txtfid,'%s',2);
        temp2 = fscanf(txtfid,'%f ',3);
    end
end
%}
 