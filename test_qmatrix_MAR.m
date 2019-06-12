%%% This script tests the calibration of the reciprocal space for the tilt
%%% and the strain field calculation

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';


del = -16.1;
gam = 13.1;%10.65;
twoTheta = 21.2;
detdist = 0.35; % in meters

eval('Init_parameters');


angles_struct = ND_analysis.qmatrix_MAR(Ekev,del,gam,0,0,detdist*1000,Ndet,pixsize*1000,1);

matrix_ccd_to_lab = [angles_struct.ivec' angles_struct.jvec' angles_struct.Detcen'];

matrix_lab_to_ccd = inv(matrix_ccd_to_lab);

for jj = 1: size(angles_struct.qmat,1)
    for kk = 1: size(angles_struct.qmat,2)
        
        qmat_singlepix_lab = [angles_struct.qmat(jj,kk,1) angles_struct.qmat(jj,kk,2) angles_struct.qmat(jj,kk,3)];
        qmat_singlepix_ccd = matrix_ccd_to_lab\qmat_singlepix_lab';%matrix_lab_to_ccd*qmat_singlepix_lab';
        
        qmat_ccd(jj,kk,1) = qmat_singlepix_ccd(1);
        qmat_ccd(jj,kk,2) = qmat_singlepix_ccd(2);
        qmat_ccd(jj,kk,3) = qmat_singlepix_ccd(3);
        qmat_ccd(jj,kk,4) = sqrt(qmat_singlepix_ccd(1)^2+qmat_singlepix_ccd(2)^2+qmat_singlepix_ccd(3)^2);
    end
end


scrsz = get(0,'ScreenSize');  % Plot size parameters
imd = scrsz(3)/6;
imb = scrsz(4)-2*imd;

figure(1);
set(1,'Position',[2*imd imb imd imd]);
clf;
imagesc(qmat_ccd(:,:,1));axis image;colorbar;colormap hot;title('magnitude qx per pixel in ccd frame');

figure(2);
set(2,'Position',[2*imd imb imd imd]);
clf;
imagesc(qmat_ccd(:,:,2));axis image;colorbar;colormap hot;title('magnitude qy per pixel in ccd frame');

figure(3);
set(3,'Position',[2*imd imb imd imd]);
clf;imagesc(qmat_ccd(:,:,3));axis image;colorbar;colormap hot;title('magnitude qz per pixel in ccd frame');

figure(4);
set(4,'Position',[2*imd imb imd imd]);
clf;imagesc(qmat_ccd(:,:,4));axis image;colorbar;colormap hot;title('magnitude q per pixel in ccd frame');

