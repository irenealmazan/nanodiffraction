% This script saves the data in excel, csv format.

clear all; close all;

addpath(genpath('./nanodiff_functions'));
addpath(genpath(['/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/m_scripts']));

datapath = './Hruszkewycz_2018Q1/Data';

Grain_array = [1 7 9 11 12];


for jj = 1:numel(Grain_array)
    
    load(['./results_paper/Figure5_grain' num2str(Grain_array(jj)) '_struct.mat'],'grain_struct');
    
    filename = ['results_paper/grain' num2str(Grain_array(jj)) '.xlsx'];
    
    
    varnames =  {['Contour'] ['Lattice_constant'] ['Lattice_constant_bottom'] ['Lattice_constant_up']};
    T = table(grain_struct.contour_values_down',grain_struct.dspace_distr',grain_struct.dspace_distr'-grain_struct.sigma_dspace_distr'./2,grain_struct.dspace_distr'+grain_struct.sigma_dspace_distr'./2,'VariableNames',varnames);
    
    writetable(T,filename,'Sheet',1,'Range','D1');
    
    varnames =  {['Contour'] ['Tilt'] ['Tilt_bottom'] ['Tilt_up']};
    T = table(grain_struct.contour_values_down',grain_struct.tilt_tot_distr',grain_struct.tilt_tot_distr'-grain_struct.sigma_tilt_tot_distr'./2,grain_struct.tilt_tot_distr'+grain_struct.sigma_tilt_tot_distr'./2,'VariableNames',varnames);
      
    writetable(T,filename,'Sheet',2,'Range','D1');
    
    varnames = {['theta'] ['intensity']};
    T = table(grain_struct.thetalist',grain_struct.rock_curve','VariableNames',varnames);
    writetable(T,filename,'Sheet',3,'Range','D1');


end


