kf = [ 0 0 1 ]';
ki = [ 0.2865 -0.1771 0.9416 ]';
% ROIxstart = 1;%50;
% ROIystart = 133;%81;
% ROIxsize = 404;%200;
% ROIysize = 380%200;
% Parameters of Cu(Ga,In)Se2
%approximately 30% Ga.
a_CuInSe2 = 5.814; c_CuInSe2 = 11.63;  %https://link.springer.com/article/10.1007/BF02654305-4;
a_CuGaSe2 = 5.60; c_CuGaSe2 = 11.0; %https://ac.els-cdn.com/0022369766901570/1-s2.0-0022369766901570-main.pdf?_tid=3919c8ac-cdfe-11e7-8f5b-00000aab0f27&acdnat=1511187922_4fa529a971815446c8e0ac247d6cc071
Ga_x = 0.3; In_x = 1 - Ga_x;
a_CuInGaSe2 = a_CuInSe2*In_x + a_CuGaSe2*Ga_x;  % using Vegard's law
c_CuInGaSe2 = c_CuInSe2*In_x + c_CuGaSe2*Ga_x;

d_111 = 3.7412e-4; %using aps sector 7 calculator, in microns

Ekev = 9.0;
detdist = 4.7268e5;
lam = etolambda(10400)*1e-4;
pixsize = 55.;
Ndet = 512; % pixels

del = -16.65;
gam = 10.65;
twoTheta = 20.997;

degperpix =0.0067;

