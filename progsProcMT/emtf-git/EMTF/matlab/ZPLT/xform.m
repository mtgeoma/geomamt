%**********************************************************************************
%
%   takes TFs input from Z file, converts to fixed coordinate system
%    (x points in direction 0 used as reference for channel orientations)
%    rotates, and computes rh, phi + standard errors
%   theta defines rotation angle, ixy defines pairs of channels to rotate
%   Usage : [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
%                     xform(Z,Sig_s,Sig_e,periods,orient,ixy,theta);

function [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,ixy,theta);

%   Extract 2x2 impedance matrices for pairs of channels specified in
%   integer array ixy ... at the same time convert all matrices and vectors into
%   common (usually geographic) coordinate system (NOTE: x axis will be 0 degrees,
%   y-axis 90 ... relative to whatever direction is used to define orientations of
%   channels as specified in array "orient").
%   In general ixy will have 2 rows and M = # of E-setups columns
%   a separate impedance will be returned for each specified E-setup
%   in this example rows  1 and 2 of the TF (as output in the Z.**** file)
%   correspond to the E outputs from the first station, while rows 3 and 4
%   correspond to the second station.  Two impedances, with matrices needed for
%   calculation of errors will be output.

[ dum,Nch] = size(orient);
Nche = Nch - 2; 
%   z_to_imp makes the 2x2 impedances, changes coordinates
[Z2x2,SIG_S,SIG_E] = z_to_imp(Z,Sig_e,Sig_s,Nche,ixy,orient);

%ON output Z2x2 is a 4xM matrix where M = Nimp*Nbt ... Nbt = # of frequency bands
% estimated, Nimp = # of pairs of impedances = 2 for the example here.
% The order of components is: Z2x2(:,1:Nbt) are impedance components for 1st site,
%    Z2x2(:,Nbt+1:2*Nbt) components for second site, etc.  The order of the 4 impedance
%     components (i.e., the 4 rows of Z2x2) is: Zxx, Zxy, Zyx, Zyy
%   SIG_S and SIG_E are of the same size, and follow the same ordering conventions

%   TO ROTATE:
%     use M-file rot_z ... give desired rotation angle for x-axis in 
%     degrees clockwise,  output matrices correspond exactly to
%     input matrices, but now everything is in the new coordinate system
%fprintf(1,'About to enter rot_z')

[Z2x2R,SIG_SR,SIG_ER]=rot_z(Z2x2,SIG_S,SIG_E,theta);

%   TRANSFROM TO rho, phi
[ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = imp_ap(Z2x2R,SIG_SR,SIG_ER,periods);
