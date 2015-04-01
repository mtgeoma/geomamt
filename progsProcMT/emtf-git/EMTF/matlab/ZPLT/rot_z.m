%**********************************************************************************
function [Z2x2R,SIG_SR,SIG_ER]=rot_z(Z2x2,SIG_S,SIG_E,theta)
%  rotates impedances, signal, error covs into new coordinate
%   system 

%  theta is in degrees (positive is clockwise)

c = cos(theta*pi/180.);
s = sin(theta*pi/180.);

U = [ c*c , c*s , c*s , s*s ;
     -c*s , c*c , -s*s , c*s ;
     -c*s , -s*s , c*c , c*s ;
     s*s  , -c*s , -c*s , c*c ];

Z2x2R = U*Z2x2;
SIG_ER = U*SIG_E;
SIG_SR = U*SIG_S;
