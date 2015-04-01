%**********************************************************************************
function [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = imp_ap(Z,SIG_S,SIG_E,periods)
%imp_ap(...) : computes app. res., phase, errors, given imped., cov.
%USAGE:  [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = imp_ap(Z,SIG_S,SIG_E,periods);
% INPUT:
%   Z(4,:)  = array of 2x2 impedance matrices
%   SIG_S(4,:)  = inverse signal covariance matrix (2 H)
%   SIG_E(4,:)  = residual covariance matrix (2 E)
%   periods = array of periods (sec)

rad_deg = 57.2958;
temp = size(periods);
nbt = temp(2);
temp = size(Z);
nt = temp(2);
nimp = nt/nbt;

rxy = zeros(nbt,nimp);
ryx = zeros(nbt,nimp);
rxy_se = zeros(nbt,nimp);
ryx_se = zeros(nbt,nimp);
pxy = zeros(nbt,nimp);
pyx = zeros(nbt,nimp);
pxy_se = zeros(nbt,nimp);
pyx_se = zeros(nbt,nimp);

for k=1:nimp
   k1 = nbt*(k-1)+1;
   k2 = nbt*k;
   ryx(:,k) = abs(Z(3,[k1:k2])').^2;
   rxy(:,k) = abs(Z(2,[k1:k2])').^2;
   pyx(:,k) = rad_deg*atan((imag(Z(3,[k1:k2]))')./real(Z(3,[k1:k2])'));
   pxy(:,k) = rad_deg*atan((imag(Z(2,[k1:k2]))')./real(Z(2,[k1:k2])'));
%   compute error variance for off-diagonal impedances (R & I sep) ...
   ryx_se(:,k) = real(SIG_S(1,[k1:k2])').*real(SIG_E(4,[k1:k2])')/2;
   rxy_se(:,k) = real(SIG_S(4,[k1:k2])').*real(SIG_E(1,[k1:k2])')/2;
end

%  one SE for phase in degrees
pyx_se = rad_deg*sqrt(ryx_se./ryx);
pxy_se = rad_deg*sqrt(rxy_se./rxy);

for l = 1:nbt
  ryx(l,:) = ryx(l,:)*periods(l)/5. ;
  rxy(l,:) = rxy(l,:)*periods(l)/5. ;
  ryx_se(l,:) = sqrt(ryx_se(l,:).*ryx(l,:)*periods(l)*4/5.);
  rxy_se(l,:) = sqrt(rxy_se(l,:).*rxy(l,:)*periods(l)*4/5.);
end
