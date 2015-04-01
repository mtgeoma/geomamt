%**********************************************************************************
function [rho,rho_se,ph,ph_se] = ap_res(z,sig_s,sig_e,periods)
%ap_res(...) : computes app. res., phase, errors, given imped., cov.
%USAGE: [rho,rho_se,ph,ph_se] = ap_res(z,sig_s,sig_e,periods) ;
% Z = array of impedances (from Z_***** file)
% sig_s = inverse signal covariance matrix (from Z_****** file)
% sig_e = residual covariance matrix (from Z_****** file)
% periods = array of periods (sec)

rad_deg = 57.2958;
rho = (abs(z).^2)/5.;
temp = size(periods);
nbt = temp(2);
temp = size(z);
nchnbt = temp(2);
nche = nchnbt/nbt;

zerr = zeros(2,nchnbt);
for ib = 1:nbt
   i1 = 2*(ib-1) + 1;
   i2 = 2*ib;
   s = real(diag(sig_s(:,i1:i2)));
   i1 = nche*(ib-1) + 1;
   i2 = nche*ib;
   e = real(diag(sig_e(:,i1:i2)));
   zerr(:,i1:i2) = s*e';
end

for l =1:2
  for k = 1:nche
    rho(l,k:nche:nche*(nbt-1)+k) = rho(l,k:nche:nche*(nbt-1)+k).*periods;
    rho_se(l,k:nche:nche*(nbt-1)+k) = periods;
  end
end

rho_se = sqrt(zerr.*rho_se.*rho*2/5.);
ph = rad_deg*atan(imag(z)./real(z));
ph_se = rad_deg*sqrt(zerr/2)./abs(z);
