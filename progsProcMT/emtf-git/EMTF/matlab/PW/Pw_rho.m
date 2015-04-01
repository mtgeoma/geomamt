% compute apparent resistivities and phases from selected impedance 
%  components in arrays Vx, Vy   with error bars in Vx_sig, Vy_sig
%
%  Usage:  [rho,rho_se,ph,ph_se] = Pw_rho(V,V_sig,ind,periods);
function  [rho,rho_se,ph,ph_se] = Pw_rho(V,V_sig,ind,periods);
nbt = length(periods);
nc = length(ind);
T = periods*ones(1,nc);
rho = (abs(V(:,ind)).^2)/5.;
rho = rho.*T;
rho_se = real(V_sig(:,ind)).*sqrt(T.*rho*2/5.);
ph = (180/pi)*atan(imag(V(:,ind))./real(V(:,ind)));
ph_se = (180/pi)*real(V_sig(:,ind))./(abs(V(:,ind))*sqrt(2));
return
