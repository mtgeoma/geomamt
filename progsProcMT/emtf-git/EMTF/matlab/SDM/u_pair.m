%  selects out and rotates pairs of data H/E data channels specified
%  in arrays Hp/Ep/Hz
%
%  Usage :  function [Uh,Ue,Uz] = u_pair(u,Hp,Ep,Hz,orient,decl,stcor,period,rho_ref)

function [Uh,Ue,Uz] = u_pair(u,Hp,Ep,Hz,orient,decl,stcor,period, ...
    rho_ref,snr_units)

nt = length(u);
temp = size(Hp);
nH = temp(1);
temp = size(Ep);
nE = temp(1);
temp = size(Hz);
nZ = temp(1);
%orient = pi*(orient - ones(2,1)*decl')/180;
orient = pi*orient(1,:)/180;
c = cos(orient);
s = sin(orient);
%escl is scaling factor to convert E into nT ... for 100 ohm-m
%  apparent resistivity this will make E and H exactly the same magnitude
%  (more generally use sqrt(period/(5*rho)) )
if(snr_units)
   escl = 1
else
  escl = sqrt(period/(5*rho_ref));
end

% H vectors go in Uh
Uh = zeros(4,nH)+i*zeros(4,nH);
Uh(1:2,:) = stcor(:,Hp(:,3));
for k=1:nH
   ROT = [ c(Hp(k,1)) c(Hp(k,2)); ...
         s(Hp(k,1)) s(Hp(k,2))];
   Uh(3:4,k) = ROT*u(Hp(k,1:2));
end
% E vectors go in Ue
Ue = zeros(4,nE)+i*zeros(4,nE);
Ue(1:2,:) = stcor(:,Ep(:,3));
for k=1:nE
   ROT = [ c(Ep(k,1)) c(Ep(k,2)); ...
         s(Ep(k,1)) s(Ep(k,2))];
   Ue(3:4,k) = escl*ROT*u(Ep(k,1:2));
end
% Hz components go in Uz
Uz = zeros(3,nZ)+i*zeros(3,nZ);
Uz(1:2,:) = stcor(:,Hz(:,2));
Uz(3,:) = u(Hz(:,1)).';
return
end
