%**********************************************************************************
function [mode] = get_mode(rho,ixy,nbt,icomp);
%USAGE:    [mode] = get_mode(rho,ixy,nbt,icomp);
%   returns the portion of the real array rho(2,nche*nbt)
%  corresponding to the mode defined by ixy,icomp  ---
%         for Hy mode ixy = 2, icomp = tm dipole #s
%         for Hx mode ixy = 1, icomp = te dipole #s
%    this works pretty much in general ... icomp
%    can have as many dipole numbers in it as desired

temp = size(rho);
nchnbt = temp(2);
nche = nchnbt/nbt;
itemp = icomp;
nmode = length(icomp);
mode = zeros(nbt,nmode);
for ib = 1:nbt
   mode(ib,:) = rho(ixy,icomp+(ib-1)*nche);
end
return
end
