%************************************************************************
%Z_in(cfile) : Reads in One Z_***** file  (TFs and signal/error cov.
%USAGE:  [z,sig_s,sig_e,periods,ndf,stcde,orient,nch,nche,nbt] = Z_in(cfile);
% nch = total # of channels ; nche = nch-2 = # of predicted channels ; nbt = # of bands
%   (NOTE: First two channels are always the "predictors"
% z(2,nche*nbt) = complex TFs
%   NOTE: Z(1,1:nche) corresponds to response to Hx sources for first band,
%         Z(2,1:nche) is Hy for for first band,
%         Z(1,nche+1:2*nche) corresponds to response to Hx sources for second band,
%         Z(2,nche+1:2*nche) is Hy for for second band,  ETC. ...
% sig_s(2,2*nbt) = complex inverse signal covariance
% sig_e(nche,nche*nbt) = complex residual error covariance
% stdec(3) = station coordinates, declination
% periods(nbt) = periods in seconds
% orient(2,nch) = orientation (degrees E of geomagnetic N) for each channel

function [z,sig_s,sig_e,periods,ndf,stdec,orient,nch,nche,nbt,chid,csta] = Z_in(cfile)

fid = fopen(cfile,'r');
if(fid < 0) 
   'error: cannot open file ',cfile
end
% skip text on first four lines
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
% skip text in coordinate line
cline = fgets(fid);
stdec = sscanf(cline,'coordinate %f %f declination %f');
cline = fgets(fid);
nchnbt = sscanf(cline,'number of channels %d number of frequencies %d');
nch = nchnbt(1);
nche = nch-2;
nbt = nchnbt(2);
cline = fgets(fid);
for k=1:nch
   idum = fscanf(fid,'%d',1);
   orient(1:2,k) = fscanf(fid,'%f',2);
   csta(k,1:3) = fscanf(fid,'%1s',3);
   chid(k,1:2) = fscanf(fid,'%1s',2);
   cline = fgets(fid);
end
cline = fgets(fid);
z = zeros(2,nche*nbt) + i*zeros(2,nche*nbt);
sig_e = zeros(nche,nche*nbt) + i*zeros(nche,nche*nbt);
sig_s = zeros(2,2*nbt) + i*zeros(2,2*nbt);
ndf = zeros(nbt,1);
for ib = 1:nbt
  cline = fgets(fid);
  periods(ib) = sscanf(cline,'period : %f');
  cline = fgets(fid);
  ndf(ib) = sscanf(cline,'number of data point %d');
  k1 = nche*(ib-1) + 1;
  k2 = nche*ib;
  cline = fgets(fid);
  ztemp = fscanf(fid,'%e',[4,nche]);
  z(1:2,k1:k2) = ztemp(1:2:3,:)+i*ztemp(2:2:4,:);
  chead = fgets(fid);
  chead = fgets(fid);
  stemp = fscanf(fid,'%e',[2,3]);
  ncht = 2;
  for k = 1:ncht
    for l = 1:ncht
      if(l < k ) 
        kl = (k*(k-1))/2+l;
        sig_s(k,2*(ib-1)+l) = stemp(1,kl)+i*stemp(2,kl);
      else
        kl = (l*(l-1))/2+k;
        sig_s(k,2*(ib-1)+l) = stemp(1,kl)-i*stemp(2,kl);
      end
    end
  end
  chead = fgets(fid);
  chead = fgets(fid);
  nse = (nche*(nche+1))/2;
  stemp = fscanf(fid,'%e',[2,nse]);
  for k = 1:nche
    for l = 1:nche
      if(l < k ) 
        kl = (k*(k-1))/2+l;
        sig_e(k,nche*(ib-1)+l) = stemp(1,kl)+i*stemp(2,kl);
      else
        kl = (l*(l-1))/2+k;
        sig_e(k,nche*(ib-1)+l) = stemp(1,kl)-i*stemp(2,kl);
      end
    end
  end
  cline = fgets(fid);
end  
