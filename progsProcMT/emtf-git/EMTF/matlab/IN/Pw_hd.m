%************************************************************************
%Usage:  [fid,recl,nbt,nt,nsta,nsig,nch,ih,stcor,decl,chid,csta,sta,orient] = Pw_hd(cfile);
%   Input:   cfile = file name
%   Returns: fid_uev = file id
%            recl = record length
%            nbt,nt,nsta,nsig = # of bands, components, stations, evecs
%            nch(nsta),ih(nsta+1),stcor(2,nsta),decl(nsta),sta(nsta),
%            orient(nsta)  = the usual
%
%*********************************************************************

function [fid,recl,nbt,nt,nsta,nsig,nch,ih,stcor,...
              decl,chid,csta,sta,orient]= Pw_hd(cfile)
%  'b' option in fopen assumes binary Pw-file was written on UNIX system
%  change to 'l' for files written on PC
fid = fopen(cfile,'r','b');
fseek(fid,4,'bof');
nt    = fread(fid,1,'long');
nsta  = fread(fid,1,'long');
nsig  = fread(fid,1,'long');
nbt   = fread(fid,1,'long');
for k=1:nsta
  fseek(fid,8,'cof');
  nch(k) = fread(fid,1,'long');
  ih(k)  = fread(fid,1,'long');
  stcor(1:2,k) = fread(fid,2,'float');
  decl(k)  = fread(fid,1,'float');
  sta(1:3,k)   = fread(fid,3,'char');
end
fseek(fid,4,'cof');
nbytes = fread(fid,1,'long');
chid_length = nbytes - 11;
for l = 1:nt
   orient(1:2,l) = fread(fid,2,'float');
   chid(1:chid_length,l) = fread(fid,chid_length,'char');
   csta(1:3,l) = fread(fid,3,'char');
   if(l < nt) fseek(fid,8,'cof'); end
end
%recl(1) = length of header block ...
recl(1) = 24+nsta*31+nt*(19+chid_length);
ns = nt*(nt+1)/2;
recl(2) = 64+8*(ns+2*nt);
