%************************************************************************
%sdm_init : %Initializes S# file containing SDMs from multiple station program
%Usage:  [fid_uev,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,decl,sta,chid,csta,orient,periods] = sdm_init(cfile);
%   Input:   cfile = file name
%   Returns: fid_uev = file id
%            irecl = record length
%            nbt,nt,nsta,nsig = # of bands, components, stations, evecs
%            nch(nsta),ih(nsta+1),stcor(2,nsta),decl(nsta),sta(nsta),
%            orient(nsta)  = the usual
%            periods(nbt)  = periods

function [fid,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,...
         decl,sta,chid,csta,orient,periods] = sdm_init(cfile)
%  'b' option in fopen assumes binary Pw-file was written on UNIX system
%  change to 'l' for files written on PC
fid = fopen(cfile,'r','b');
irecl = fread(fid,1,'long');nbt   = fread(fid,1,'long');
nt    = fread(fid,1,'long');
nsta  = fread(fid,1,'long');
nsig  = fread(fid,1,'long');
csta = [];
for k=1:nsta
  nch(k) = fread(fid,1,'long');
  ih(k)  = fread(fid,1,'long');
  stcor(1:2,k) = fread(fid,2,'float');
  decl(k)  = fread(fid,1,'float');
  sta(1:3,k)   = fread(fid,3,'char');
  for l=1:nch(k)
    csta = [csta [sta(1:3,k)]];
  end
end
%%%  This is for NEW sdm files ... 6 character channel IDs
for l = 1:nt
   orient(1:2,l) = fread(fid,2,'float');
   chid(1:6,l) = fread(fid,6,'char');
end
%  Check to see if  orientations are reasonable:
if any( orient(2,:) ~= 0 )  
%    assume that this is an old file ... get rid of this if you
%   have no old sdm files, and channels with non=zero tilts
fseek(fid,-nt*14,'cof')
clear chid
%%%  This is for OLD sdm files ... 2 character channel IDs
   for l = 1:nt
      orient(1:2,l) = fread(fid,2,'float');
      chid(1:2,l) = fread(fid,2,'char')
   end
end

% go through file and get all of the periods ...
for ib = 1:nbt
   status = fseek(fid,irecl*ib,'bof');
   periods(ib) = fread(fid,1,'float') ;
end

%  File header written by this:
%       write(ivunits(ix),rec=1) irecl,nbt,nt,nstau,nsig,
%     &      (nchu(k),ih(k),stcor(1,k),stcor(2,k),decl(k),sta(k),
%     &      k=1,nstau),(orient(1,l),orient(2,l),l=1,nt)

%  Record for each band written by this
%      subroutine wrt_uev(iouev,irec,u,ev,var,nt,nsig,nf,period)
%      complex u(nt,nsig)
%      real ev(nsig),var(nt),period
%      integer nf
%      write(iouev,rec=irec) period,nf,ev,var,u
%      return
%      end
