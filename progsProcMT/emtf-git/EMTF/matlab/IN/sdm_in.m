%************************************************************************
%uev_in : Reads in sdm etc. for a single band (S0**** file)
% Usage: [period,nf,var,S] = sdm_in(fid_uev,nt,ib,irecl)
% Inputs:
%     fid_uev = uev file id
%     nt, = number of components
%     ib = period band sought
%     irecl = band record length
% Outputs:
%     period = actual period
%     nf = # of data vectors
%     ev = eigenvalues
%     var = error variances
%     u = eigenvectors

function [period,nf,var,S] =uev_in(fid,nt,ib,irecl)

status = fseek(fid,ib*irecl,'bof');
period = fread(fid,1,'float');
nf = fread(fid,1,'long');
ns = (nt*(nt+1))/2;
var = fread(fid,nt,'float');
s = fread(fid,[2,ns],'float');
s = s(1,:)+i*s(2,:);
%   convert S to full storage mode ... save programming
%    hassles later
col = [];
row = [];
for k=1:nt
  row = [ row , k*ones(1,k) ]; col = [ col  , [1:k]];
end
S = sparse(row,col,s) ;
S = full(S);  S =  S + S' - diag(diag(S));

%  Record for each band written by this
%      subroutine wrt_uev(iouev,irec,s,ns,var,nt,nf,period)
%      complex s(ns)
%      real var(nt),period
%      integer nf
%      write(iouev,rec=irec) period,nf,var,s
%      return
%      end
