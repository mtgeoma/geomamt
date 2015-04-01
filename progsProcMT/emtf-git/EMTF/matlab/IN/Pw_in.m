%***********************************************************************
% Pw_in reads in array TFs from Pw**** file for one band
%
% USAGE: [period,nf,tf,xxinv,cov] = Pw_in(fid,recl,ib);
%
%  INPUTS:
%    fid = file id for input Pw**** file 
%    recl(2) = array of header, data record lengths returned
%              by Pw_hd
%    ib =    frequency band desired
%
%  OUTPUTS:
%    period, nf period, # of data points used for estimate
%    tf(2,nt)   array TF
%    xxinv(2,2) inverse signal power matrix
%    cov(nt,nt) complex Hermitian residual covariance (full matrix)
%
%***********************************************************************

function [period,nf,tf,xxinv,cov] = Pw_in(fid,recl,ib)
nskip = recl(1) + (ib-1)*recl(2) + 4;
fseek(fid,nskip,'bof');
period = fread(fid,1,'float');
nf = fread(fid,1,'long');
nt = fread(fid,1,'long');
ns = fread(fid,1,'long');
fseek(fid,8,'cof');
tf = fread(fid,[4,nt],'float');
tf= tf(1:2:3,:)+i*tf(2:2:4,:);
xxinv = fread(fid,8,'float');
xxinv = xxinv(1:2:7)+i*xxinv(2:2:8);
xxinv = reshape(xxinv,2,2);
cov = fread(fid,2*ns,'float');
cov = cov(1:2:2*ns-1)+i*cov(2:2:2*ns);
%   convert cov to full storage mode ... save programming
%    hassles later
row = [];
col = [];
for k=1:nt
  row = [ row , k*ones(1,k) ]; col = [ col  , [1:k]];
end
S = sparse(row,col,cov) ;
S = full(S);  cov =  S + S' - diag(diag(S));
return
