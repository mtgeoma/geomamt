%  using two "array TF" vectors in Pw.tf(2,nt)
%  computes TFs for all components relative to the two reference
%  components at a specified site. (For all frequncies
%   as opposed to ref_site which just does this for a single freq.)
%  Also outputs one standard error for these linear combinations
%  using full error covariance info as input in xxinv, cov
%
%  Usage: [v,sig_v] = pwrfsite(Pw,ref);
%
%   ref(2) = integers giving component numbers of reference channels
%    Pw = structure containing TFs and error covariance matrices
%   Pw.tf(2,nt,nbt) = TFs ouput in Pw**** file
%   Pw.xxinv(2,2,nbt) = inverse signal power matrix
%   Pw.cov(nt,nt.nbt)  = covariance of residuals (all components)

function [v,sig_v] = pwrfsite(Pw,ref)
nm = size(Pw.tf); nt = nm(2); nbt = nm(3);

v = zeros(2,nt,nbt) + zeros(2,nt,nbt);
sig_v = zeros(2,nt,nbt) + zeros(2,nt,nbt);
for ib = 1:nbt
% compute TFs relative to reference channels
  B = Pw.tf(:,ref,ib); Binv = inv(B); 
  v(:,:,ib) = Binv*Pw.tf(:,:,ib);

% error calculation
  M1 = Binv*Pw.xxinv(:,:,ib)*Binv';

% term corresponding to noise in predicted channels
  err = diag(Pw.cov(:,:,ib))';

%  now term corresponding to noise in reference channels
  err = err + sum(v(:,:,ib).*(Pw.cov(ref,ref,ib)*conj(v(:,:,ib))));
  sig_v(:,:,ib) = sqrt(real([M1(1,1)*err;M1(2,2)*err]));
  sig_v(:,ref,ib) = zeros(2,2);
end

return
