%**********************************************************************
%  sing_s  : reconstructs approximation to sdm using dominant evec,
%
%  Usage: S = sing_s(u,ev)

function S = sing_s(u,ev)
S = u*diag(ev)*u';
return
end
