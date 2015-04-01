%************************************************************************
%sn_in : reads in a SN_**** file output by multmtrn
% USAGE:  [ndf,E,S,NI] = sn_in(cfile);
%    cfile = name for SN****** file to read
%    ndf(nbt) =  # of data vectors used in each band
%    E(nt+1,nbt)  = eigenvalues in noise units
%    S(nt+1,nbt)  = signal power array
%    NI(nt+1,nbt)  = incoherent power array
%
%*********************************************************************

function [ndf,E,S,NI] = sn_in(cfile)
fid = fopen(cfile,'r');
nt = fscanf(fid,'%d',1);
nt1 = nt+1;
nbt = fscanf(fid,'%d',1);
ndf = fscanf(fid,'%d',nbt);
title = fgets(fid);
title = fgets(fid);
E = fscanf(fid,'%e',[nt1,nbt]);
title = fgets(fid);
title = fgets(fid);
S = fscanf(fid,'%e',[nt1,nbt]);
title = fgets(fid);
title = fgets(fid);
NI = fscanf(fid,'%e',[nt1,nbt]);
end
