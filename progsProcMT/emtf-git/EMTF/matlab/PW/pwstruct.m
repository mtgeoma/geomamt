%  load in contents of Pw-file, put header info into pwhd structure
%   array TFs and error covariances into pw structure
%     Usage: [pw,pwhd] = pwstruct(cfile)

function [pw,pwhd] = pwstruct(cfile)

%  open file, read header
[fid,recl,nbt,nt,nsta,nsig,nch,ih,stcor,...
                 decl,chid,csta,sta,orient] = Pw_hd(cfile);

period = zeros(nbt,1);
nf = period;
tf = zeros(2,nt,nbt);
xxinv = zeros(2,2,nbt);
cov = zeros(nt,nt,nbt);
csta = [];
for k = 1:nsta
  for l=1:nch(k)
    csta = [csta [sta(1:3,k)]];
  end
end
% standardize case for first two characters in channel ids
chid(1,:) = upper(setstr(chid(1,:)));
chid(2,:) = lower(setstr(chid(2,:)));
ch_name = [];
i2 = 11
for k = 1:nt
   corient = num2str(fix(orient(1,k)));
   ctemp = [ chid(:,k)' ':' csta(:,k)' ':' blanks(4) ];
   i1 = i2 - length(corient) + 1;
   ctemp(i1:i2) = corient;
   ch_name = [ch_name ; ctemp];
end
   
for ib = 1:nbt
  [period(ib),nf(ib),tf(:,:,ib),xxinv(:,:,ib),cov(:,:,ib)] =...
               Pw_in(fid,recl,ib);
end
pw = struct('T',period,'nf',nf,'tf',tf,'xxinv',xxinv,'cov',cov);
pwhd = struct('nbt',nbt,'nt',nt,'nsta',nsta,'nsig',nsig,'nch',nch,'ih',ih,...
   'stcor',stcor,'decl',decl,'chid',chid,'sta',sta,...
   'orient',orient,'ch_name',ch_name);
