%  ch_pair pairs off E and H channels at a single site
%  following default rules ...
%      i.e., for each station (as defined by the array nch)
%     the first Hx is sought and the first Hy ... if both are
%      found a pair is created (all other Hx/Hy in the "station"
%      are ignored); The same is done with Ex/Ey ... so there are
%      0-2 rotatable channel pairs identified for each station.
%      Other sorts of pairings which might be desired will require
%      other routines/direct user specification
% Usage: [Hp,Ep,Hz] = ch_pair(nch,chid);
%     Hp/Ep give pairs of H/E channels, Hz gives first Hz at
%        each station if any; for Hp/Ep(Hz) column 3(2) gives
%        the corresponding station #, 

function [Hp,Ep,Hz] = ch_pair(nch,chid);

nsta = length(nch);
ih = [1];
for ista=1:nsta
   ih = [ih,ih(ista)+nch(ista)];
end
Hp = [];
Ep = [];
Hz = [];
for ista=1:nsta
%   Hx/Hy pair
   Hxy = zeros(1,3);
   Hxy(3) = ista;
   for k=ih(ista+1)-1:-1:ih(ista)
      if(upper(setstr(chid(1:2,k)))' == 'HX') Hxy(1) = k; end
      if(upper(setstr(chid(1:2,k)))' == 'HY') Hxy(2) = k; end
   end
   if(min(Hxy) > 0 ) 
      Hp = [Hp ; Hxy ];
   end        
%   Ex/Ey pair
   Exy = zeros(1,3);
   Exy(3) = ista;
   for k=ih(ista+1)-1:-1:ih(ista)
      if(upper(setstr(chid(1:2,k)))' == 'EX') Exy(1) = k; end
      if(upper(setstr(chid(1:2,k)))' == 'EY') Exy(2) = k; end
   end
   if(min(Exy) > 0 ) 
      Ep = [Ep ; Exy ];
   end        
%  Hz 
   for k=ih(ista):ih(ista+1)-1
      if(upper(setstr(chid(1:2,k)))' == 'HZ') Hz = [ Hz ; [k,ista] ]; end
   end
end
return
