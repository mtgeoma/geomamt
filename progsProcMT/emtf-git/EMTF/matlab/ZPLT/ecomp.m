%   ecomp sorts out and pairs off Ex, Ey components for construction
%    of conventional impedance tensors
%
%   USAGE: [xy,yx,hz,xypairs,DipoleSetup] = ecomp(nche,chid);
function [xy,yx,hz,xypairs,DipoleSetup] = ecomp(nche,chid)

nch = nche+2;

%   select out appropriate components
%   x/y mode (call it TE) ... corresponds to second H TF
xypairs = [];
xy = []; yx = []; hz = [];
for k = 3:nch
   if(upper(setstr(chid(k,1:2))) == 'EX')
      xy = [xy  k-2];
   elseif(upper(setstr(chid(k,1:2))) == 'EY')
      yx = [yx k-2];
   elseif(upper(setstr(chid(k,1))) == 'H')
      hz = [hz k-2 ];
   end
end
nxy = length(xy); nyx = length(yx);

if min([nxy nyx]) == 0
   %  all dipoles in one line ... for now just issue error msg 
   DipoleSetup = 'EMAP ';
   fprintf('Error: at least one of Ex/Ey not found in this file\n')
elseif (nxy ==1) & (nyx == 1)
   % conventional MT
   DipoleSetup = 'MT   ';
   xypairs = [ xy ; yx]; 
else
   % multiple dipoles in one direction ... and at least some cross-dipoles
   %  "tensor emap"
   DipoleSetup = 'TEMAP';
   %  find "nearest" pairing of available Ex/Ey components
   dist = abs(ones(nyx,1)*xy-yx'*ones(1,nxy));
   if nxy >= nyx
      mindist = min(dist);
      for k = 1:nxy
         for l = 1:nyx
            if dist(l,k) == mindist(k)
               yxu(k) = yx(l);
            end
         end
      end
      xypairs = [xy; yxu];
   else
      xyu = zeros(1,nyx);
      mindist = min(dist');
      for k = 1:nyx
         for l = 1:nxy
            if dist(k,l) == mindist(k)
               xyu(k) = xy(l);
            end
         end
      end
      xypairs = [xyu; yx];
   end
end
