%  this m-file sets up a 0/1 array indicating which
%  frequencies should be plotted for each band,
%  after merging with one or more other bands.
%
%  The idea behind this routine is that from the file name (cfile)
%  the routine can figure out which band the file corresponds to,
%  and then assign a default set of 0/1 indicators for that band.  
%
%  This particular script is used by apresplt.m to define the 
%   character string  MKPLTIND = ['pltind = mkpltind(cfile);'];
%   Then, in the main plotting routine z_mtem  :
%                   eval(MKPLTIND)
%  is executed for each file to find which frequencies should be displayed
%
%  This routine is all very hard-wired for a specific
%  band setup, particular sampling rates, file naming conventions, etc.
%  (in particular, for the Parkfield 97 EMAP survey ; more details below)
%
%  The idea is that you should tailor the definition of string MKPLTIN
%   to your particular application  (e.g., to plot all frequencies define
%          MKPLTIND = ['pltind = ones(nbt,1);'];
%   in z_mtem;  For other more complicated applications a new M-file 
%   like this one would be required.

function pltind = mkpltind(cfile)

%  for this example there are three bands:  mid (M = 960 hz)
%  low (L=120 hz) and very low (V=3.125 hz)  ;  there are estimates at
%   30 frequencies for each sampling band;  The range of frequencies
%  to plot are given by limits M1:M2, etc.

M1 = 2; M2 = 15;
L1 = 5;L2 = 16;
V1 = 1;V2 = 30;
pltind_M = zeros(30,1);
pltind_L = pltind_M;
pltind_V = pltind_M;

pltind_M(M1:M2) = ones(M2-M1+1,1);
pltind_L(L1:L2) = ones(L2-L1+1,1);
pltind_V(V1:V2) = ones(V2-V1+1,1);

%  also skip band M(3) ... (60 hz?)
pltind_M(3) = 0;

%  7th character of file name indicates band for Z-files

c = upper(cfile(length(cfile)-7));
eval(['pltind = pltind_' c ';']);
