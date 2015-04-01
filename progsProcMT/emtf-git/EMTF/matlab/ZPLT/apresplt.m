%  this script changes into the survey directory (to save
%  effort in browsing for files), defines MKPLTIND
%  (a character string which tells how to choose which
%  frequencies to take from which sampling band for plotting
%  combined results from several sampling bands), 
%  and then runs z_mtem which does the actual plotting

%  define filter for listing Z-files
zid = '*.z*' ;

%   cd to survey directory
%cd c:\mt24\acq24\park1

%  define routine for making array pltind
%
%  EXAMPLE 1 : Parkfield 1997 array:
%          MKPLTIND = ['pltind = [ pltind; mkpltind(cfile)];'];
%
%  Example 2 : Plot all frequencies : this is the "default"
    MKPLTIND = ['pltind = [pltind; ones(nbt,1)];'];

%  run z_mtem
z_mtem
