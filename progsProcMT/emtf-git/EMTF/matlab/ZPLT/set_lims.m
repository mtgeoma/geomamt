function [lims,orient] = set_lims(dir,periods,rho)
 
fid = fopen([dir '/limits.mat'],'r');
if(fid > 0 )
%   a limits file is found in the current directory; use this to define 
%   initial limits
   fclose(fid)
   eval( [ 'load ' dir '/limits.mat'] ) ;
elseif(nargin == 3 )
% use max/min limits of periods, rho to set limits
  x_min = min(periods);
  x_max = max(periods);
  xx_min = 10^(floor(log10(x_min)));
  if ((log10(x_min)-log10(xx_min)) < 0.15)
          xx_min = 10^(log10(xx_min)-0.3);
  end
  xx_max = 10^(ceil(log10(x_max)));
  if ((log10(xx_max)-log10(x_max)) < 0.15)
          xx_max = 10^(log10(xx_max)+0.3);
  end
  y_min = min(min(rho));
  y_max = max(max(rho));
  yy_min = 10^(floor(log10(y_min)));
  if ((log10(y_min)-log10(yy_min)) < 0.15)
          yy_min = 10^(log10(yy_min)-0.3);
  end
  yy_max = 10^(ceil(log10(y_max)));
  if ((log10(yy_max)-log10(y_max)) < 0.15)
          yy_max = 10^(log10(yy_max)+0.3);
  end
  lims = [xx_min, xx_max, yy_min, yy_max];
  orient = 0.;
else
% use some general limits
  lims = [.01,1000,1,10000];
  orient = 0.;
end
