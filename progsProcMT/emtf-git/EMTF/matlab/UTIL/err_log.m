%***********************************************************************
function [xb,yb] = err_log(x,y,yerr,ll,lims)
%err_log : used for plotting error bars with a y-axis log scale
%  takes VECTORS x and y and outputs matrices (one row per
% data point) for plotting error bars
%    ll = 'XLOG' for log X axis

n = size(x);
xb = zeros(6,n(2));
yb = zeros(6,n(2));
barsize = .0075;
if ((ll == 'XLOG') | (ll == 'xlog'));
  dx = log(lims(2)/lims(1))*barsize;
  xb(3,:)  = log(x);
  xb(4,:)  = xb(3,:);
  xb(1,:) = xb(3,:) - dx;
  xb(2,:) = xb(3,:) + dx;
  xb(5,:) = xb(3,:) - dx;
  xb(6,:) = xb(3,:) + dx;
  xb = exp(xb);
else
  dx = (lims(2)-lims(1))*barsize;
  xb(3,:)  = x;
  xb(4,:)  = xb(3,:);
  xb(1,:) = xb(3,:) - dx;
  xb(2,:) = xb(3,:) + dx;
  xb(5,:) = xb(3,:) - dx;
  xb(6,:) = xb(3,:) + dx;
end
yb(1,:) = (y - yerr)';
yb(2,:) = (y - yerr)';
yb(3,:) = (y - yerr)';
yb(4,:) = (y + yerr)';
yb(5,:) = (y + yerr)';
yb(6,:) = (y + yerr)';
