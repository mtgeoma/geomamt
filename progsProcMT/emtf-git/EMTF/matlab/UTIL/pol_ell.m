function pol_ell(x,y,dxr,dyr,dxi,dyi,clr)

% This function plots polarization ellipse centered
% at (x(n),y(n)) corresponding to the complex vectors
%  ( dxr(n) + i* dxi(n) , dyr(n) + i*dyi(n)); 
%   n = 1:N ;  clr is the color of the line used for 
%  the polarization ellipse  
%  The routine calls ellipse.m to compute the ellipse
%  Scaling of the complex vector into the (x,y) space
%   must be done before calling this routine

hold on
N = length(x);
clrsym = [clr '-'];
for n = 1:N
  u = dxr(n)+i*dxi(n);
  v = dyr(n)+i*dyi(n);
  cuv = ellipse(u,v,x(n),y(n));
  plot(cuv(1,:),cuv(2,:),clrsym)
end
set(gca,'Box','on');
