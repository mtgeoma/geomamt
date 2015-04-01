%  construct curve to plot polarization ellipse
function [cuv] = ellipse(u0,v0,x,y)
scale = .5;
ncuv = 30 ;
theta = (2*pi/ncuv)*[0:ncuv];
theta = exp(i*theta);
cuv(1,:) = .6*real(u0*theta)+x;
cuv(2,:) = .6*real(v0*theta)+y;
cuv(:,31)  = [x;y]; 
