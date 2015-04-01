%**********************************************************************************
function y_smooth = smooth(y,nsm)
%smooth : uses spline interpolation to make plots look prettier
%USAGE:  y_smooth = smooth(y,nsm)
nm = size(y);
n = nm(1);
m = nm(2);
msm = nsm*(m-1)+1;
y_smooth = zeros(n,msm);
temp = 0:nsm:nsm*(m-1) ;
x = temp';
temp = 0:nsm*(m-1);
xi = temp';
for i=1:n
    yin=y(i,:)';
    yout = spline(x,yin,xi)  ;
    y_smooth(i,:) = yout' ;
end

