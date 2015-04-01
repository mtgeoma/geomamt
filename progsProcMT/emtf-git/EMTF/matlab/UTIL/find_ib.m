%*********************************************************************
%  find_ib : given a list of period/freq. bands, 
%            and  a target period/freq (t0)
%            finds number of closest band
%
% Usage  [ib] = find_ib(Nbands,periods,t0)
%
function [ib] = find_ib(Nbands,periods,t0)
ib = 1;
distance = abs(log(periods(1)/t0));
for i = 2:Nbands
   test = abs(log(periods(i)/t0));
   if(distance > test );
      distance = test;
      ib = i;
   end
end
