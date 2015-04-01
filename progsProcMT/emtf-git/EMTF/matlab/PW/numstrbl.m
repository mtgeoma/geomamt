function [c] = numstrbl(n);
%converts an integer to a string, with zero set to blank;
if(n == 0 )
   c = ' ';
else
   c = num2str(n);
end
