%	MATLAB function replace_
%
%	function to replace the character '_' in a string by '\_'.
%	This is necesary if '_' is used in strings to be printed on a 
%	plot, because the '_' itself makes the following charcter
%	being a subscript.
%
%	USAGE: [cstring] = replace_(cstring);

function [cstring] = replace_(cstring)

lenc = length(cstring);
L = findstr(cstring,'_');
%  char(59) is a semicolon 
cbl = ones(1,length(cstring))*59;
cbl = char(cbl);
cbl(L) = '\';
cstring = reshape([cbl; cstring],1,2*length(cstring));
L = findstr(cstring,';');
LL = ones(1,length(cstring));
LL(L) = 0;
LL = logical(LL);
cstring = cstring(LL);