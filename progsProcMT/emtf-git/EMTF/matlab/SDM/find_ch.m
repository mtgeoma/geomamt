function find_ch(h_check);

% makes lists ind1 and ind2 of channels in each of the two
% groups: ind1 = those channels checked, ind2 = the rest

global ind1 ind2
nt = length(h_check);
n1 = 0 ; n2 = 0;
ind1 = []; ind2 = [];
for k=1:nt
   if( get(h_check(k),'value') == 1 )
      n1 = n1 + 1;
      ind1(n1) = k;
    else
      n2 = n2 + 1;
      ind2(n2) = k;
    end
end
