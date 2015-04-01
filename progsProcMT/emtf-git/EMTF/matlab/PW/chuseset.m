%  chuseset sets up the polarization and predicted channels to plot
%   by settin poluse and chuse

lk = get(gco,'UserData');
l = lk(1); k = lk(2);
if( k == 1 )
   poluse = l;
   set(chhand(1,l,k),'Value',1);
   set(chhand(1,3-l,k),'Value',0);
else
   chuse(l,k) = get(chhand(1,l,k),'Value');
end