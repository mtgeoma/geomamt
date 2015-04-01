function chk_clr(h)
%  clears a figure, after making sure it exists
hwin = get(0,'Children');
if( sum(hwin == h) > 0 )
   delete(h);
end
