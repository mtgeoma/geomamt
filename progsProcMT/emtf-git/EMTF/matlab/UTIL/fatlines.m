function fatlines(h_ax,thick)
%  makes all lines thicker in axis with handle h_ax
hh = get(h_ax,'Children');
for k=1:length(hh)
   if(get(hh(k),'Type') == 'line')
      set(hh(k),'LineWidth',thick)
   end
end
