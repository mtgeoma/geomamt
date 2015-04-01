% finds band which is closest to frequncy where slider is released,
%  then displays in small box above plot button

T = get(h_slider,'Value');
T = 10^T;
ib = find_ib(nbt,periods,T);
if(exist('h_ib'))
   delete(h_ib);
   clear h_ib;
end
h_ib = axes('Position',[.05,.34,.05,.04],'XTickLabelMode','manual', ...
          'YtickLabelMode','manual',...
            'box','on','Xlim',[0,1],'Ylim',[0,1]);
text('Parent',h_ib,'Position',[.18,.32],'String',num2str(ib),...
     'FontWeight','bold');
