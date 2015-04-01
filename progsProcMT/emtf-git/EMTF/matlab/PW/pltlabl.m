%   font size for labeling curves 
if(nplot <= 3 ) 
   fs = fs3;
   dy_n = .12
elseif (nplot <= 5)
   fs = fs5;
   dy_n = .1;
else
   fs = fs6;
   dy_n = .08;
end
xlab_n = 1.02;

nlab = length(labels(:,1));
y0_n = 1;
c1 = 0;
for c = 1:nlab
   
   if  ifref(c)
      y0_n = y0_n -1.3*dy_n;
		text('position',[xlab_n,y0_n],'Units','Normalized',...
      	'string','Ref. Ch.:','color','k',...
         'FontWeight','bold','FontSize',fs-1);
      y0_n = y0_n - .7*dy_n;
		text('position',[xlab_n,y0_n],'Units','Normalized',...
      	'string',labels(c,:),'color','k',...
         'FontWeight','bold','FontSize',fs);
      y0_n = y0_n - .1*dy_n;
		text('position',[xlab_n,y0_n],'Units','Normalized',...
      	'string',' ______ ','color','k',...
      	'FontWeight','bold','FontSize',fs-1,...
      	'VerticalAlignment','middle',...
         'Interpreter','none');
      y0_n = y0_n - .1*dy_n
		text('position',[xlab_n,y0_n],'Units','Normalized',...
      	'string',' ______ ','color','k',...
      	'FontWeight','bold','FontSize',fs-1,...
      	'VerticalAlignment','middle',...
         'Interpreter','none');
      y0_n = y0_n - .15*dy_n
   else
      y0_n = y0_n - dy_n;
      c1 = c1 + 1;
      text('position',[xlab_n,y0_n],'Units','Normalized',...
      	'string',labels(c,:),'color',line_styles(c1,1),...
         'FontWeight','bold','FontSize',fs);
   end
end

