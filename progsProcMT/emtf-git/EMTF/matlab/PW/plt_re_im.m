%***********************************************************************
%plt_re_im : PLOTS real and imaginary parts of a series of TF curves, with error bars
%   on log-linear scale
%  %  USAGE [hfig] = plt_re_im(T,tr,ti,tr_se,lims,c_title,pos_rect,labels,ifref)
%   returns figure handle in hfig
%   Inputs : T(nbt,1) = periods
%            tr(nbt,nplot), ti(nbt,nplot) = real and imaginary parts of TF
%            tr_se(nbt,nplot) = s.e. of real/imag parts (same)
%            lims(6) = plotting limts for T, tr, im
%            c_title = plot title (character string)
%            pos_rect = position of figure on screen
%            lables(nplot+nref,6) = strings to lable each curve with

function hfig = plt_tr_im(T,tr,ti,tr_se,lims,c_title,pos_rect,labels,ifref)

pwstmon;
%       number of sets of curves to plot
nplot = length(tr(1,:));
%nbt = nm(1);

%	Set up some figure parameters
sym_lin_width = 1.0;
line_styles = ['b-';'r-';'g-';'c-';'m-';'k-';'y-';'b-';'r-';'g-';'c-';'m-'];
symbol_styles = ['bo';'rx';'go';'cx';'mo';'kx';'yo';'bx';'ro';'gx';'co';'mx'];
paper_width = 6.5;
paper_height = 10;
width_n = 0.7;
height_tr_n = .35;
x0_n = .15;
y0_ti_n = .1;
y0_tr_n = .55;
 
rect_paper = [1.,1.,paper_width,paper_height];
rect_tr = [x0_n,y0_tr_n,width_n,height_tr_n];
rect_ti = [x0_n,y0_ti_n,width_n,height_tr_n];
hfig = figure('Position',pos_rect,'PaperPosition',rect_paper);

%       plot real parts
tr_axes = axes('position',rect_tr) ;
set(gca,'Xgrid','on');
for c = 1:nplot ;
   [xb,yb] = err_log(T',tr(:,c),tr_se(:,c),'XLOG',lims);
   lines = semilogx(xb,yb,line_styles(c,:),T,tr(:,c),symbol_styles(c,:));
   set(lines,'LineWidth',sym_lin_width);
   hold on ;
end ;

pltlabl;

hold off ;
title_pos_y = lims(3) + .2*(lims(4)-lims(3)) ;
title_pos_x = log(lims(1)) + .1*(log(lims(2)/lims(1))); ;
title_pos_x = exp(title_pos_x);
set(gca,'FontSize',12,'FontWeight','bold')
text(title_pos_x,title_pos_y,'Real Part','FontWeight','bold', ...
     'FontSize',12);
axis(lims(1:4));
title(c_title);
set(get(gco,'Title'),'FontSize',14,'FontWeight','bold')

%      plot imaginary part
ti_axes = axes('position',rect_ti);
for c = 1:nplot;
   [xb,yb] = err_log(T',ti(:,c),tr_se(:,c),'XLOG',lims);
   lines = semilogx(xb,yb,line_styles(c,:),T,ti(:,c),symbol_styles(c,:));
	set(lines,'LineWidth',sym_lin_width);
   hold on;
end;
axis(lims(1:4));
title_pos_y = lims(3) + .8* (lims(4)-lims(3));
set(gca,'FontSize',12,'FontWeight','bold')
text(title_pos_x,title_pos_y,'Imaginary Part','FontWeight','bold', ...
  'FontSize',12);

xlabel('Period (s)');
hold off;
%
%  Make the plot nice
%
set(tr_axes,'ygrid','on','xgrid','on');
set(tr_axes,'FontWeight','bold');
%phase_ticks = [-40,-30,-20,-10,0,10.,20,30.,40];
%set(ph_axes,'ygrid','on','xgrid','on','YTickMode','manual',...
%    'Ytick',phase_ticks);
set(ti_axes,'ygrid','on','xgrid','on')
set(ti_axes,'FontWeight','bold');
%h_menu = uimenu('Label','Print','Parent',hfig)
%         uimenu(h_menu,'Label','Printer','Callback','print -fgcf');
%         uimenu(h_menu,'Label','File','Callback',['print -fgcf ' uigetfile]);
