%***********************************************************************
%plt_amp_ph : PLOTS a series of amplitude and phase curves, with error bars
%   on log-log or log-linear scale
%  lims = [xmax,xmin,ymax,ymin] are plot limnits; empty array => automatic
%  USAGE [hfig] = plt_amp_ph(T,amp,amp_se,ph,ph_se,lims,...
%				ll,c_title,pos_rect,labels,ifref)
%   returns figure handle in hfig
%   Inputs : T(nbt,1) = periods
%            amp(nbt,nplot) = amplitudes
%            amp_se(nbt,nplot) = s.e. of amplitudes
%            ph(nbt,nplot) = phases
%            ph_se(nbt,nplot) = s.e. of phases
%            lims(6) = plotting limts for T, amp, ph
%            ll = 'XLOG' to plot amplitudes on log-linear scale
%            c_title = plot title (character string)
%            pos_rect = position of figure on screen
%            lables(nplot+nref,6) = strings to label each curve with
%					(including reference channels)
%			    ifref = 1/0 array = 1 for reference labels

function hfig = plt_amp_ph(T,amp,amp_se,ph,ph_se,lims,...
   plt_type,c_title,pos_rect,labels,ifref)

pwstmon;
%       number of sets of curves to plot
nplot = length(amp(1,:));
%nbt = nm(1);

%	Set up some figure parameters
sym_lin_width = 1.0;
line_styles = ['b-';'r-';'g-';'c-';'m-';'k-';'y-';'b-';'r-';'g-';'c-';'m-'];
symbol_styles = ['bo';'rx';'go';'cx';'mo';'kx';'yo';'bx';'ro';'gx';'co';'mx'];
paper_width = 6.5;
paper_height = 10;
width_n = 0.7;
height_ph_n = .3;
height_amp_n = .4;
x0_n = .15;
y0_ph_n = .1;
y0_amp_n = .5;

rect_paper = [1.,1.,paper_width,paper_height];
rect_amp = [x0_n,y0_amp_n,width_n,height_amp_n];
rect_ph = [x0_n,y0_ph_n,width_n,height_ph_n];
hfig = figure('Position',pos_rect,'PaperPosition',rect_paper);

%       plot amplitudes
amp_axes = axes('position',rect_amp) ;
set(gca,'Xgrid','on');
for c = 1:nplot ;
   [xb,yb] = err_log(T',amp(:,c),amp_se(:,c),'XLOG',lims);
   if(plt_type == 'AMPH')
      lines = semilogx(xb,yb,line_styles(c,:),T,amp(:,c),symbol_styles(c,:));
   else
      lines = loglog(xb,yb,line_styles(c,:),T,amp(:,c),symbol_styles(c,:));
   end
   set(lines,'LineWidth',sym_lin_width);
   hold on ;
end ;

pltlabl;

hold off ;
title_pos_y = lims(3) + .2*(lims(4)-lims(3)) ;
title_pos_x = log(lims(1)) + .1*(log(lims(2)/lims(1))); ;
title_pos_x = exp(title_pos_x);
set(gca,'FontSize',12,'FontWeight','bold')
text(title_pos_x,title_pos_y,'Amplitude','FontWeight','bold', ...
     'FontSize',12);
axis(lims(1:4));
title(c_title);
set(get(gco,'Title'),...
   'FontSize',14,...
   'FontWeight','bold',...
   'Interpreter', 'none');

%      plot phases
ph_axes = axes('position',rect_ph);
lims_ph = lims([1 2 5 6 ]);
for c = 1:nplot;
   [xb,yb] = err_log(T',ph(:,c),ph_se(:,c),'XLOG',lims);
   lines = semilogx(xb,yb,line_styles(c,:),T,ph(:,c),symbol_styles(c,:));
	set(lines,'LineWidth',sym_lin_width);
   hold on;
end;
axis(lims_ph);
title_pos_y = lims_ph(3) + .8* (lims_ph(4)-lims_ph(3));
set(gca,'FontSize',12,'FontWeight','bold')
text(title_pos_x,title_pos_y,'Phase','FontWeight','bold', ...
  'FontSize',12);

xlabel('Period (s)');
ylabel('Degrees');
hold off;
%
%  Make the plot nice
%
set(amp_axes,'ygrid','on','xgrid','on');
set(amp_axes,'FontWeight','bold');
%phase_ticks = [-40,-30,-20,-10,0,10.,20,30.,40];
%set(ph_axes,'ygrid','on','xgrid','on','YTickMode','manual',...
%    'Ytick',phase_ticks);
set(ph_axes,'ygrid','on','xgrid','on')
set(ph_axes,'FontWeight','bold');
%h_menu = uimenu('Label','Print','Parent',hfig)
%         uimenu(h_menu,'Label','Printer','Callback','print -fgcf');
%         uimenu(h_menu,'Label','File','Callback',['print -fgcf ' uigetfile]);
