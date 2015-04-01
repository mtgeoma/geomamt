%***********************************************************************
%pltrhom : PLOTS a series of rho curves, with error bars on log-log
%  this version is a variant on plot_rho.m that allows for all arrays
%   to be divided into a series of bands ... so high/mid/low bands can
%   be plotted

function [rho_axes,ph_axes] = ...
  pltrhom(NBT,pltind,periods,rho,rho_err,ph,ph_err,lims,c_title,hfig)

globinc
%  global sym_line_width line_styles symbol_styles...
%       rect rect_paper rect_rho rect_ph
%

figure(hfig);
%  NB is number of "sampling bands"
NB = length(NBT);
indsB = zeros(2,NB);
indsB(:,1) = [1;NBT(1)];
for b=2:NB
  indsB(1,b) = indsB(2,b-1)+1;
  indsB(2,b) = indsB(1,b) + NBT(b) - 1;
end
  
nm = size(rho);
nperiodsT = nm(1);
ncurves = nm(2);
rho_axes = axes('position',rect_rho) ;
set(gca,'Xgrid','on');

xt1 = ceil(log10(lims(1))); xt2 = floor(log10(lims(2)));
xt2 = max(xt1,xt2);
xt = 10.^[xt1:1:xt2];

for b = 1:NB
  nperiods = NBT(b);
  i1 = indsB(1,b);i2 = indsB(2,b);
  if(2*fix(b/2) == b) 
     c0 = Nsym;
  else
     c0  = 0;
  end
  nper = length(periods);
  ind = zeros(nper,1);
  ind(i1:i2) = ones(i2-i1+1,1);
  ind = find(pltind(1:nper) & ind);
  for c = 1:ncurves ;
    [xb,yb] = err_log(periods(ind),rho(ind,c),rho_err(ind,c),'XLOG',lims);
        lines = loglog(xb,yb,line_styles(c0+c,:),periods(ind),rho(ind,c),...
                    symbol_styles(c0+c,:));
        set(lines,'LineWidth',sym_line_width);
    hold on ;
  end;
end ;
hold off ;
title_pos_y = log(lims(3)) + .2*(log(lims(4)/lims(3))) ;
title_pos_y = exp(title_pos_y) ;
title_pos_x = log(lims(1)) + .1*(log(lims(2)/lims(1))); 
title_pos_x = exp(title_pos_x);
text(title_pos_x,title_pos_y,'\rho_a','FontSize',14);
set(gca,'FontSize',11,...
	'FontWeight','bold',...
	'Xlim',lims(1:2),...
	'Ylim',lims(3:4));
ylabel('\Omega m');
title(c_title);
set(get(gca,'title'),'FontWeight','bold')
%
ph_axes = axes('position',rect_ph);
lims_ph = lims;
lims_ph(3) = 0.0;
lims_ph(4) = 90.0;
for b = 1:NB
  nperiods = NBT(b);
  i1 = indsB(1,b);i2 = indsB(2,b);
  if(2*fix(b/2) == b) 
     c0 = Nsym;
  else
     c0  = 0;
  end
  ind = zeros(nper,1);
  ind(i1:i2) = ones(i2-i1+1,1);
  ind = find(pltind(1:nper) & ind);
  for c = 1:ncurves;
    [xb,yb] = err_log(periods(ind),ph(ind,c),ph_err(ind,c),'XLOG',lims);
    lines = semilogx(xb,yb,line_styles(c0+c,:),periods(ind),ph(ind,c),...
              symbol_styles(c0+c,:));
         set(lines,'LineWidth',sym_line_width);
    hold on;
  end
end;
axis(lims_ph);
title_pos_y = lims_ph(3) + .8* (lims_ph(4)-lims_ph(3));
text(title_pos_x,title_pos_y,'\phi','FontSize',14);
set(gca,'FontWeight','bold',...
	'FontSize',11)
xlabel('Period (s)');
ylabel('Degrees');
end;
hold off;
%
%  Make the plot nice
%
set(rho_axes,'ygrid','on','xgrid','on');
set(rho_axes,'FontWeight','bold');
phase_ticks = [0,15.,30.,45.,60.,75.,90.];
set(ph_axes,'ygrid','on','xgrid','on','YTickMode','manual',...
    'Ytick',phase_ticks);
set(ph_axes,'FontWeight','bold');
