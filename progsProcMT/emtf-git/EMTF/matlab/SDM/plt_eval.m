%***********************************************************************
%plt_eval :  plots normalized eigenvalues of SDM
%USAGE:  [hfig] = plt_eval(E,ndf,ismooth,c_title)
%      INPUTS:
%        E(nt+1,nbt) = eigenvalue vs. period array
%           E(1,:) contains period
%           nt, nbt are total # of components, # of bands
%        ndf(nbt) = # of data points used for each band
%           (if called with an empty string ndf, this is not plotted)
%        ismooth > 0 will cause curves to be interpolated with
%           a cubic spline before plotting
%        c_title = title
%      OUTPUTS:
%        hfig = graphics handle for figure
%
%**********************************************************************

function [hfig] = plt_eval(E,ndf,ismooth,c_title)

%  get location, size of eigenvalue plot window (i.e., rect_win)
%     from stmonitr.m
%rect_win =  [20,200,450,600];
stmonitr

rect_paper = [1.,1.,6.5,9.];
rect_eval = [.15,.4,.8,.5];
rect_ndf = [.15,.08,.8,.18];
hfig = figure('Position',rect_win,'PaperPosition',rect_paper);
nm = size(E);
nt1 = nm(1);
nbt = nm(2);
period = E(1,:);
xmin = E(1,1)/1.1 ;
xmax = E(1,nbt)*1.1 ;
E = 10*log10(E) ;
ymin = -10;
ymax = 10*(fix(max(E(2,:))/10)+1);
h_ax_eval = axes('Position',[rect_eval]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(E(2:nt1,:),3);
   for c = 1:nt1-1 ;
      if( c <= 2 )
          linestyle = 'r';
      else
          linestyle = 'b--';
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      if(c <= 2) 
        line_handles = get(gca,'Children');
        set(line_handles(1),'LineWidth',line_thick);
      end
      hold on ;
   end
else
   for c = 1:nt1-1 ;
      if( c <= 2 )
          linestyle = 'r';
      else
          linestyle = 'b--';
      end
      semilogx(period,E(c+1,:),linestyle);
      if(c <= 2) 
        line_handles = get(gca,'Children');
        set(line_handles(1),'LineWidth',line_thick);
      end
      hold on ;
   end
end

lims = [ xmin,xmax,ymin,ymax];
axis(lims);
tit_all = ['Eigenvalues of SDM : SNR Units  : ',c_title];
title(tit_all);
set(get(gca,'Title'),'FontWeight','bold','FontSize',12)
set(gca,'FontWeight','bold','FontSize',11)
xlabel('PERIOD (s)')
ylabel('SNR : dB');

nm = size(ndf);
if(nm(1) == nbt )
% also plot sample sizes
   nmax = log10(max(ndf));
   nmax = 10^(fix(nmax)+1); 
   nmin = 1;
   lims = [ xmin,xmax,nmin,nmax];
   h_ax_ndf = axes('Position',[rect_ndf]);
   loglog(period,ndf,'*');
   axis(lims);
   title_pos_y = log10(nmin) + .2*(log10(nmax/nmin));
   title_pos_x = log10(xmin) + .05*(log10(xmax/xmin));
   title_pos_y = 10^title_pos_y;
   title_pos_x = 10^title_pos_x;
   text(title_pos_x,title_pos_y,'# of Data Vectors','FontWeight','bold');
   set(gca,'FontWeight','bold','FontSize',11)
end
end
