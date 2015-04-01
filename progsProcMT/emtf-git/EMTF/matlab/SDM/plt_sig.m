%***********************************************************************
%plt_sig : plots Signal or Noise poweer for array data
%USAGE  [hfig_sig] = plt_sig(S,ismooth,chid,stn,stid,c_title);
%
%INPUTS:
%    S(nt+1,nbt) is signal (or noise) power array
%       where : nt is total # of components in array
%               nbt is # of frequency bands
%    ismooth > 0 ==> interopolate curves with splines
%                    before plotting
%    chid  is list of channel IDs (i.e. Hx, Ey etc)
%    stn  is list identifyning station # for each channel 
%    stid  is list of station IDs  (3 character strings)
%          used for labeling plot
%    ctitle : title for plots
%  OUTPUTS:
%    hfig_sig is graphics handle for figure

function [hfig_sig] = plt_snr(S,ismooth,chid,stn,stid,c_title)

%  set figure window location and size (rect_snr) in stmonitr.m
stmonitr

scale = scale_snr;
rect = rect_snr;
rect = rect*scale;
rect_paper = [1.,1.,6.5,9.];
hfig_sig = figure('Position',rect,'PaperPosition',rect_paper);
ntemp = size(S);
nt1 = ntemp(1);
nt = nt1 - 1;
nbt = ntemp(2);
Stemp = log10(S(2:nt1,:)) ;
if(ismooth > 0)
  Ssm = smooth(Stemp,3);
  Ssm = 10.^Ssm;
  S1sm = smooth(S(1,:),3);
  S1sm = 10.^S1sm;
end
xmin = S(1,1)/1.1;
xmax = S(1,nbt)*1.1;
xt1 = ceil(log10(xmin)); xt2 = floor(log10(xmax)); 
xt2 = max(xt1,xt2); 
xt = 10.^[xt1:1:xt2]; 

lymin = floor(min(min(Stemp)));
lymax = ceil(max(max(Stemp)));
ymin = 10^lymin;ymax = 10^lymax;
Stemp = 10.^Stemp;

title_pos_y = lymin + .9*(lymax-lymin);
title_pos_y = 10^title_pos_y;
title_pos_x = log10(xmin) + .07*(log10(xmax/xmin));
title_pos_x = 10^title_pos_x;
line_styles = ['r- ';'b--';'g- ';'m- ';'c--';'kx ';'y+ '];
comp_plot = ['EX';'EY';'HX';'HY';'HZ'];

lims = [ xmin,xmax,ymin,ymax];

for iplot = 1:5;
  subplot(3,2,iplot);
  nlines = 0;
  for icomp = 1:nt;
    sta = stn(icomp);
    if( comp_plot(iplot,:)' == upper(setstr(chid(1:2,icomp))) )
      nlines = nlines + 1;
      if(ismooth > 0)
        loglog(S1sm,Ssm(icomp,:),line_styles(sta,:));
      else
        loglog(S(1,:),Stemp(icomp,:),line_styles(sta,:));
      end
      hold on
    end
  end
  fatlines(gca,line_thick)
  axis(lims)
  set(gca,'FontWeight','bold','FontSize',11)
  set(gca,'Xtick',xt);
  hold off
  text(title_pos_x,title_pos_y,comp_plot(iplot,1:2),'FontWeight','bold')
  if(iplot > 3) 
    xlabel('PERIOD (s)')
  end
  if(iplot < 3 )
    title(c_title)
  end
  if( (iplot == 1 ) | (iplot == 3 )| (iplot == 5 )  )
    ylabel('Power (nT^2 /Hz or (mV/km)^2 / Hz')
  end
end
%h_menu = uimenu('Label','Print','Parent',hfig_snr)
%         uimenu(h_menu,'Label','Printer','Callback','print -fgcf');
%         uimenu(h_menu,'Label','File','Callback',['print -fgcf ' uigetfile]);

%  lable box ... identify line style/colors with stations
nn = size(stid);
nsta = nn(2);
x = [0,1];
y = [1:nsta]'*ones(1,2);
axes('position',[.6,.1,.3,.25])
for k=1:nsta
   plot(x,y(k,:),line_styles(k,:));
   hold on
end
fatlines(gca,line_thick)
set(gca,'Visible','off','Xlim',[0,4],...
    'Ylim',[0,nsta+2])
for k=1:nsta
   text(2,k,setstr(stid(:,k)'),'FontWeight','bold','FontSize',14);
end
