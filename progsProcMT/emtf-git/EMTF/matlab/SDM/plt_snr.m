%***********************************************************************
%plt_snr : plots Signal-to-Noise ratios for array data
%USAGE  [hfig_snr] = plt_snr(S,NI,ismooth,chid,stn,stid,c_title);
%
%INPUTS:
%    S(nt+1,nbt) is signal power array
%    NI(nt+1,nbt) is noise power array
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
%    hfig_snr is graphics handle for figure

function [hfig_snr] = plt_snr(S,N,ismooth,chid,stn,stid,c_title)

%  set figure window location and size (rect_snr) in stmonitr.m
%scale = .9
%rect = [550,100,550,800];
stmonitr

scale = scale_snr;
rect = rect_snr;
rect = rect*scale;
rect_paper = [1.,1.,6.5,9.];
hfig_snr = figure('Position',rect,'PaperPosition',rect_paper);
ntemp = size(S);
nt1 = ntemp(1);
nt = nt1 - 1;
nbt = ntemp(2);
SNR = 10*(log10(S(2:nt1,:)./N(2:nt1,:))) ;
if(ismooth > 0)
  SNRsm = smooth(SNR,3);
  S1sm = smooth(S(1,:),3);
end
xmin = S(1,1)/1.1;
xmax = S(1,nbt)*1.1;
ymin = -10;
ymax = 10*(fix(max(max(SNR))/10)+1);
title_pos_y = ymin + .9*(ymax-ymin);
title_pos_x = log10(xmin) + .07*(log10(xmax/xmin));
title_pos_x = 10^title_pos_x;
line_styles = ['r- ';'b--';'g- ';'m- ';'c--';'kx ';'y+ '];
comp_plot = ['EX';'EY';'HX';'HY';'HZ'];

xt1 = ceil(log10(xmin)); xt2 = floor(log10(xmax));
xt2 = max(xt1,xt2);
xt = 10.^[xt1:1:xt2];

lims = [ xmin,xmax,ymin,ymax];

for iplot = 1:5;
  subplot(3,2,iplot);
  nlines = 0;
  for icomp = 1:nt;
    sta = stn(icomp);
    if( comp_plot(iplot,:)' == upper(setstr(chid(1:2,icomp))) )
      nlines = nlines + 1;
      if(ismooth > 0)
        semilogx(S1sm,SNRsm(icomp,:),line_styles(sta,:));
      else
        semilogx(S(1,:),SNR(icomp,:),line_styles(sta,:));
      end
      hold on
    end
  end
  fatlines(gca,line_thick)
  set(gca,'Xlim',lims(1:2),'Ylim',lims(3:4))
  set(gca,'FontWeight','bold','FontSize',11)
  set(gca,'XTick',xt)
  hold off
  text(title_pos_x,title_pos_y,comp_plot(iplot,1:2),'FontWeight','bold')
  if(iplot > 3) 
    xlabel('PERIOD (s)')
  end
  if(iplot < 3 )
    title(c_title)
  end
  if( (iplot == 1 ) | (iplot == 3 )| (iplot == 5 )  )
    ylabel('SNR : dB')
  end
end
h_menu = uimenu('Label','Print','Parent',hfig_snr)
         uimenu(h_menu,'Label','Printer','Callback','print -fgcf');
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
