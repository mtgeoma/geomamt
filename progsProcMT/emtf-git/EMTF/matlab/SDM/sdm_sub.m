%*********************************************************************
%  sdm_sub  script to load in evecs, evals, reconstruct singular sdm
%           and then extract submatrices for 2 sets of channels
%           and computes and plot eigenvalues for each diagonal submarix
%           plus singular values for cross-product matrix (i.e.,
%           canonical covariances;
%  NOTE:  need to set ind1, ind2, and initialize input file with
%          call to sdm_init

line_sty1 = 'r-';
line_sty2 = 'b--';

n1 = length(ind1);
n2 = length(ind2);
TITLE1 = ['GROUP 1 =  '];
for k=1:n1-1
   TITLE1 = [ TITLE1 csta(:,ind1(k))' ' ' chid(:,ind1(k))' ' : ' ];
end
TITLE1 = [ TITLE1 csta(:,ind1(n1))' ' ' chid(:,ind1(n1))' ];
TITLE2 = [ 'GROUP 2 =  '];
for k=1:n2-1
   TITLE2 = [ TITLE2 csta(:,ind2(k))' ' ' chid(:,ind2(k))' ' : ' ];
end
TITLE2 = [ TITLE2 csta(:,ind2(n2))' ' ' chid(:,ind2(n2))' ];

nc = min(n1,n2);
ev1 = zeros(n1,nbt);
ev2 = zeros(n2,nbt);
ccov = zeros(nc,nbt);
ccor = zeros(nc,nbt);
for ib = 1:nbt
   [period(ib),nf(ib),var,S] = sdm_in(fid_sdm,nt,ib,irecl);
   sig = 1./sqrt(var);
   S11 = S(ind1,ind1);S22 = S(ind2,ind2); S12 = S(ind1,ind2);
   temp = real(eig((S12'/S11)*S12,S22));
   temp = sort(temp);
   ccor(:,ib) = temp(n2:-1:n2-nc+1);
   S11 = diag(sig(ind1))*S11*diag(sig(ind1));
   S22 = diag(sig(ind2))*S22*diag(sig(ind2));
   S12 = diag(sig(ind1))*S12*diag(sig(ind2));
   temp = abs(real(eig(S11)));
   temp = sort(temp);
   ev1(:,ib) = temp(n1:-1:1);
   temp = abs(real(eig(S22)));
   temp = sort(temp);
   ev2(:,ib) = temp(n2:-1:1);
   temp = svd(S12);
   temp = -sort(-temp);
   ccov(:,ib) = temp(1:nc);
end


%  now plot

stmonitr
rect_win = rect_cc;

%rect_win =  [50,50,700,700];
rect_paper = [1.,1.,6.5,6.5];
rect_11 = [.1,.52,.35,.35];
rect_22 = [.55,.1,.35,.35];
rect_21 = [.55,.52,.35,.35];
rect_12 = [.1,.1,.35,.35];
rect_tit = [ .2 .93 .6 .05 ];
if( exist('hfig_cc') )
   delete(hfig_cc);
end
hfig_cc = figure('Position',rect_win,'PaperPosition',rect_paper);
xmin = period(1)/1.1 ;
xmax = period(nbt)*1.1 ;
ev1 = 10*log10(ev1);
ev2 = 10*log10(ev2);
emax = max(max(ev1));
emax2= max(max(ev2));
emax = max([emax emax2]);
ccov = 10*log10(ccov);
ymin = -10
ymax = 10*(ceil(max(emax)/10)+1);
ymax = max(ymax,10)
x_lab = (log10(xmax)-log10(xmin))*.1+log10(xmin);
x_lab = 10^x_lab;
y_lab = .85*ymax;

h11 = axes('Position',[rect_11]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ev1,3);
   for c = 1:n1 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:n1 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,ev1(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end
tit_all = 'Eigenvalues of S11';
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)

lims = [ xmin,xmax,ymin,ymax];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13)
ylabel('SNR : dB');
text(-.1,1.25,TITLE1,'Units','Normalized','FontSize',10,...
            'FontWeight','bold','Color','red');
text(-.1,1.15,TITLE2,'Units','Normalized','FontSize',10,...
             'FontWeight','bold','Color','blue');

h22 = axes('Position',[rect_22]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ev2,3);
   for c = 1:n2 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:n2 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,ev2(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

tit_all = 'Eigenvalues of S22';
text(x_lab,y_lab,tit_all,'FontSize',11,'FontWeight','bold');
lims = [ xmin,xmax,ymin,ymax];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13)
xlabel('PERIOD (s)');
ylabel('SNR : dB');

h12 = axes('Position',[rect_21]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ccov,3);
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1; 
      else
          linestyle = line_sty2;
      end
      semilogx(period,ccov(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

%   last canonical correlations
tit_all = 'Canonical Covariances' ;
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)
lims = [ xmin,xmax,ymin,ymax];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13)
ylabel('SNR : dB');

h12 = axes('Position',[rect_12]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ccor,3);
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,ccor(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

y_lab = .85;
tit_all = 'Canonical Coherence';
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)
lims = [ xmin,xmax,0,1];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13)
xlabel('PERIOD (s)');
ylabel('Coherence ');
%h_menu = uimenu('Label','Print','Parent',hfig_cc) 
%         uimenu(h_menu,'Label','Printer','Callback','print -fgcf'); 
%         uimenu(h_menu,'Label','File','Callback',['print -fgcf ' uigetfile]);
