%   plots some eigenvectors for a single specified frequency band
%
%  Usage: 

function [hfig_evec] = evec_plt(ib,ivec,rho_ref,snr_units,l_ellipse)

global ll_lim fid_sdm irecl nbt nt nsta nch ...
                 stcor decl sta chid orient periods asp
global Hp Ep Hz

%  get this figure size info out of stmonitr.m now
%lower_left = [50,50];
%pix_x = 200;
%extra = 100;
%space = 25;
%sfac = .7;
stmonitr

%  read in sdm for band ib
[period,nf,var,S] = sdm_in(fid_sdm,nt,ib,irecl);
%  solve generalized eigenvalue problem
[u,eval] = eig(S,diag(var));
if( 1-snr_units) 
   u = diag(var)*u;
else
   u = sqrt(diag(var))*u;
end
%  make sure eigenvalues are in correct order ...
[temp,ind] = sort(diag(eval));
nvec = length(ivec);
%  and select out the desired vectors
u(:,1:nvec) = u(:,ind(nt-ivec+1));
nn = size(Hp);
Hind = reshape(Hp(:,1:2),nn(1)*2,1);

% calculate window sizes, figure scalings
pix_y = pix_x/asp;
width = (pix_x+space)*nvec+extra;
height = pix_y + extra;
norm_width = pix_x/width;
x_step = (pix_x+space)/width;
x0 = extra/width;
y0 = .5*extra/height;
norm_height = pix_y/height;
rect_fig = [ lower_left [width height ]]; 
if(snr_units) 
  figname = [ 'Band =  ' num2str(ib) ...
             '     ::   Period =  ' num2str(period)  ' sec.      '...
                          'SNR Units']
else
  figname = [ 'Band = ' num2str(ib) ...
             '     ::   Period =  ' num2str(period)  ' sec.      '...
             ' ::   Ref rho =  ' num2str(rho_ref) ];
end
hfig_evec=figure('Position',rect_fig,'Name',figname,'NumberTitle','off');
rect_plt = [ x0 , y0 , norm_width, norm_height ];

% loop over desired eigenvectors
axlab = [1,1]; 
for k=1:nvec
   ctit = ['Eigenvector #', num2str(ivec(k))];
   u(:,k) = chngph(u(:,k),Hind);
   [uH,uE,uZ] = u_pair(u(:,k),Hp,Ep,Hz,orient,decl,stcor,period,...
                            rho_ref,snr_units);
   hfig = evplt_HE(rect_plt,ll_lim,sfac,uH,uE,axlab,ctit,l_ellipse);
   if(k == 1 ) 
       xtxt = ll_lim(1) + .2*(ll_lim(2)-ll_lim(1));
       ytxt = ll_lim(3) - .15*(ll_lim(4)-ll_lim(3));
       text('Position',[xtxt,ytxt],'string','Green: Magnetics',...
           'Color',[0,.7,0],'FontSize',12,'FontWeight','bold');
   elseif (k == nvec)
       xtxt = ll_lim(1) + .2*(ll_lim(2)-ll_lim(1));
       ytxt = ll_lim(3) - .15*(ll_lim(4)-ll_lim(3));
       text('Position',[xtxt,ytxt],'string','Red: Electrics',...
            'Color','r','FontSize',12,'FontWeight','bold');
   end
   rect_plt(1)  = rect_plt(1) + x_step;
   axlab = [1,0];
end
return
end
