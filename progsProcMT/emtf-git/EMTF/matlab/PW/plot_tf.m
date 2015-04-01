%  rotates Pw response, computes TFs relative to chosen reference,
%  add TFs to array v1, and plots

[pwRot] = rotatePw(pw,pwhd,rot_ch,sing_ch,theta_rot);
[v,sig_v] = pwrfsite(pwRot,[1 2]);

nt_rot = 2*length(rot_ch(:,1))+length(sing_ch);
XY = ['xy'];
l = poluse; 
labels = [ labels ; [chname(l,1,1) XY(l) squeeze(chname(l,1,3:6))']];
ifref = [ ifref 1 ];
kk = 0;
nrot = length(rot_ch(:,1));   
nsing = length(sing_ch);   
chplot = [];
for k = 1:nrot 
   for l = 1:2 
      kk = kk + 1;
      if(chuse(l,k))  
         chplot = [ chplot kk ]; 
         labels  = [ labels ; ...
         [ chname(l,k,1)  XY(l) squeeze(chname(l,k,3:6))' ]];
         ifref = [ ifref 0 ];
      end
   end
end
for k = 1:nsing
   kk = kk + 1;
   if(chuse(1,k+nrot))  
      chplot = [ chplot kk ];
      labels  = [labels; squeeze(chname(1,nrot+k,1:6))'];
      ifref = [ ifref 0 ];
   end
end
nplot = length(chplot);
labels = char(labels);

temp = squeeze(v(poluse,chplot,:));
if(nplot ~= 1) temp = temp.'; end
v1 = [ v1 temp ];
amp = abs(v1);
temp = squeeze(sig_v(poluse,chplot,:));
if(nplot ~= 1) temp=temp'; end
sig_v1 = [ temp sig_v1 ] ;

if(nplot == 0)
   fprintf(2,'%s','No channels selected for plotting');
else
   if plt_type == 'AMPH'
%  	plot selected channels as amplitudes and phases
		ph = (180/pi)*atan2(imag(v1),real(v1));
		tr = [];ti = [];
      lims = pwstlims(T,amp,ph,tr,ti,plt_type);
   	amp_se = 2*real(sig_v1)/sqrt(2);
   	ph_se = 2*(180/pi)*real(sig_v1)./(amp*sqrt(2));
   	c_title = [ 'Rotation Angle = ' num2str(fix(theta_rot)) ];
   	hfig_plt = plt_amp_ph(T,amp,amp_se,ph,ph_se,lims,plt_type,...
         c_title,fig_pos,labels,ifref);
      h_plots = [ h_plots hfig_plt];
	elseif plt_type == 'REIM'
%  	plot real and imaginary parts of selected channels 
   	ph = []; amp = [];
   	tr = real(v1) ;ti = imag(v1);
   	lims = pwstlims(T,amp,ph,tr,ti,plt_type);
      tr_se = 2*real(sig_v1)/sqrt(2);
      c_title = [ 'Rotation Angle = ' num2str(fix(theta_rot)) ];
      hfig_plt = plt_re_im(T,tr,ti,tr_se,lims,c_title,fig_pos,labels,ifref);
      h_plots = [ h_plots hfig_plt];
	else
  	 %  plot selected channels as amplitudes and phases with log-log scale
   	ph = (180/pi)*atan2(imag(v1),real(v1));
   	tr = [];ti = [];
   	lims = pwstlims(T,amp,ph,tr,ti,plt_type);
   	amp_se = 2*real(sig_v1)/sqrt(2);
   	ph_se = 2*(180/pi)*real(sig_v1)./(amp*sqrt(2));
   	c_title = ''
   	hfig_plt = plt_amp_ph(T,amp,amp_se,ph,ph_se,lims,plt_type,...
         c_title,fig_pos,labels,ifref);
      h_plots = [ h_plots hfig_plt];
   end
end
