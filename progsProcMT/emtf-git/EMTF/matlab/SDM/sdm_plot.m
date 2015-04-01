 %  matlab script (not function) for looking at S0**** files
%  Run this script from the directory containing the sdm files (S0***)
%   that you want to look at
%  After choosing a file from the menu, a plot of sdm eigenvalues
%  vs. period will appear.  There are two additional sorts of plots
%  so far: (1) there is a slider bar beneath the period axis on the
%   eigenvector plot ... with the slider, and the "PLOT" pushbutton
%   you can plot the eigenvectors.  These appear as a plot of
%  magnetic and electric vectors on a map of station locations 
%  (note that this assumes station locations are in the sdm files!)
%  (2) you can plot eigenvalues and canonical coherences for subgroups
%  of data channels.  You get to this via the pull down menu called
%   "Plot Options".  A set of pushbuttons listing all components in 
%  the array pops up ... you use this to pick the components in "group 1";
%  the remainder are in group 2.  The program then plots eigenvalues
%  for each group separately + canonical coherences and covariance.
%   If there is coherent noise present which only occurs at some sites,
%   or which only occurs in E components (e.g.), then this should
%   help you figure out which sites/components are not contaminated.

clear ll_lim fid_sdm irecl nbt nt nsta nch stcor decl sta chid ...
            orient periods asp Hp Ep Hz hfig_eval

global ll_lim fid_sdm irecl nbt nt nsta nch ...
                 stcor decl sta chid orient periods asp
global Hp Ep Hz
global hfig_eval

hfig_evec = [];
clear hfig_cc;
clear hfig_snr;
clear h_ib;
[cfile,cpath] = uigetfile('*.S0');
eval(['cd ' cpath]);
arrayid = cfile(1:length(cfile)-3);
sn_file  = [ cpath arrayid '.SN'];
s0_file = [ cpath cfile ];
[ndf,E,Sig,Noise] = sn_in(sn_file);
[fid_sdm,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,decl,sta,chid,csta, ...
                                     orient,periods] = sdm_init(s0_file);

stn = [];
for ista = 1:nsta
   stn = [ stn ista*ones(1,nch(ista)) ];
end

%  rho_ref is the apparent resistivity assumed for scaling E-field vectors
%  into nT ... if actual apparent resistivity = rho_ref, then H and E vectors
%  will have the same length on eigenvector plots 
%  ivec is an array which tells which eigenvectors of the SDM should
%  be plotted for each band; 1 = eigenvalue associated with largest eigenvector, etc.
%           ( both of these should ultimately be changeable from menu ...)
rho_ref = 100;
n_evec = 3; ivec = [1:n_evec];
snr_units = 0; l_ellipse = 1;

title_sig = [arrayid ': Signal Power'];
title_noise = [arrayid ': Noise Power'];

ib = 1;
ismooth = 1;
hfig_eval = plt_eval(E,ndf,ismooth,arrayid);
h_menu = uimenu('Parent',hfig_eval,'Label','Signal/Noise');
         uimenu(h_menu,'Label','Signal Power','Callback', ...
          'hfig_sig = plt_sig(Sig,0,chid,stn,sta,title_sig)');
         uimenu(h_menu,'Label','Noise Power','Callback', ...
          'hfig_noise = plt_sig(Noise,0,chid,stn,sta,title_noise)');
         uimenu(h_menu,'Label','SNR','Callback', ...
           'hfig_snr = plt_snr(Sig,Noise,1,chid,stn,sta,arrayid);');
         uimenu(h_menu,'Label','Can. Coh','Callback','can_coh');
         uimenu(h_menu,'Label','Print','Callback','print');
         uimenu(h_menu,'Label','Quit','Callback','quit_sdm');
h_menu_ev_opt = uimenu('Parent',hfig_eval,'Label','Eigenvectors',...
         'Callback','evecmenu');
slide_min = log10(periods(1)); 
slide_max = log10(periods(nbt));
slide_val = ceil(slide_min);
h_slider = uicontrol('Parent',hfig_eval,'Style','slider',...
          'Units','normalized','Position',[.13,.30,.84,.03], ...   
          'Max',slide_max,'Min',slide_min,'Value',slide_val, ...
          'CallBack','get_ib');
h_evec_plot = uicontrol('Parent',hfig_eval,'Style','pushbutton', ...
          'Units','normalized','Position',[.025,.29,.10,.04], ...
         'string','PLOT','CallBack', ...
          'hev1=evec_plt(ib,ivec,rho_ref,snr_units,l_ellipse);hfig_evec=[hfig_evec hev1];');
[fid_sdm,irecl,nbt,nt,nsta,nsig,nch,ih,stcor ...
         ,decl,sta,chid,csta,orient,periods] = sdm_init(s0_file);
%evec_set;
%  sets up for eigenvector plotting ... call once to open file/intialize

%  testing 1,2,3
%    stcor = stcor+ones(2,1)*[0:3];
%    chid = ['Hx';'Hy';'Hz';'Ex';'Ey'; ...
%        'Hx';'Hy';'Hz';'Ex';'Ey'; ...
%        'Hx';'Hy';'Hz';'Ex';'Ey'; ...
%        'Hx';'Hy';'Hz';'Ex';'Ey'];
%   chid = chid';

[Hp,Ep,Hz] = ch_pair(nch,chid);

% some calculations for scaling plot sizes ...
lat_max =max(max(stcor(1,:)));
lat_min =min(min(stcor(1,:)));
lon_max =max(max(stcor(2,:)));
lon_min =min(min(stcor(2,:)));
lat_av = (lat_max+lat_min)/2.;
lat_range = (lat_max-lat_min);
lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);

marg = 1.5/(sqrt(nsta))*max([lat_range,lon_range]);
if(marg == 0)
% %&*&^(*(7%$%$#  user forgot to define site locations in sp file
%              ...  yet again !!!
%   so define dummy site locations ...
   nn = size(stcor);nsta = nn(2);
   stcor = [ 1:nsta; 1:nsta ];
   lat_max =max(max(stcor(1,:)));
   lat_min =min(min(stcor(1,:)));
   lon_max =max(max(stcor(2,:)));
   lon_min =min(min(stcor(2,:)));
   lat_av = (lat_max+lat_min)/2.;
   lat_range = (lat_max-lat_min);
   lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);
   marg = .6*max([lat_range,lon_range]);
end
if(marg == 0) marg ==1 ; end

asp = (2*marg+lon_range)/(2*marg+lat_range);
lat_max = lat_max+marg;
lon_max = lon_max+marg;
lat_min = lat_min-marg;
lon_min = lon_min-marg;
ll_lim = [ lat_min,lat_max,lon_min,lon_max];
