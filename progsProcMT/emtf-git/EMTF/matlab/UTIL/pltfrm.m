% include file to set values of paramters appropriate to
%  running on specific computers (mostly all positions/sizes
%  of figure windows, to accomodate different computer screens

%  DELL NOTEBOOK

%  eval_plt (main window)
size_fac2 = 50;
size_fac1 = 60;
rect_win =  [20,40,400,520];

%  plt_snr  (Signal-to-noise plot)
scale_snr = 1;
rect_snr = [300,40,500,520];

%  canonical coherences
rect_cc =  [50,50,500,500];
width_cc = 125;
height_cc_1 = 25;
loc_cc  = [600,50];

%  evec_plt  (eigenvector plots)
lower_left = [50,50];
pix_x = 180;
extra = 100;
space = 25;


%  SUN workstation
%rect_win =  [20,200,450,600];
%scale_snr = .9
%rect_snr = [550,100,550,800];

%  evec_plt  (eigenvector plots)
%lower_left = [50,50];
%pix_x = 200;
%extra = 100;
%space = 25;

%  canonical coherences
% rect_cc =  [50,50,700,700];

