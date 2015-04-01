%  This script reads in reuslts output in the array TF
%  file Pw**** output by multmtrn, and plots any chosen interstation and/or
%  intercomponent transfer functions desired, in any chosen coordinate
%  system

%  GLOBAL VARIABLES
%  channels chosen for plotting, rotation are global variables
%  ref_ch  =  reference channels (currently has to be an actual channel #)
%  rot_ch  =  channel numbers for rotation pairs
%  sing_ch =  channel numbers for single (unpaired for rotation) channels
%  plot_ch =  channels to plot TFs for
%  theta_rot = rotation angle
%  poluse = polarization to plot (1 = Hx , 2 = Hy)
%  chplot = channel number in rotated channel list
global ref_ch sing_ch rot_ch plot_ch theta_rot poluse chplot

%  all of these are associated with dialogue box which controls
%  channel rotations etc.
global chnum chuse chname chhand name_list nlines plt_type

%  default plotting types
ptyps = ['AMPH';'LAMP';'REIM'];
plt_type = ptyps(1,:);

%   get pathname for Pw file to plot
[cfile,cpath] = uigetfile('*.Pw','Choose Pw-File To Plot');
cfile = [cpath cfile]; 

%  load in Pw file and header into Pw and Pwhd data structures
[pw,pwhd] = pwstruct(cfile);

% open dialog box for choosing rotation pairs, polarizations,
[h_pwdialog] = pwdialog(pwhd);

h_plots = [];
