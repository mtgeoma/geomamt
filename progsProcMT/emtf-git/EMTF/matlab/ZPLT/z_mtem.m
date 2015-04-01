%**********************************************************************************
%
%  Example of how to load in results from one Z-file, convert all quantities
%   (including error covariance matrices) into a common (usually geographic) coordinate
%    system, rotate if desired, compute rho and phi for off-diagonal components, and plot
%    the two modes. 

globinc
pltind = [];

[cfile,dir] = uigetfile(zid);
eval(['cd ' dir]);
l = length(cfile);
c_title = [replace_(cfile(1:l))];
cfile = [dir cfile ];
MeasurementCoordinates = 1;


%   read in transfer functions, residual covariance, inverse signal power
%   from Z.****** file
[Z,Sig_s,Sig_e,periods,ndf,stdec,orient,Nch,Nche,nbt,chid,csta] = Z_in(cfile);
eval(MKPLTIND);

NBT = [nbt];

% get info about dipole setup, channel ordering
[xy,yx,hz,xypairs,DipoleSetup] = ecomp(Nche,chid);

% set rotation angle 
th_x = orient(1,xy(1)+2);
theta = th_x;

if DipoleSetup == 'MT   '
%   just two E channels ... plot on one page

   [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,xypairs,theta);
                    
   %  xform  takes TFs input from Z file, converts to fixed coordinate system
   %    (x points in direction 0 used as reference for channel orientations)
   %    rotates, and computes rh, phi + standard errors
   %   theta defines rotation angle, ixy defines pairs of channels to rotate

   rho = [ rxy   ryx ];
   ph = [ pxy   pyx ];
   rho_se = 2*[ rxy_se   ryx_se ];
   ph_se = 2*[ pxy_se   pyx_se ];

   % reset plotting limits now that rho is known
   [lims] = set_lims(dir,periods,rho);

   %  set up figure
   [hfig] = set_fig(lims);
   uimenu('Parent',hfig,'Label','Plot Options',...
         'Callback','mt_menu')

   %location and size of plotting window on screen
   pltall = ones(length(periods),1);
   [rho_axes,ph_axes] = pltrhom(NBT,pltall,periods,rho,rho_se,ph,ph_se,lims,...
           c_title,hfig);
 elseif DipoleSetup == 'TEMAP'
    %   "tensor EMAP" ... multiple E channels with some cross lines
    %  here just plot in measurement coordinate system initially
 
   %   tranfrom to rho, phi, error bars (for all impedance elements
   %   in measurement coordinate system ... including diagonals)
   [rho,rho_se,ph,ph_se] = ap_res(Z,Sig_s,Sig_e,periods) ;
   %  WANT 2 se for plots ...
   rho_se = 2*rho_se;
   ph_se = 2*ph_se;

%   select out appropriate components for two modes
   rhoxy = get_mode(rho,2,nbt,xy);
   rhoxy_se = get_mode(rho_se,2,nbt,xy);
   phixy = get_mode(ph,2,nbt,xy);
   phixy_se = get_mode(ph_se,2,nbt,xy);
   rhoyx = get_mode(rho,1,nbt,yx);
   rhoyx_se = get_mode(rho_se,1,nbt,yx);
   phiyx = get_mode(ph,1,nbt,yx);
   phiyx_se = get_mode(ph_se,1,nbt,yx);
   
   % reset plotting limits now that rho is known
   [lims_xy] = set_lims(dir,periods,rhoxy);
 
   % reset plotting limits now that rho is known
   [lims_yx] = set_lims(dir,periods,rhoyx);
   
   lims(3) = min([lims_xy(3),lims_yx(3)]);
   lims(4) = max([lims_xy(4),lims_yx(4)]);
   lims(1) = lims_xy(1);
   lims(2) = lims_xy(2);
  
   %  set up figure
   [hfig_xy] = set_fig(lims,1);
   c_title_xy = [c_title ' :: XY Mode'];
   uimenu('Parent',hfig_xy,'Label','Plot Options',...
         'Callback','mt_menu')
 
   %location and size of plotting window on screen
   pltall = ones(length(periods),1);
   [rhoxy_axes,phixy_axes] = pltrhom(NBT,pltall,periods,rhoxy,rhoxy_se,phixy,...
           phixy_se,lims_xy,c_title_xy,hfig_xy);
   
  %  set up figure
   [hfig_yx] = set_fig(lims,2);
   c_title_yx = [c_title ' :: YX Mode'];
   
   %location and size of plotting window on screen
   [rhoyx_axes,phiyx_axes] = pltrhom(NBT,pltall,periods,rhoyx,rhoyx_se,phiyx,...
           phiyx_se,lims_yx,c_title_yx,hfig_yx);
           
else
   fprintf(1,'Not coded for this case')
end
