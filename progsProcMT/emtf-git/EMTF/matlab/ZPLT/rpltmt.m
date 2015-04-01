%mtrplt  replots MT figure
   if 1-min(lims_old == lims)
      %  set up figure with new limits
      delete(hfig);
      [hfig] = set_fig(lims);
      uimenu('Parent',hfig,'Label','Plot Options',...
         'Callback','mt_menu')
   else
      delete(rho_axes);
      delete(ph_axes);
   end   

   if theta ~= theta_old
      % transform/rotate
      [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,xypairs,theta);
      rho = [ rxy   ryx ];
      ph = [ pxy   pyx ];
      rho_se = 2*[ rxy_se   ryx_se ];
      ph_se = 2*[ pxy_se   pyx_se ];
   end
         
%   replot in any case
 
   %location and size of plotting window on screen
   [rho_axes,ph_axes] = pltrhom(NBT,pltind,periods,rho,rho_se,ph,ph_se,lims,...
           c_title,hfig);
