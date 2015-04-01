%rpltem  replots TEMAP figures
   if 1-min(lims_old == lims)
      %  set up figures with new limits
      delete(hfig_xy);
      [hfig_xy] = set_fig(lims,1);
      delete(hfig_yx);
      [hfig_yx] = set_fig(lims,2);
      
      uimenu('Parent',hfig_xy,'Label','Plot Options',...
         'Callback','mt_menu')
   else
      delete(rhoxy_axes);
      delete(phixy_axes);
      delete(rhoyx_axes);
      delete(phiyx_axes);
   end   

   if theta ~= theta_old
      % transform/rotate
      [rhoyx,rhoxy,phiyx,phixy,rhoyx_se,rhoxy_se,phiyx_se,phixy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,xypairs,theta);
      rhoxy_se = 2*rhoxy_se;
      phixy_se = 2*phixy_se;
      rhoyx_se = 2*rhoyx_se;
      phiyx_se = 2*phiyx_se;
   end
         
%   replot in any case

hfig_xy
  [rhoxy_axes,phixy_axes] = pltrhom(NBT,pltind,periods,rhoxy,rhoxy_se,phixy,...
           phixy_se,lims,c_title_xy,hfig_xy);

  [rhoyx_axes,phiyx_axes] = pltrhom(NBT,pltind,periods,rhoyx,rhoyx_se,phiyx,...
           phiyx_se,lims,c_title_yx,hfig_yx);
