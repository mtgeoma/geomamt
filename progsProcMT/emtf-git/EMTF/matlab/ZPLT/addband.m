[cfile,dir] = uigetfile(zid);
eval(['cd ' dir]);
l = length(cfile);
c_title = ['Site: ' cfile(3:l)];
cfile = [dir cfile ];

%    addband   adds results from an additional sampling band to existing plot

%   read in transfer functions, residudual covariance, inverse signal power
%   from Z.****** file
[Z1,Sig_s1,Sig_e1,periods1,ndf,stdec,orient,Nch1,Nche,nbt,chid1,csta] = Z_in(cfile);
if Nch == Nch1 
  if min(min((upper(chid1) == upper(chid) )))
%    number of channels, IDs, agree with previous ... merge new Z_*** plot with old  
    NBT = [NBT nbt];
    eval(MKPLTIND);
    Z = [Z Z1];
    Sig_s = [ Sig_s Sig_s1 ];
    Sig_e = [ Sig_e Sig_e1 ];
    periods = [ periods periods1 ];
    NB = sum(NBT);
 
    if(DipoleSetup == 'MT   ')
      delete(hfig);
    
     % transform/rotate
      [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,xypairs,theta);
      rho = [ rxy   ryx ];
      ph = [ pxy   pyx ];
      rho_se = 2*[ rxy_se   ryx_se ];
      ph_se = 2*[ pxy_se   pyx_se ]; 
 
      ind = find(pltind);
      [lims] = set_lims(dir,periods(ind),rho(ind,:));
      [hfig] = set_fig(lims);
      uimenu('Parent',hfig,'Label','Plot Options',...
         'Callback','mt_menu')
      %location and size of plotting window on screen
      [rho_axes,ph_axes] = pltrhom(NBT,pltind,periods,rho,rho_se,ph,ph_se,lims,...
           c_title,hfig);   
    else
      if MeasurementCoordinates

        %   tranfrom to rho, phi, error bars (for all impedance elements
        %   in measurement coordinate system ... including diagonals)
        [rho,rho_se,ph,ph_se] = ap_res(Z,Sig_s,Sig_e,periods) ;
        %  WANT 2 se for plots ...
        rho_se = 2*rho_se;
        ph_se = 2*ph_se;

        %   select out appropriate components for two modes
        rhoxy = get_mode(rho,2,NB,xy);
        rhoxy_se = get_mode(rho_se,2,NB,xy);
        phixy = get_mode(ph,2,NB,xy);
        phixy_se = get_mode(ph_se,2,NB,xy);
        rhoyx = get_mode(rho,1,NB,yx);
        rhoyx_se = get_mode(rho_se,1,NB,yx);
        phiyx = get_mode(ph,1,NB,yx);
        phiyx_se = get_mode(ph_se,1,NB,yx);
        ind = find(pltind);
        [lims_xy] = set_lims(dir,periods(ind),rhoxy(ind)); 
        [lims_yx] = set_lims(dir,periods(ind),rhoyx(ind)); 
        lims(3) = min([lims_xy(3),lims_yx(3)]);
        lims(4) = max([lims_xy(4),lims_yx(4)]);
        lims(1) = lims_xy(1);
        lims(2) = lims_xy(2); 
      else  
        % transform/rotate
        [rhoyx,rhoxy,phiyx,phixy,rhoyx_se,rhoxy_se,phiyx_se,phixy_se] = ...
                     xform(Z,Sig_s,Sig_e,periods,orient,xypairs,theta);
 
        ind = find(pltind);
        [lims_xy] = set_lims(dir,periods(ind),rhoxy(ind)); 
        [lims_yx] = set_lims(dir,periods(ind),rhoyx(ind)); 
        lims(3) = min([lims_xy(3),lims_yx(3)]);
        lims(4) = max([lims_xy(4),lims_yx(4)]);
        lims(1) = lims_xy(1);
        lims(2) = lims_xy(2); 
      end
      rhoxy_se = 2*rhoxy_se;
      phixy_se = 2*phixy_se;
      rhoyx_se = 2*rhoyx_se;
      phiyx_se = 2*phiyx_se;
      % set plotting limits 
      delete(hfig_xy);
      delete(hfig_yx); 
      [hfig_xy] = set_fig(lims,1);
      uimenu('Parent',hfig_xy,'Label','Plot Options',...
         'Callback','mt_menu')
      c_title_xy = [c_title ' :: XY Mode'];
 
      %location and size of plotting window on screen
      [rhoxy_axes,phixy_axes] = pltrhom(NBT,pltind,periods,rhoxy,rhoxy_se,phixy,...
           phixy_se,lims,c_title_xy,hfig_xy);
           
      [hfig_yx] = set_fig(lims,2);
      c_title_yx = [c_title ' :: YX Mode']
 
      %location and size of plotting window on screen
      [rhoyx_axes,phiyx_axes] = pltrhom(NBT,pltind,periods,rhoyx,rhoyx_se,phiyx,...
           phiyx_se,lims,c_title_yx,hfig_yx);         
    end  
  else
    fprintf(1,'Channel ID''s for new file don''t agree with existing plot\n');
  end
else
  fprintf(1,'Number of channels in new file don''t agree with existing plot\n');
end
