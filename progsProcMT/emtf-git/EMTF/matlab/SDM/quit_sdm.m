global hfig_can_coh_menu
delete(hfig_eval); clear hfig_eval;

 
if(exist('hfig_snr'))
   chk_clr(hfig_snr);
   clear hfig_snr;
end
if(exist('hfig_noise'))
   chk_clr(hfig_noise);
   clear hfig_noise;
end
if(exist('hfig_sig'))
   chk_clr(hfig_sig);
   clear hfig_sig;
end
if(exist('hfig_cc')) 
  chk_clr(hfig_cc);
  clear hfig_cc;
end
if(exist('hfig_evec_menu')) 
  chk_clr(hfig_evec_menu);
  clear hfig_evec_menu;
end

if(exist('hfig_can_coh_menu') & ~ isempty(hfig_can_coh_menu))
  chk_clr(hfig_can_coh_menu);
  clear hfig_can_coh_menu;
end
for h=hfig_evec
  chk_clr(h); 
end
