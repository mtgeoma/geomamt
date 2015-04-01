function h_check = ch_menu(chid,csta,loc);
global hfig_can_coh_menu
[nchid,nt] = size(chid);
stmonitr
width = width_cc;
height_1 = height_cc_1;
height = (nt+3)*height_1;
if(nargin < 3)
  loc = loc_cc;
end
rect_fig = [loc width height ];
hfig_can_coh_menu = figure('Position',rect_fig);
for k=1:nt
  rect = [ width*.1 (nt-k)*height_1+5 width*.8 height_1*.8];
  name = [csta(:,k)' '  ' chid(:,k)' ];
  h_check(k) = uicontrol('style','checkbox', ...
               'parent',hfig_can_coh_menu,'position',  ...
         rect,'string',name,'callback','find_ch(h_check);');
end
rect = [ width*.1 (nt+1.6)*height_1+5 width*.8 height_1*.8];
h_plot = uicontrol('style','pushbutton', ...
           'parent',hfig_can_coh_menu,'position',  ...
       rect,'string','PLOT','callback','sdm_sub');
rect = [ width*.1 (nt+.5)*height_1+5 width*.8 height_1*.8];
h_quit = uicontrol('style','pushbutton','parent',hfig_can_coh_menu,'position', ...
    rect,'string','QUIT','callback','cc_quit');
set(gcf,'Name','Sta./Comp.');
return
