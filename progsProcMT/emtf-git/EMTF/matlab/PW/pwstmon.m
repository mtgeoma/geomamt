% pwstmon sets up window sizes appropriate to a particular platform

ctemp = computer;
if(ctemp(1:4)  == 'PCWI' )
   %  ASSUME THIS IS THE DELL NOTEBOOK
   %  fig_pos (TF figure window)
   fig_pos = [ 20,30,450,520];
   fs3 = 11; fs5 = 10; fs6 = 9;
   fs_chname = 10;
   width = 400;
   height_1 = 20;
%  location of lower left corner for Pw_plot dialogue box
   loc = [ 300 50 ];
   nextra = 7;
else
   %  ASSUME THIS IS SUN WORKSTATION 
   %  fig_pos (TF figure window)
   fig_pos = [ 20,30,450,600];
   fs3 = 13; fs5 = 12; fs6 = 11
   fs_chname = 11;
   width = 400;
   height_1 = 25;
%  location of lower left corner for Pw_plot dialogue box
   loc = [ 600 200 ];
   nextra = 7;
end
