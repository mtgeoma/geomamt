nplots = length(h_plots);
if(nplots > 0) 
   delete(h_plots(nplots));
   if (nplots > 1)
      h_plots = h_plots(1:nplots-1);
   else
      h_plots = [];
   end
end
plot_tf;
