%**********************************************************************************
function clear_win
%clear_win: Deletes all nmatlab graphics windows open on screen
hwin = get(0,'Children')
for win=hwin'
   delete(win)
end
end
