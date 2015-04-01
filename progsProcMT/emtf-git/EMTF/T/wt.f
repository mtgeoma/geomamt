c
c***********************
c
        subroutine wt(w,psiprime)
 
        if(w.lt.(1.5)) then
           w=1.0
           psiprime=1.0
           return
        end if
 
        w=1.5/w
        psiprime=0.0
        return
        end
