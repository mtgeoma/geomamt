 
c***************************************                
        subroutine mult(r,t1,t2,t3,t4,det)
        complex t1,t2,t3,t4,det,r
        r=t1*t2-t3*t4
        r=r/det
        return
        end
