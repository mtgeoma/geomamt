c
c*****************************
c
        subroutine rdcndwt(w,psiprime)
 
        parameter (u0=2.8)
 
        if (w.gt. 4) then
        w=0.0
        psiprime=0.0
        return
        end if
 
        t=(-exp(u0*(w-u0)))
        u=w
        w=exp(t)
        t=u0*t*u
 
        if((t+1.).lt.(1.e-20)) then
           psiprime=0.
           return
        end if
 
        psiprime=(w*(1.+t))
 
        return
        end
