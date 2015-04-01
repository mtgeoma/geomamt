c
c*********************
c
        subroutine modobs(xr,xo,wtfunc,tfunc,errscale,ru,rt)
 
        complex xr(*),xo,tfunc(3),xp
        real ru,rt
 
        external wtfunc
 
        xp=tfunc(1)*xr(1)+tfunc(2)*xr(2)
 
        pr=real(xp)
        pi=aimag(xp)
        zr=real(xo)
        zi=aimag(xo)
 
        w=abs(zr-pr)/errscale
        rt = w*w
        call wtfunc(w,psiprime)
        ru=psiprime
 
        zr=w*zr+(1.-w)*pr
 
        w=abs(zi-pi)/errscale
        rt = rt + w*w
        call wt(w,psiprime)
        ru=ru+psiprime
 
        zi=w*zi+(1.-w)*pi
 
        xo=cmplx(zr,zi)
 
        return
        end
