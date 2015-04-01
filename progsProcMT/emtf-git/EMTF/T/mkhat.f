c
c********************************
c
        subroutine mkhat(s,nu,h)

        complex s(3),h(3)
        integer nu

        real*4 ss(4),det

        ss(1) = s(1)/nu
        ss(2) = s(3)/nu
        ss(3) = real(s(2))/nu
        ss(4) = aimag(s(2))/nu

        det = ss(1)*ss(2) - ss(3)*ss(3) - ss(4)*ss(4)
          
        h(1) = cmplx((ss(2)/det),0.)
        h(3) = cmplx((ss(1)/det),0.)
        h(2) = cmplx((-ss(3)/det),-ss(4)/det)

        return
        end
