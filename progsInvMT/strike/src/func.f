      real function func(x,p,nu)
      implicit none
      real p, gammq, a, b, x
      integer nu
      a = nu/2.0
      b = x/2.0
      func = gammq(a,b) - (1.0 - p)
      return
      end
