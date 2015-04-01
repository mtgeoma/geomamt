      real function chiprob(nu,p)
      implicit none
      real p, x1, x2, xacc, rtbis, func, test
      integer nu
	external func
      x1 = 1.0
      test = func(x1,p,nu)
      x2 = 1000.0
      xacc = 0.00001
      chiprob = rtbis(func,x1,x2,xacc,p,nu)
      return
      end
