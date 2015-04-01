c
c************************************
c 
        subroutine demean(x,m,n)
        real x(m,n),xbar(25)
 
        do 5 j=1,m
        xbar(j)=0.0
        do 2 i=1,n
        xbar(j)=xbar(j)+x(j,i)
2       continue
        xbar(j)=xbar(j)/n
5       continue
 
        do 10 j=1,m
        do 10 i=1,n
        x(j,i)=x(j,i)-xbar(j)
10      continue
 
        return
        end
