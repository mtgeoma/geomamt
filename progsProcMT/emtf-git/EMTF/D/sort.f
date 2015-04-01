c
c***************************************
c 
        subroutine sort(x,n)
 
c****** subroutine sort sorts array x of length n from smallest to
c       largest, in place
c********************************************************************
 
        real x(n)
 
        do 50 i=2,n
        temp=x(i)
 
            do 20 j=1,i-1
 
           if(temp.le.x(j)) then
 
              do 10 k=j,i-1
              l=i+j-k
10            x(l)=x(l-1)
 
              x(j)=temp
              go to 50
              end if
 
20         continue
 
50      continue
        return
        end
