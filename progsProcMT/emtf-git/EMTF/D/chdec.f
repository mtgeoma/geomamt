
        subroutine chdec(b,a,n,ier)
 
c       cholesky decomposition of real symmetric matrix; b is input
c       matrix of order n in symmetric storage mode (requires array of
c       length n*(n+1)/2 in calling routine); a is output (in same form)
c       of lower triangular matrix
 
        real a(*),b(*)
 
        ier=0
 
        do 30 i=1,n
        ii=(i-1)*i/2
 
           do 15 j=1,i-1
           jj=(j-1)*j/2
           t=b(ii+j)
           i1=ii
 
              do 10 k=1,j-1
              i1=i1+1
              jj=jj+1
10            t=t-a(i1)*a(jj)
 
15         a(ii+j)=t/a(jj+1)
 
        i1=ii
        t=b(ii+i)
 
           do 20 k=1,i-1
           i1=i1+1
20         t=t-a(i1)*a(i1)
 
        if(t.le.0) then
           ier=-1
           return
           end if
 
30      a(ii+i)=sqrt(t)
 
        return
        end
