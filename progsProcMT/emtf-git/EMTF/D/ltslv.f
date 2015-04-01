c
c**************************
c
        subroutine ltslv(a,n,b,m)
 
c       linear equation solver for real lower triangular system
c       the equation to solve is ax=b where a is n x n and
c       b is n x m; a is lower triangular in symmetric storage mode
c       result is returned in b
 
        real a(6),b(n,m)
 
        do 50 i=1,m
 
           do 20 j=1,n
 
           j1=j*(j-1)/2+1
           jl=j1+j-2
 
             temp=0.0
             do 10 k=j1,jl
             temp=temp+a(k)*b(k-j1+1,i)
10           continue
 
           b(j,i)=(b(j,i)-temp)/a(jl+1)
20         continue
 
50      continue
        return
        end
