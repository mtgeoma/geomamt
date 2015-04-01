c
c********************************
c
        subroutine isort(ix,n)
 
c****** subroutine isort sorts array ix  of length n from smallest to
c       largest, in place
c********************************************************************
 
        integer ix(n),itemp
 
        do 50 i=2,n
 
        itemp=ix(i)
 
           do 20 j=1,i-1
 
           if(itemp.le.ix(j)) then
              do 10 k=j,i-1
              l=i+j-k
10            ix(l)=ix(l-1)
              ix(j)=itemp
 
              go to 50
              end if
 
20         continue
 
50      continue
        return
        end
