c_____________________________________________________________
c
        subroutine stack(x,nch,s)
 
        complex x(*),s(*)
 
        ij=0
           do 140 i=1,nch
           do 140 j=1,i
           ij=ij+1
140        s(ij)=s(ij)+x(i)*conjg(x(j))
        return
        end
c______________________________________________________________
c
        subroutine wstack(x,w,nch,s)
 
        complex x(*),s(*)
c***************adds weighted cross product to  SDM
 
        ij=0
           do 140 i=1,nch
           do 140 j=1,i
           ij=ij+1
140        s(ij)=s(ij)+x(i)*conjg(x(j))*w
        return
        end

