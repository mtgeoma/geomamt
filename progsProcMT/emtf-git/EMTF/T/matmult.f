c_______________________________________________________________________
c
        subroutine matmult(a,b,c,nn,np,nq)
 
c       multiplies n x p matrix a by p x q matrix b and puts result in c
 
        real a(nn,np),b(np,nq),c(nn,nq)
 
        do 10 in=1,nn
        do 10 iq=1,nq
        temp=0.
 
           do  5 ip=1,np
5          temp=temp+a(in,ip)*b(ip,iq)
 
10      c(in,iq)=temp
        return
        end
c______________________________________________________________________
c
        subroutine matmultc(a,b,c,nn,np,nq)
 
c       multiplies n x p matrix a by p x q matrix b and puts result in c
c       complex version
 
        complex a(nn,np),b(np,nq),c(nn,nq),temp
 
        do 10 in=1,nn
        do 10 iq=1,nq
        temp=0.
 
           do  5 ip=1,np
5          temp=temp+a(in,ip)*b(ip,iq)
 
10      c(in,iq)=temp
        return
        end
c______________________________________________________________________
c
        subroutine mmacb(a,b,c,nn,np,nq)
 
c       multiplies conjugate transpose of
c         p x n matrix a by p x q matrix b and puts result in c
c       complex version
 
        complex a(np,nn),b(np,nq),c(nn,nq),temp
 
        do 10 in=1,nn
        do 10 iq=1,nq
        temp=0.
 
           do  5 ip=1,np
5          temp=temp+conjg(a(ip,in))*b(ip,iq)
 
10      c(in,iq)=temp
        return
        end
c______________________________________________________________________
c
        subroutine mmabc(a,b,c,nn,np,nq)
 
c       multiplies n x p matrix a by conjugate transpose of
c           q x p matrix b and puts result in c
c       complex version
 
        complex a(nn,np),b(nq,np),c(nn,nq),temp
 
        do 10 in=1,nn
        do 10 iq=1,nq
        temp=0.
 
           do  5 ip=1,np
5          temp=temp+a(in,ip)*conjg(b(iq,ip))
 
10      c(in,iq)=temp
        return
        end
