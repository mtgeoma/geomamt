      subroutine refproj(n1,i1,nt,npcmx,u,b,nf,eval)

      complex u(n1,npcmx),b(nt,2),uref(10,2),utemp(2)
      real eval(*)
      integer i1(n1)

      npcu = n_eig_sig(nf,eval,n1,npcmx)
      if(npcu.eq.2) return 
      do i = 1,npcu
         do j = 1,2
            uref(i,j) = (0.,0.)
            do l = 1,n1
               uref(i,j) = uref(i,j) + conjg(u(l,i))*b(i1(l),j)
            enddo
         enddo
      enddo

      do i = 1,n1
         do j = 1,2
            utemp(j) = (0.,0.)
            do l = 1,npcu
               utemp(j) = utemp(j) + u(i,l)*uref(l,j)
            enddo
         enddo
         do j = 1,2
            u(i,j) = utemp(j)
         enddo
      enddo
      return 
      end
ccc_____________________________________________________________________
ccc
      subroutine refproj_cor(n1,i1,nt,u,pctf,sig_s,nu)

      complex sig_s(2,2),pctf(2,nt),u(n1,2),uc
      integer nu,l,j,i,i1(n1)

      do i = 1,n1
         do j = 1,2 
            uc = 0.0
            do l = 1,2
               uc = uc+sig_s(l,j)*u(i,l)
            enddo
            pctf(j,i1(i)) = pctf(j,i1(i)) - uc
         enddo
      enddo
      return
      end
