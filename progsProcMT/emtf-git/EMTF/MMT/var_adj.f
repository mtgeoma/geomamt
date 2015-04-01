      subroutine var_adj(nt,a,var,work)
      real var(nt),work(nt,*),rmu,dmu
      complex a(nt,nt)

      parameter (dmu = .9,pmin=.4,pmax=1.0001)

c      write(*,*) 'nt = ',nt
c      write(*,*) 'var =',var

      do k = 1,nt
         work(k,nt+1) = var(k)
      enddo
      rmu = 1.0
10    continue
      do k = 1,nt
         work(k,nt+2) = 1.0
      enddo

      do l = 1,nt
         do k = 1,nt
            if(k.eq.l) then
               work(k,l) = abs(1.-rmu*a(k,k))**2. 
            else
               work(k,l) = abs(rmu*a(k,l))**2.*var(l)/var(k)
            endif
         enddo
      enddo
c      write(*,*) 'I-mB : m =  ',rmu
      do i = 1,nt
c         write(*,'(15f8.4)') (work(i,k),k=1,nt)
      enddo
c      write(*,*)

      iw1 = nt+2
      iw2 = nt+3
      call sgesv(nt,1,work,nt,work(1,iw2),work(1,iw1),nt,info)
c      write(*,*) 'info = ',info
      do i = 1,nt
         if((work(i,iw1).le.pmin).or.(work(i,iw1).gt.pmax)) then
            rmu = rmu*dmu
            go to 10
         endif
      enddo
      do i = 1,nt
         var(i) = var(i)*work(i,iw1)
      enddo
      return
      end
