      subroutine minpwr(pc,nd,npc,nt,x,pwrmin,nomit)

      integer nd,npc,i,j,k,nt,nomit
      real pwrmin,pwr
      complex pc(nd,2),x(nd,nt)
      do i = 1,nd
         pwr = 0
         do j = 1,npc
            pwr = pwr + pc(i,j)*conjg(pc(i,j))
         enddo
         if(pwr .lt. pwrmin) then
            do k = 1,nt
               x(j,k) = x(j,k)*1.e-10
            enddo
            nomit = nomit+1
         endif
      enddo
      return
      end
