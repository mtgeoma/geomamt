      subroutine pc_coeff(xr,nd,nt,u,npc,var,xpc)
      real var(nt)
      complex u(nt,npc),xr(nd,nt),xpc(nd,npc)
      include 'nstamx.inc'
      integer i1(ntmx)
      real usc(ntmx)

ccc   given NT total data channels with estimated noise
ccc   variances of VAR, and NPC principal components U,
ccc   from the "PC coefficients" ...

      do i = 1,nt
         i1(i) = i
      enddo 

ccc   scale PCs to allow for noise variations
      do j =  1,nt
         sc = sqrt(var(j))
         do ip = 1,npc
            u(j,ip) = u(j,ip)/sc
         enddo
      enddo
      do ip = 1,npc
         usc(ip) = 0.0
         do j = 1,nt
            usc(ip) = usc(ip) + u(j,ip)**2.
         enddo
      enddo
      do j =  1,nt
         sc = sqrt(var(j))
         do ip = 1,npc
            u(j,ip) = u(j,ip)/(sc*usc(ip))
         enddo
      enddo
    
      call lc_dat(xr,nd,nt,u,nt,npc,i1,xpc)

ccc   return PCs to original scaling
      do ip = 1,npc
         do j =  1,nt
            u(j,ip) = u(j,ip)*var(j)*usc(ip)
         enddo
      enddo
      return
      end
