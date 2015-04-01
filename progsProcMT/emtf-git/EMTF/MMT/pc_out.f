      subroutine pc_outinit(ioupc,cfout,nsta,nt,npc,nbt,stcor,decl,
     &           orient,ih)
      character*80 cfout
      real stcor(2,nsta),decl(nsta),orient(2,nt)
      integer ih(0:nsta)

c     open file for output of ordered FCs
      open(unit=ioupc,file=cfout,form='unformatted')
      write(ioupc) nsta,nt,npc,nbt
      write(ioupc) ih
      write(ioupc) stcor
      write(ioupc) decl
      write(ioupc) orient
      return
      end
c______________________________________________________________________
c
      subroutine pc_out(ioupc,x,xpc,nt,nd,nfb,idl,iband,ixs,
     &   u,ev,npc,var,period,lraw_out)

      complex  x(nd,nt),u(nt,npc),xpc(nd,npc)
      real ev(npc),var(nt),period
      integer idl,iband(2),ixs(nd)
      logical lraw_out
 
ccc   multiply xpc = x*u
ccc   (modulo a bunch of scaling ....)
      if(.not.lraw_out) then
         call pc_coeff(x,nd,nt,u,npc,var,xpc)
      endif

      write(ioupc) period,nd,nfb,idl,iband
      write(ioupc) ixs
      write(ioupc) var

      if(lraw_out) then
         call wrt_i(ioupc,nt,xpc,ev)
         write(ioupc) x
      else
         write(ioupc) ev
         write(ioupc) u
         write(ioupc) xpc
      endif

      return
      end
c______________________________________________________________________
c
      subroutine pc_coeff(xr,nd,nt,u,npc,var,xpc)
      real var(nt),sc
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
            usc(ip) = usc(ip) + abs(u(j,ip))**2.
         enddo
         usc(ip) = sqrt(usc(ip))
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
ccc_____________________________________________________________________
ccc
      subroutine wrt_i(ioupc,nt,u,ev)
      real ev(nt)
      complex u(nt,nt)
      do i = 1,nt
         do j = 1,nt
            u(i,j) = (0.,0.)
         enddo
         u(i,i) = (1.,0.)
      enddo
      write(ioupc) ev
      write(ioupc) u
      return
      end
