c______________________________________________________________________
c
      subroutine geogcor(u,nt,nvec,nsta,ih,orient,decl,theta,irot)

c           Given : NVEC complex input vectors in U (total # of components
c                is NT, # of stations is NSTA, NCHMX is second dimension of
c                   array orient, IH is array giving index for first H component
c                    for each station (so ih(ista+1)-ih(ista) =nch(ista)))
c                and some directions giving the coordinate system of the components
c                  of the vector:
c                   ORIENT = orientations of field component coordinte system
c                      relative to geomagnetic north
c                   DECL = degrees E of geographic N of the H-axis at each site
c         THis routine rotates U so that on output the components are
c          all expressed in a fixed coordinate system, with x-axis pointing
c                   THETA degrees E of Geographic N
c       (this is the way it works when irot=+1 ; when irot = -1, the coordinate
c            transformation is reversed  )

      real decl(nsta),orient(2,nt),theta
      integer ih(*)
      include 'nstamx.inc'
      complex u(nt,nvec),temp(nchmx,nchmx)

      write(*,*) 'nt,nvec,nsta',nt,nvec,nsta
      write(*,*) 'ih',(ih(k),k=1,nsta+1)

      do ista = 1,nsta
         nch = ih(ista+1)-ih(ista)
         call mkccmt(theta,temp,orient(1,ih(ista)),nch)

         do j=1,nvec
            call rotmult(temp,u(ih(ista),j),nch)
         enddo
      enddo
      return
      end
c_____________________________________________________________________
c
      subroutine mkccmt(theta,u,orient,nch)
c       makes rotation matrix (u) for rotating up to 5 channel data an angle
c       theta (in degrees); the assumed order of the 5 channels is
c       Hx, Hy, Hz, Ex, Ey the rotation direction is clockwise for positive

      include 'nstamx.inc'
      parameter (pc=.0174533)
      real u(nchmx,nchmx),orient(2,*),theta
      do i=1,nch
         do j=1,nch
            u(i,j)=0.0
         enddo
      enddo
      if(nch.eq.5) then
         u(1,1)=cos((orient(1,1)-theta)*pc)
         u(1,2)=cos((orient(1,2)-theta)*pc)
         u(2,1)=sin((orient(1,1)-theta)*pc)
         u(2,2)=sin((orient(1,2)-theta)*pc)
         u(4,4)=cos((orient(1,4)-theta)*pc)
         u(4,5)=cos((orient(1,5)-theta)*pc)
         u(5,4)=sin((orient(1,4)-theta)*pc)
         u(5,5)=sin((orient(1,5)-theta)*pc)
         u(3,3)=1.0
         return
      else if(nch.eq.4) then
         u(1,1)=cos((orient(1,1)-theta)*pc)
         u(1,2)=cos((orient(1,2)-theta)*pc)
         u(2,1)=sin((orient(1,1)-theta)*pc)
         u(2,2)=sin((orient(1,2)-theta)*pc)
         u(3,3)=cos((orient(1,3)-theta)*pc)
         u(3,4)=cos((orient(1,4)-theta)*pc)
         u(4,3)=sin((orient(1,3)-theta)*pc)
         u(4,4)=sin((orient(1,4)-theta)*pc)
         return
      else if(nch.eq.3) then
         u(1,1)=cos((orient(1,1)-theta)*pc)
         u(1,2)=cos((orient(1,2)-theta)*pc)
         u(2,1)=sin((orient(1,1)-theta)*pc)
         u(2,2)=sin((orient(1,2)-theta)*pc)
         u(3,3) = 1.0
         return
      else
         if(nch .ne. 2) then
           write(*,*) 'In mkccmt : Only changing coordinates for ',
     &        ' for first two channels'
         endif
         do i=1,nch
            u(i,i) = 1.0
         enddo
         u(1,1)=cos((orient(1,1)-theta)*pc)
         u(1,2)=cos((orient(1,2)-theta)*pc)
         u(2,1)=sin((orient(1,1)-theta)*pc)
         u(2,2)=sin((orient(1,2)-theta)*pc)
      end if
      end
c________________________________________________________________
c
      subroutine rotmult(u,x,nch)
      include 'nstamx.inc'
      real u(nchmx,nchmx)
      complex work(nchmx),x(nch)

      do i=1,nch
      work(i) = (0.,0.)
         do j=1,nch
            work(i) = work(i) + u(i,j)*x(j)
         enddo
      enddo
      do i = 1,nch
         x(i) = work(i)
      enddo
      return
      end
