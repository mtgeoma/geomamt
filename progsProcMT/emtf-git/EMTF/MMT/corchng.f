      subroutine ccmat(cc,orient1,theta)

ccc   makes 2x2 matrix cc to convert a vector pair of field measurements
ccc   made in directions orient1(1),orient1(2) into a standard (right handed,
ccc   with z pointing down) system of orthogonal coordinates with the x-axis
ccc   pointing in the direction theta.  All angles are given in degrees
ccc   relative to a common reference direction (usually geographic or
ccc   geomagnetic north).  Positive angles correspond to clockwise rotation
ccc   This works even if the original measurements are not orthogonal

      real pi
      parameter(pi = 3.141592654)
      real cc(2,2),theta,orient1(2),pc,det,temp

      pc = pi/180
ccc   matrix to convert vector components from standard (rotated) orthogonal
ccc   coordinate system to arbitrary (not necessarily orthogonal)
ccc   measurement coordinate system

      cc(1,1)=cos((orient1(1)-theta)*pc)
      cc(2,1)=cos((orient1(2)-theta)*pc)
      cc(1,2)=sin((orient1(1)-theta)*pc)
      cc(2,2)=sin((orient1(2)-theta)*pc)

ccc   we want the inverse of this ... Use Kramers rule

      det = cc(1,1)*cc(2,2)-cc(1,2)*cc(2,1)
      temp = cc(1,1)
      cc(1,1) = cc(2,2)/det
      cc(2,2) = temp/det
      cc(1,2) = -cc(1,2)/det
      cc(2,1) = -cc(2,1)/det

      return
      end 
ccc_____________________________________________________________________
ccc
      subroutine cctf(u,nt,nvec,ipair,npair,orient,theta)
ccc   change coordinates of components of complex vectors in U(nt,nvec)
ccc   for NPAIR pairs of channels listed in array IPAIR
ccc   orientations of all channels are in ORIENT as usual
ccc   and the x-axis of the new (orthogonal) output coordinate system
ccc   is in the direction THETA

      real orient(2,nt),theta,cc(2,2),orient1(2)
      integer ipair(2,npair),nt,nvec,npair,i,k
      complex u(nt,nvec),u1(2)

      do i =1,npair
         orient1(1) = orient(1,ipair(1,i))
         orient1(2) = orient(1,ipair(2,i))
         call ccmat(cc,orient1,theta)
         do k = 1,nvec
            u1(1) = cc(1,1)*u(ipair(1,i),k)+cc(1,2)*u(ipair(2,i),k)
            u1(2) = cc(2,1)*u(ipair(1,i),k)+cc(2,2)*u(ipair(2,i),k)
            u(ipair(1,i),k) = u1(1)
            u(ipair(2,i),k) = u1(2)
         enddo
      enddo 

      return
      end
c______________________________________________________________________
c
      subroutine cc_geog(u,nt,nvec,nsta,ih,ie,orient,theta,e_and_h)

ccc   convert pairs of Hx/Hy and Ex/EY eigenvector components
ccc      into a standard coordinate system  (replaces old geogcor)
ccc    THIS assumes that orient is expressed relative to the same
ccc    FIXED direction as theta 
ccc   If e_and_h is .true. , rotate both E and H channels (so ie needs
ccc    to be defined)  If .false.  ie is not referred to

      real orient(2,nt),theta
      integer ih(*),ie(*),ipair(2,500),npair,nt,nvec,ista,nsta
      complex u(nt,nvec)
      logical e_and_h

      npair = 0
      do ista = 1,nsta
         if(ie(ista) - ih(ista) .ge. 2) then
ccc         at least two H at station ista ... assume first two are Hx/Hy
            npair = npair + 1
            ipair(1,npair) = ih(ista) 
            ipair(2,npair) = ih(ista)+1 
         endif
         if(e_and_h) then
            if(ih(ista+1) - ie(ista) .ge. 2) then
ccc            at least two E at station ista ... assume first two are Ex/Ey
               npair = npair + 1
               ipair(1,npair) = ie(ista) 
               ipair(2,npair) = ie(ista)+1 
            endif
         endif
      enddo
      call cctf(u,nt,nvec,ipair,npair,orient,theta)

      return
      end
ccc_____________________________________________________________________
ccc
      subroutine ccsdm(s,nt,ipair,npair,orient,theta)
ccc   change coordinates of elements of complex SDM in S
ccc   for NPAIR pairs of channels listed in array IPAIR
ccc   orientations of all channels are in ORIENT as usual
ccc   and the x-axis of the new (orthogonal) output coordinate system
ccc   is in the direction THETA

      real orient(2,nt),theta,cc(2,2),orient1(2)
      integer ipair(2,npair),ii(2),sindex,nt,npair,i,k,l
      complex S(*),s1(2),sc1(2)

      do i =1,npair
         orient1(1) = orient(1,ipair(1,i))
         orient1(2) = orient(1,ipair(2,i))
         call ccmat(cc,orient1,theta)

ccc      multiply on left by 2x2 transformation matrix for each pair
ccc      loop over columns  (NOTE : loop needs to go backwards
ccc       wto keep from overwriting input SDM prematrurely
         do k = nt,1,-1
            do l = 1,2
               ii(l) = ipair(l,i)
               if(ii(l) .ge. k) then
                  s1(l) = s(sindex(ii(l),k))
               else
                  s1(l) = conjg(s(sindex(k,ii(l))))
               endif 
            enddo
            
            sc1(1) = cc(1,1)*s1(1)+cc(1,2)*s1(2)
            sc1(2) = cc(2,1)*s1(1)+cc(2,2)*s1(2)
            do l = 1,2
               if(ii(l) .ge. k) then
                  s(sindex(ii(l),k)) = sc1(l)
               endif 
            enddo
         enddo

ccc      multiply on right by transposes of 2x2 transformation matrix 
ccc      for each pair 
ccc      loop over rows
         do k = 1,nt
            do l = 1,2
               ii(l) = ipair(l,i)
               if(k.ge. ii(l)) then
                  s1(l) = s(sindex(k,ii(l)))
               else
                  s1(l) = conjg(s(sindex(ii(l),k)))
               endif 
            enddo
            
            sc1(1) = cc(1,1)*s1(1)+cc(1,2)*s1(2)
            sc1(2) = cc(2,1)*s1(1)+cc(2,2)*s1(2)
            do l = 1,2
               if(k.ge.ii(l)) then
                  s(sindex(k,ii(l))) = sc1(l)
               endif 
            enddo
         enddo
      enddo 

      return
      end
