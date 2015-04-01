c
c**********************
c
        subroutine anfld(x,n,kmx,b,c,a,y)
        
        complex x(kmx,2),y(n),b(n,2),c(2),a(2),wrk(2,2)
        
c       this routine takes i complex n-vectors (x) and transforms them
c       to the output complex n-vector (y) which satisfies the i linear 
c       constraints   <y,b(i)> = c(i); a contains the coefficients of the
c       representation of y in terms of the x's 
 
c       special case i = 2

        i = 2

        do 10 ii=1,i
        do 10 jj=1,i
        wrk(ii,jj)=(0.0,0.0)
        
           do 5 kk=1,n
           wrk(ii,jj)=wrk(ii,jj)+x(kk,jj)*b(kk,ii)
5          continue
10      continue

        call cinv2(wrk)
       
        call cmmult(wrk,c,a,2,2,1)
        
        do 15 in=1,n
        y(in)=(0.0,0.0)
        
           do 15 ii=1,i
           y(in)=y(in)+a(ii)*x(in,ii)
15      continue

        return
        end
c______________________________________________________________________
c
      subroutine mkbih(b,nt,ih,ie,nsta,ista)

      complex b(nt,2)
      integer ih(*),ie(*)

      if(ista.ge.1) then
c       use station ISTA as "normal"
         do j = 1,nt
            do k = 1,2
               b(j,k) = (0.,0.)
            enddo
         enddo
         b(ih(ista),1) = (1.,0.)
         b(ih(ista)+1,2) = (1.,0.)
      else
c        use average as normal
         do j = 1,nt
            do k = 1,2
               b(j,k) = (0.,0.)
            enddo
         enddo
         rn = 0.
         do i = 1,nsta 
            if(ih(i)+1 .lt. ie(i)) then
               b(ih(i),1) = cmplx(1.,0.)
               b(ih(i)+1,2) = cmplx(1.,0.)
               rn = rn+1.
            endif
         enddo
         do i = 1,nt
            do k = 1,2
               b(i,k) = b(i,k)/rn
            enddo 
         enddo 
      end if

      return
      end
ccc_____________________________________________________________________
ccc
      subroutine pw_err(a,sig_r,sig_ref,nsta,nt,ih,res_diag,pw,iref,
     &    pw_se)
      complex a(2,2),sig_r(2,2,*),ac(2),temp(2),pw(nt,2),sig_ref(2,2),
     &    temp1
      integer ih(*),iref,nt
      real res_diag(*),pw_se(nt,2),err_ref(2)

      do k = 1,2
         ac(1) = conjg(a(1,k))
         ac(2) = conjg(a(2,k))
         do ista = 1,nsta
           call cmmult(sig_r(1,1,ista),ac,temp,2,2,1)
           err = real(temp(1)*a(1,k)+temp(2)*a(2,k))
           if(ista.eq.iref) err_ref(k) = err
           do i = ih(ista),ih(ista+1)-1
              pw_se(i,k) = res_diag(i)*err
           enddo
         enddo
      enddo

      if(iref .gt. 0 )then
ccc      propogate error in reference component into other components
         do i = 1,nt
            if((i.ne.ih(iref)) .and.(i.ne.ih(iref)+1)) then
                temp(1) = pw(i,1)
                temp(2) = pw(i,2)
                call cmmult(sig_ref,temp,ac,2,2,1)
                ac(1) = conjg(ac(1))
                ac(2) = conjg(ac(2))
                call cmmult(temp,ac,temp1,1,2,1)
                do k = 1,2
                   pw_se(i,k) = pw_se(i,k) + real(temp1)*err_ref(k)
                enddo
            else
                do k = 1,2
                   pw_se(i,k) = 0.0
                enddo
            endif
         enddo
      endif
      do i = 1,nt
         do k = 1,2
            pw_se(i,k) = sqrt(pw_se(i,k))
         enddo
      enddo
      return
      end 
