c
c******************************
c
      subroutine mkrhat(edfl1,nf,x,nch,h,ierr,lw)

      integer*2 icheck(5000)
      complex s(3),x(nch,nf),h(3)
      logical lw(nf)

      edfl2 = edfl1/2.
      ierr = 0
      no = 0
      nos = 0

      do 1 i = 1,3
1     s(i) = (0.,0.)

      do 5 i = 1,nf
      lw(i) = .false.
5     call stack(x(1,i),2,s)
      iter = 0

10    continue
      nt = nf - no
      call mkhat(s,nt,h)

      ncheck = 0

         do 20 i = 1,nf

         if(lw(i)) go to 20

         edf = x(1,i)*conjg(x(1,i))*h(1) + x(2,i)*conjg(x(2,i))*h(3)
     &             +2.*real(conjg(x(2,i))*x(1,i)*h(2))

         if(edf .gt. edfl1) then
            call unstack(x(1,i),2,s)
            lw(i) = .true.
    
            no = no + 1         
         else if(edf .gt. edfl2) then
            ncheck = ncheck + 1
            
            icheck(ncheck) = i
         end if

20       continue
     
      nt = nf - no
      if(nt .lt. 2) then
         ierr = -1
         return
      end if

      call mkhat(s,nt,h) 
                     
      if(ncheck .gt. 0) then
         nm = 0

25       continue

         nn = 0
         do 30 i = 1,ncheck
 
         ii = icheck(i)
         if(lw(ii)) go to 30

        edf = x(1,ii)*conjg(x(1,ii))*h(1) + x(2,ii)*conjg(x(2,ii))*h(3)
     &             +2.*real(conjg(x(2,ii))*x(1,ii)*h(2))
                                                 
         if(edf .gt. edfl1) then
            nn = nn + 1
            nm = nm + 1
            no = no + 1
                    
            call unstack(x(1,ii),2,s)

            lw(ii) = .true.

         end if

30       continue

         if(nm .eq. ncheck) then
            nos = no
            go to 10
         end if
          
         if(nn .gt. 0) then
            
            nt = nf - no
            if(nt .lt. 2) then
               ierr = -1
               return
            end if

            call mkhat(s,nt,h)
            go to  25

         end if
         
         return

      else if(no .eq. nos) then
         return

      else
         nos = no
         go to 10

      end if

      end 
