c
c***************************************
c
        subroutine mkset(xx,nchp1,nmax,
     1   ifirst,x,nwmx,iset,id,irec,lrec,rmsval,lgood,roff)

        include 'decimate.inc'

        integer ifirst(nd)
        real xx(nd,nchp1,0:nmax),x(nch,nwmx),temp(50)

        logical lgood

        lgood = .true.

c       find set number

        nmod = 1
        do 5 i = 2,id
5       nmod = nmod * idec(i)

        irec = nint(xx(id,1,ifirst(id)))
        r = (nwin(id)-olap(id))*nmod
        iset = nint( (irec-ioffs(id))/r ) + 1
        roff = r * (iset-1) + ioffs(id) - irec

c     positive roff is equivalent to a fast clock
          
c       copy set into array x; fill in missing data; check to see if
c       too much is missing

          do 30 j = 1,nch
          i1 = nwin(id)-2 + npwmx(id)
          ii1 = ifirst(id)+i1
          ii1 = mod(ii1,nmax)
          lrec = nint(xx(id,1,ii1))

           do 30 i = 0,i1

           k = ifirst(id) + i
           k = mod(k,nmax)

           if (xx(id,j+1,k) .gt. rmsval ) then

c       data point is missing; fill in using a running median of
c       missmx(1,id) points centered at current missing point; this
c       number must be odd

              m = missmx(1,id)/2
              temp(m+1) = xx(id,j+1,k)
              do 4 l = 1,m
                 mk1 = k - m + l - 1
                 mk1 = mod(mk1,nmax)
                 mk1 = mod(mk1+nmax,nmax)
                 temp(l)  = xx(id,j+1,mk1)
                 mk2 = k + l
                 mk2 = mod(mk2,nmax)
                 ml1 = m+l+1
                 temp(ml1) = xx(id,j+1,mk2)
4                continue

              call sort(temp,missmx(1,id))

              nmiss = 0
                 do 6 l = 1,missmx(1,id)
6                if( temp(l).gt.rmsval ) nmiss = nmiss + 1

              if(nmiss .gt. missmx(2,id)) then
c              too many points missing to fill in bad point
 
                 lgood = .false.
                 go to 40
              end if

c       (detail below takes care of the distinction between the definition
c               of median for even or odd n)
              n1 =  (missmx(1,id)+1-nmiss)/2
              n2 = n1
              if(mod(nmiss,2).ne.0) n2 = n2+1
              x(j,i+1) = (temp(n1)+temp(n2))/2.

           else
c       else data point is good; just copy into x

              x(j,i+1) = xx(id,j+1,k)

           end if

30      continue

c       reset pointer for start of next set

40      continue
        rec = iset * ( nwin(id) - olap(id) ) * nmod + ioffs(id)
        ir = nint( (rec - irec)/nmod )
        ifirst(id) = mod(ifirst(id) + ir, nmax)

        return
        end
