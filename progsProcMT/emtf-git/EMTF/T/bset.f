        subroutine bset(nb,ibandlim,idl,period,nd,nbmax,
     1  ierr,samprate,npts,cfile)
     
        integer ibandlim(2,nbmax),idl(nbmax),npts(nd)
        
        real period(nbmax),samprate(nd)
        character*40 cfile
        
        open(unit=59,file=cfile,status = 'old', err = 1000)
        
        
        ierr = 0
        read(59,*) nb
        if(nb.gt.nbmax) then
           print*,'error: nb exceeds nbmax; nb = ',nb,'nbmax=',nbmax
           ierr = 1
           return
           end if
           
        do 10 ib=1,nb
        read(59,*) idl(ib),ibandlim(1,ib),ibandlim(2,ib)
        if((idl(ib).gt.nd).or.(idl(ib).lt.0)) then
           print*,'error: decimation level invalid; idl=',idl(ib)
           ierr = 2
           return
           end if
           
        period(ib)=samprate(idl(ib))*npts(idl(ib))*2./
     1                            (ibandlim(1,ib)+ibandlim(2,ib))
     
     
10      continue
        close(59)
        return

 1000   print*,' ERROR OPENING BS_NOD FILE: ', cfile
        end
