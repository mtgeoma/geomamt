c
c************************************
c
        subroutine badrec(nd,ndp,idp,id,irec,lrec,lgood)
        
        integer ndp(0:nd),idp(0:nd,2,*)
        logical lgood
        
        lgood = .true.
        
        ndpid = ndp(0)
        do 5 i = 1,ndpid
        if((idp(0,1,i).le.lrec) .and. (idp(0,2,i).ge.irec)) then
           lgood = .false.
           return
           end if
5       continue

        ndpid = ndp(id)
        do 15 i = 1,ndpid
        if((idp(id,1,i).le.lrec) .and. (idp(id,2,i).ge.irec)) then
           lgood = .false.
           return
           end if
15      continue

        return
        end
