c
c**********************************
c
        subroutine dcfilt(temp,rmsval,nfc,fc,nd,ndmx,id2,xd)
 
        real fc(ndmx,*),temp(*)
 
        integer nfc(nd)
 
        maxnmiss = nfc(id2)-1
        nfcid = nfc(id2)
        nmiss = 2*nfcid-1
        wt = 0.
        t = 0.
 
           do 10 i = 1,nfcid
           if( temp(i) .lt. rmsval) then
              wt = wt + fc(id2,nfcid - i + 1)
              t = t + fc(id2,nfcid - i + 1) * temp(i)
              nmiss = nmiss - 1
              end if
10         continue
 
           do 20 i = 2,nfcid
           ii = nfcid+i-1
           if( temp(ii) .lt. rmsval) then
              wt = wt + fc(id2,i)
              t = t + fc(id2,i) * temp(ii)
              nmiss = nmiss - 1
              end if
20         continue
 
        if (nmiss .gt. maxnmiss) then
           xd = rmsval + 10000.
           return
        end if
 
        xd = t/wt
        return
        end
