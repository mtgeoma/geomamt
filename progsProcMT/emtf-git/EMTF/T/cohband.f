      subroutine cohband(iband,ibcoh,nfreqmax,nfb,ncb,cohw,icohbm)

c      Given limits of narrow processing band (IBAND) returns wider
c      limits of coherence sort band (IBCOH);
c      COHW is input parameter controlling relative width of
c      coherence presort band

      integer iband(2),ibcoh(2),nfreqmax,nfb,ncb,icohbm
      ibcoh(1) = nint(iband(1)/cohw)
      ibcoh(1) = max(ibcoh(1),1)
      ibcoh(2) = nint(iband(2)*cohw)
      ibcoh(2) = min(ibcoh(2),nfreqmax)
      nfb = iband(2)-iband(1)+1
      ncb = ibcoh(2)-ibcoh(1)+1
      if(ncb.lt.icohbm) then
         c = sqrt(float(icohbm))
         ibcoh(1) = nint(iband(1)-c)
         ibcoh(1) = max(1,ibcoh(1))
         ibcoh(2) = ibcoh(1)-1+icohbm
         ncb = ibcoh(2)-ibcoh(1)+1
      end if
      return
      end
