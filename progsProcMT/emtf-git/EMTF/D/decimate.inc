C          include file decimate.inc
C          includes variables used in decimation process
C          and related subroutines

C       routines including decimate.inc:
C       dcimte.f
C       decset.f
C       dnff.f
C       fcorsu.f
C       mk_offst.f
C       mkset.f
C       pterst.f


      include 'params2.inc'

      integer         nfc(ndmx)
      integer         ioffr(ndmx)
      integer         ioffs(ndmx)
      integer         npwmx(ndmx)
      integer         missmx(2,ndmx)
      integer         idec(ndmx)
      integer         npw(nchmx,ndmx)
      integer         nfuse(ndmx)
      integer         nwin(ndmx)
      integer         nd, nch
        
      real            olap(ndmx)
      real            fc(ndmx,nfcmx)
        
      common /DECIMATE/ nfc, ioffr, ioffs, npwmx, missmx, idec, npw,
     &                  nfuse, nwin, nd, nch, olap, fc
