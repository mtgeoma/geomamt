C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     DHAT (Nx1) = DATA - FF + CMGT*log10(RH)
C     CDHAT = Cd^{-1/2}*DHAT;
     
      SUBROUTINE CompDHAT(Nzb,Ny,NMode,NRes,NPer,NSta,NNT,
     >           Period,SknDepth,StaPos,
     >           DatRes,DatErr,FF,DatInx,SenInx,
     >           CRho,PriRho,DFStatus,
     >           CmHtIndx,CmHt,CDhat)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny,NMode,NRes(*),NPer(*),NSta(*),NNT(*),DFStatus(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  SknDepth(NMODMX,NPERMX)
      REAL*8  StaPos(NMODMX,NSTAMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  FF(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  CDhat(*),CmHt(MM0MX,LL0MX),Period(NMODMX,NPERMX)
      REAL*8  PriRho(NZ0MX,NY0MX),CRho(NZ0MX,NY0MX)
      
      INTEGER jj,iy,iz,idx,im,ir,ip,is,mmt
      REAL*8  dRho(MM0MX),rho_dot_g
      REAL*8  Dhat,G(MM0MX)

      REAL*8  DDOT

      mmt = Nzb*Ny
      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = jj+1
          IF (DFStatus(jj).EQ.0) THEN
            dRho(jj) = D0
          ELSE
            dRho(jj) = DLOG10(CRho(iz,iy)) - DLOG10(PriRho(iz,iy))
          ENDIF
        ENDDO ! iz
      ENDDO ! iy

      idx = 0
      DO im = 1,NMODE
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im)
              IF (DatInx(im,ir,ip,is).EQ.1) THEN
                idx = idx + 1
                IF (SenInx(im,ir,ip,is).EQ.1) THEN
                  CALL ExtractSens(0,im,ir,ip,is,mmt,
     >                 SenInx,CmHtIndx,CmHt,DFStatus,D1,G)
                  rho_dot_g = DDOT(mmt,dRho,1,G,1)
                ELSE
                  CALL InterpolateSens(im,ir,ip,is,Nzb,Ny,
     >                 NMode,NRes,NPer,NSta,Period,SknDepth,StaPos,
     >                 DatInx,SenInx,CmHtIndx,CmHt,DFStatus,G)
                  rho_dot_g = DDOT(mmt,dRho,1,G,1)
                ENDIF ! SenInx

                Dhat = DatRes(im,ir,ip,is)-FF(im,ir,ip,is) 
     >               + rho_dot_g
                IF (NNT(3).LT.NNT(2)) THEN
                  CDhat(idx) = Dhat/DatErr(im,ir,ip,is) 
                ENDIF
                IF (NNT(3).EQ.NNT(2)) THEN
                  CDhat(idx) = Dhat
                ENDIF
              ENDIF ! DatInx
            ENDDO ! is
          ENDDO ! ip
        ENDDO ! ir
      ENDDO ! im

      RETURN
      END ! CompDHAT()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C<<<< FORMING REPRESENTATIVE MATRIX (R) (LxL)
C<<<< R = REPRESENTATIVE MATRIX AND H = SENSITIVITY MATRIX  (LxM)
C<<<< TO SAVE MEMORY SPACE, WE SAVE H IN THE H'  (MxL) FORMAT, THEN
C<<<< AFTER COMPUTING R, WE DON'T NEED H.  WE WILL SAVE CM*H'
C<<<< IN THAT SAME VARIABLE (CmHt)

C<<<< DETAILS
C<<<<     R = H*CM*H' = H*(D0 + D^N)*H' = H*D0*H' + H*D^N*H'
C<<<<     D0 IS DEFINED AS ONES(M,M) MATRIX
C<<<< [1] H*D0*H' = (H*D0)*(D0*H')
C<<<<     DEFINE DD (1xL) =  D0*H' (MxL) = sum(H') OF EACH COLUMN
C<<<<     THE OTHER ROW OF DD (M ROWS) ARE THE SAME AS THE FIRST ONE
C<<<< [2] H*D^N*H' = (H*D^N/2)*(D^N/2*H') = (D^N/2*H')'*(D^N/2*H')
C<<<<     REPLACE H' WITH D^(N/2)*H'


      SUBROUTINE CompSubRepm(LOGFILE_SCREEN,
     >           Nzb,Ny,Dzb,Dy,NMode,NRes,NPer,NSta,Period,
     >           NNT,DatInx,SenInx,DFStatus,PDifTime,HDiff,VDiff,
     >           CmNorm,CmHtIndx,CmHt,DD,REPM)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      INTEGER DFStatus(*),LOGFILE_SCREEN
      INTEGER PDifTime
      REAL*8  Dzb(*),Dy(*),CmNorm(*)
      REAL*8  CmHt(MM0MX,LL0MX),Period(NMODMX,NPERMX)
      REAL*8  DD(*),REPM(*)
      REAL*8  HDiff(2,NZ0MX,NY0MX),VDiff(2,NZ0MX,NY0MX)

      INTEGER ii,jj,kk,mmt,idx
      REAL*8  G(MM0MX),dotproduct

      REAL*8  sumx

      mmt = Ny*Nzb

C     Calculate DD (Lx1) = D0*H' = sum(H') and D^(N/2)*H'
      DO ii = 1,NNT(3)
        sumx = D0
        DO jj = 1,mmt
         sumx = sumx + CmHt(jj,ii)
        ENDDO
        DD(ii) = sumx

        IF (DD(ii).EQ.D0) THEN
          WRITE(6,5000) 
          WRITE(6,5010) 
          STOP
        ENDIF
      ENDDO

      IF (LOGFILE_SCREEN.EQ.1) WRITE(6,6000)  
      DO ii = 1,NNT(3)
       DO jj = 1,mmt
         G(jj) = CmHt(jj,ii)*DFStatus(jj)/DSQRT(CmNorm(jj))
       ENDDO
       CALL SolveDiff(1,PDifTime,HDiff,VDiff,Nzb,Ny,Dzb,Dy,DFStatus,G)
       DO jj = 1,mmt
         CmHt(jj,ii) = DFStatus(jj)*G(jj)
       ENDDO
      ENDDO


      IF (LOGFILE_SCREEN.EQ.1) WRITE(6,6010)  
C     Repm = CmHt'*CmHt + DD'*DD 
c     Dot Product Matrix Multiplication
      do ii = 1,NNT(3)
        do kk = 1,ii
          dotproduct = D0
          do jj=1,mmt
            dotproduct = dotproduct + CmHt(jj,ii)*CmHt(jj,kk)
          enddo
          idx = kk+(ii-1)*ii/2
          Repm(idx) = dotproduct + DD(kk)*DD(ii)
        enddo ! kk
      enddo ! ii


      IF (LOGFILE_SCREEN.EQ.1) WRITE(6,6020)  

C     Calculate Cm*H' = (D0+D^N)*H' = D0*H' + D^(N/2)*CmHt
      DO ii = 1,NNT(3)
       DO jj = 1,mmt
         G(jj) = CmHt(jj,ii)*DFStatus(jj)
       ENDDO

       CALL SolveDiff(2,PDifTime,HDiff,VDiff,Nzb,Ny,Dzb,Dy,DFStatus,G)
       DO jj = 1,mmt
        CmHt(jj,ii) =  DFStatus(jj)*(G(jj)/DSQRT(CmNorm(jj)) +DD(ii))
       ENDDO
      ENDDO


5000  FORMAT('!!! ATTENTION, ERROR WHILE COMPUTING CROSS-PRODUCT !!!')
5010  FORMAT('!!!  Please, check data input and restart          !!!')

6000  FORMAT('*** forward smoothing sensitivity matrix ***')
6010  FORMAT('*** forming cross-product matrix ***')
6020  FORMAT('*** backward smoothing sensitivity matrix ***')


      RETURN
      END ! CompSubRepm()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     extract sensitivity value from CmHt at a given irl, ipl and isl.

    
      SUBROUTINE ExtractSens(emode,iml,irl,ipl,isl,mmt,
     >           SenInx,CmHtIndx,CmHt,DFStatus,wt,G)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER iml,irl,ipl,isl,mmt,emode,DFStatus(*)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX),G(*),wt

      INTEGER igx,jj,iss,is

      IF (SenInx(iml,irl,ipl,isl).EQ.0) THEN
         WRITE(6,*) '!!! ATTENTION, ERROR IN SENS. CALCULATION !!!'
         WRITE(6,*) 
     >   '!!! Please check sensitivity inclusion and restart !!!'
         STOP
      ENDIF

C     emode = 0  : generate G in every blocks
C     emode = 1  : generate G in every blocks and multiply with DFStatus
C     emode = 2  : generate G in every blocks and multiply with constant

      igx = CmHtIndx(iml,irl,ipl)
      iss = 0
      DO is = 1,isl
        IF (SenInx(iml,irl,ipl,is).EQ.1) iss = iss + 1
      ENDDO

      IF (emode.EQ.1) THEN
       DO jj = 1,mmt
         G(jj) = CmHt(jj,igx+iss)*DFStatus(jj)
       ENDDO
       GOTO 100
      ENDIF

      IF (emode.EQ.2) THEN
       DO jj = 1,mmt
         G(jj) = CmHt(jj,igx+iss)*wt
       ENDDO
       GOTO 100
      ENDIF

      IF (emode.EQ.0) THEN
       DO jj = 1,mmt
        G(jj) = CmHt(jj,igx+iss)
       ENDDO
      ENDIF

100   CONTINUE


      RETURN
      END ! ExtractSens()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     approximate sensitivity value at a given irl, ipl and isl.

      SUBROUTINE InterpolateSens(iml,irl,ipl,isl,Nzb,Ny,
     >           NMode,NRes,NPer,NSta,Period,SknDepth,StaPos,
     >           DatInx,SenInx,CmHtIndx,CmHt,DFStatus,
     >           G)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER iml,irl,ipl,isl,Nzb,Ny
      INTEGER NMode,NRes(*),NPer(*),NSta(*),DFStatus(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX),G(*),Period(NMODMX,NPERMX)
      REAL*8  SknDepth(NMODMX,NPERMX)
      REAL*8  StaPos(NMODMX,NSTAMX)

      INTEGER mmt,np1,np2,np3,np4,jj,npp,nss
      INTEGER io(8),jo(8),ii
      REAL*8  go(MM0MX),ds(8),dp(8),wt(8),wt1,wt2,wtt,wf
      REAL*8  sumds,sumdp


      mmt = Nzb*Ny

      IF (DatInx(iml,irl,ipl,isl).EQ.0) THEN
         WRITE(6,*) '!!! ATTENTION, ERROR IN SENS. CALCULATION !!!'
         WRITE(6,*) 
     >   '!!! Please check sensitivity inclusion and restart !!!'
         STOP
      ENDIF

      np1 = NMODMX
      np2 = NRESMX
      np3 = NPERMX
      np4 = NSTAMX

      npp = NPer(iml)
      nss = NSta(iml)

C     approximate Sens from nearest Sensitivity value 
      CALL FindNearest(iml,irl,ipl,isl,npp,nss,DatInx,SenInx,
     >     SknDepth,Period,StaPos,np1,np2,np3,np4,io,jo,ds,dp)
      sumds = D0
      sumdp = D0
      DO ii = 1,8
        IF (io(ii).NE.0) THEN
          sumds  = sumds + ds(ii)
          sumdp  = sumdp + dp(ii)
        ENDIF
      ENDDO

      wtt = D0
      DO ii = 1,8
        wt(ii)   = D0
        IF (io(ii).NE.0) THEN
          wt1    = DEXP(-ds(ii)/sumds)
          wt2    = DEXP(-dp(ii)/sumdp)
          wt(ii) = wt1*wt2
          wtt    = wtt + wt(ii) 
        ENDIF
      ENDDO

      CALL ConstantVectorR8(G,mmt,D0)
      DO ii = 1,8
        IF (io(ii).NE.0) THEN
          wf = wt(ii)/wtt
          IF (wf.GE.1E-02) THEN
            CALL ExtractSens(2,iml,irl,io(ii),jo(ii),mmt,
     >           SenInx,CmHtIndx,CmHt,DFStatus,wf,go)
            DO jj = 1,mmt
              G(jj) = G(jj) + go(jj)
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END ! InterpolateSens()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Repm <= R*Repm*R'   Where R is from QR-decomposition

      SUBROUTINE MulRepm(NNT,BB,Repm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
      INTEGER NNT(*)
      REAL*8  BB(NN0MX,LL0MX),Repm(*)
    
      REAL*8  R1(LLHMX)

      INTEGER jj,ii,kk,ll,rr
      REAL*8  sumx

      DO ii = 1,NNT(3)
        DO jj = ii,NNT(3)
          sumx = D0
          DO ll = ii,jj-1
            kk = ll + (jj-1)*jj/2 
            sumx = sumx + BB(ii,ll)*Repm(kk)
          ENDDO ! ll
          DO ll = jj,NNT(3)
            kk = jj + (ll-1)*ll/2 
            sumx = sumx + BB(ii,ll)*Repm(kk)
          ENDDO ! ll
          rr     = ii + (jj-1)*jj/2
          R1(rr) = sumx
        ENDDO ! jj
      ENDDO ! ii

C     Repm <= R1*R'
      DO ii = 1,NNT(3)
        DO jj = ii,NNT(3)
          sumx = D0
          DO ll = jj,NNT(3)
            rr   = ii + (ll-1)*ll/2 
            sumx = sumx + R1(rr)*BB(jj,ll)
          ENDDO ! ll
          kk = ii + (jj-1)*jj/2 
          Repm(kk) = sumx
        ENDDO ! jj
      ENDDO ! ii

      RETURN
      END ! MulRepm

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
