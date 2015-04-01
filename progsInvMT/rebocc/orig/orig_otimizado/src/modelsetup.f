C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Locate Position of the station on the model.
C     If the station is not locate at the boundary of the model,
C     the program will automatically generate a pseudo-block.
C
      SUBROUTINE LocateModelPosition(NMode,NPer,NSta,Period,
     >           LStart,StaLoc,Nzb,Nz,Ny0,Dy0,
     >           Ny,Dy,Cy,YDis,IniRho,PriRho,StaPos,StaNod,
     >           ModGrd,DFStatus)
      INCLUDE 'parameter.h'

C        IO Variables
      INTEGER NMode,NPer(*),NSta(*),Ny0,Ny,Nzb,Nz
      INTEGER ModGrd(NZ0MX,NY0MX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER DFStatus(*)
      REAL*8  Dy0(*),Dy(*),YDis(*),Cy(*)
      REAL*8  LStart(*)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  StaLoc(NMODMX,NSTAMX)
      REAL*8  IniRho(NZ0MX,NY0MX),PriRho(NZ0MX,NY0MX)
      REAL*8  StaPos(NMODMX,NSTAMX)

C       Local Variables
      INTEGER nst,im,is,js,iy,iyy
      INTEGER pnode(NMODMX*NSTAMX),oldny
      INTEGER mgrdtmp(NZ0MX,NY0MX),PStart(NMODMX)
      REAL*8  spos(NMODMX*NSTAMX),staspace
      REAL*8  dytmp(NY0MX),rhotmp(NZ0MX,NY0MX)

C     Assign variables
      Ny   = Ny0
      CALL CopyVectorR8(1,Ny0,Dy0,1,Ny0,Dy)

C     Locate position of stations on the initial model
      DO im = 1,NMODE
        LStart(im) = LStart(im) + StaLoc(im,1)
        PStart(im) = 0
        DO iy = 2,Ny0+1
         IF (LStart(im).EQ.YDis(iy)) THEN
            PStart(im) = iy
            GOTO 100
         ELSE
          IF ((LStart(im).GT.YDis(iy-1)).AND.
     >        (LStart(im).LT.YDis(iy))) THEN
             PStart(im) = iy-1
             WRITE(6,*) 
             WRITE(6,1000) im
             WRITE(6,1010) LStart(im),YDis(iy-1)
             WRITE(99,1000) im
             WRITE(99,1010) LStart(im),YDis(iy-1)
             GOTO 100
          ENDIF
         ENDIF
        ENDDO

   
100     CONTINUE
        StaPos(im,1) = YDis(PStart(im))
        DO is = 2,NSta(im)
          staspace = StaLoc(im,is) - StaLoc(im,is-1)
          StaPos(im,is) = StaPos(im,is-1) + staspace
        ENDDO
      ENDDO ! im

C     Determine total no. of distinct stations of every modes
C     and sort them out in ascending order
      CALL FindDistinctStation(NMode,NSta,StaPos,nst,spos)

C     Start locate site position on model
      DO is = 1,nst
        IF (spos(is).EQ.YDis(1)) THEN
          pnode(is) = 1
        ENDIF
        DO iy = 2,Ny+1
          IF (spos(is).EQ.YDis(iy)) THEN
            pnode(is) = iy
            GOTO 600
          ELSE
            IF ((spos(is).GT.YDis(iy-1)).AND.
     >          (spos(is).LT.YDis(iy))) THEN

              oldny     = Ny
              Ny        = Ny+1
              pnode(is) = iy

              IF (Ny.GT.NY0MX) THEN
                WRITE(6,1200)
                WRITE(6,1210) 
                WRITE(6,1220) 
                STOP
              ENDIF

              CALL CopyVectorR8(1,oldny,Dy,1,oldny,dytmp)
              Dy(iy-1) = spos(is)-YDis(iy-1)
              Dy(iy)   = YDis(iy)-spos(is)
              CALL CopyVectorR8(iy,Ny-1,dytmp,iy+1,Ny,Dy)
    
              CALL CopyMatrixR8(1,Nzb,1,oldny,NZ0MX,NY0MX,IniRho,
     >                          1,Nzb,1,oldny,NZ0MX,NY0MX,rhotmp)
              CALL CopyMatrixR8(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,IniRho)
              CALL CopyMatrixR8(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy,iy,NZ0MX,NY0MX,IniRho)
              CALL CopyMatrixR8(1,Nzb,iy,Ny-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy+1,Ny,NZ0MX,NY0MX,IniRho)

              CALL CopyMatrixR8(1,Nzb,1,oldny,NZ0MX,NY0MX,PriRho,
     >                          1,Nzb,1,oldny,NZ0MX,NY0MX,rhotmp)
              CALL CopyMatrixR8(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,PriRho)
              CALL CopyMatrixR8(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy,iy,NZ0MX,NY0MX,PriRho)
              CALL CopyMatrixR8(1,Nzb,iy,Ny-1,NZ0MX,NY0MX,rhotmp,
     >                          1,Nzb,iy+1,Ny,NZ0MX,NY0MX,PriRho)

              CALL CopyMatrixI4(1,Nzb,1,oldny,NZ0MX,NY0MX,ModGrd,
     >                          1,Nzb,1,oldny,NZ0MX,NY0MX,mgrdtmp)
              CALL CopyMatrixI4(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,mgrdtmp,
     >                          1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,ModGrd)
              CALL CopyMatrixI4(1,Nzb,iy-1,iy-1,NZ0MX,NY0MX,mgrdtmp,
     >                          1,Nzb,iy,iy,NZ0MX,NY0MX,ModGrd)
              CALL CopyMatrixI4(1,Nzb,iy,Ny-1,NZ0MX,NY0MX,mgrdtmp,
     >                          1,Nzb,iy+1,Ny,NZ0MX,NY0MX,ModGrd)

              DO iyy = 2,Ny+1
                YDis(iyy) = YDis(iyy-1) + Dy(iyy-1)
              ENDDO
              GOTO 600
            ENDIF
          ENDIF
        ENDDO
600     CONTINUE
      ENDDO ! is

      CALL DistanceBetweenBlocks(Ny,Dy,Cy)
      CALL CheckStatus(Nzb,Ny,ModGrd,DFStatus)

C     Check station location      
      DO im = 1,NMode
        DO is = 1,NSta(im) 
          DO js = 1,nst
            IF (StaPos(im,is).EQ.spos(js)) THEN
              StaNod(im,is) = pnode(js)
              GOTO 700
            ENDIF
          ENDDO ! js
700       CONTINUE
        ENDDO ! is
      ENDDO ! im

      IF (Ny.GT.Ny0) THEN
        WRITE(6,*) 
        WRITE(6,1100) Ny-Ny0
        WRITE(6,1100) Ny,Ny0
        WRITE(99,1100) Ny-Ny0
        WRITE(99,1100) Ny,Ny0
      ENDIF


1000  FORMAT('???? Move <STARTING_POS> of mode ',i1,
     >       ' to boundary between blocks ????') 
1010  FORMAT('????  from ',f12.2,' to new location ',f12.2,' ????')

1100  FORMAT('????    Add ',i3,' NY to the STARTING MODEL ????')
1110  FORMAT('????    NEW NY = ',i3,'   OLD NY = ',i3,'  ????')

1200  FORMAT('!!! ATTENTION, ERROR IN LOCATING MODEL POSITION !!!')
1210  FORMAT('!!!   NEW <NY> IS EXCEED NY0MX in parameter.h !!!')
1220  FORMAT('!!!      Please, correct and restart  !!!')


      RETURN
      END ! SUBROUTINE LocateModelPosition

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

C     Determine total no. of distinct stations of every modes
C     and sort them out in ascending order

      SUBROUTINE FindDistinctStation(NMode,NSta,StaPos,nst,spos)
      INCLUDE 'parameter.h'

C        IO Variables
      INTEGER NMode,NSta(*),nst
      REAL*8  StaPos(NMODMX,NSTAMX),spos(*)

C       Local Variables
      INTEGER im,is,js,repeatsite
      REAL*8  spostmp(NMODMX*NSTAMX),cpos

      im  = 1
      nst = NSta(im)
      DO js = 1,nst
        spos(js) = StaPos(im,js)
      ENDDO

      DO im = 2,NMode
        DO is = 1,NSta(im)
          repeatsite = 0
          cpos = StaPos(im,is)
          DO js = 1,nst
            IF (cpos.EQ.spos(js)) THEN
               repeatsite = 1
               GOTO 500
            ENDIF
          ENDDO ! js
          IF (repeatsite.EQ.0) THEN
            IF (cpos.LT.spos(1)) THEN
              CALL CopyVectorR8(1,nst,spos,1,nst,spostmp)
              spos(1) = cpos
              CALL CopyVectorR8(1,nst,spostmp,2,nst+1,spos)
              nst = nst + 1
            ENDIF
            DO js = 1,nst-1
              IF ((cpos.GT.spos(js)).AND.(cpos.LT.spos(js+1))) THEN
                CALL CopyVectorR8(1,nst,spos,1,nst,spostmp)
                spos(js+1) = cpos
                CALL CopyVectorR8(js+1,nst,spostmp,js+2,nst+1,spos)
                nst = nst + 1
              ENDIF
            ENDDO ! js
            IF (cpos.GT.spos(nst)) THEN
              spos(nst+1) = cpos
              nst = nst + 1
            ENDIF
          ENDIF

          IF (nst.GT.NMODMX*NSTAMX) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN INPUT FILE !!!'
            WRITE(6,*) '!!! NUMBER OF DISTINCT STATIONS IS MORE THAN',
     >                 'NMODMX*NSTAMX !!!' 
            WRITE(6,*) '!!! Please, correct and restart !!!'
            STOP
          ENDIF 

500       CONTINUE
        ENDDO ! is
      ENDDO ! im

      RETURN
      END ! FindDistinctStation

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
    
      SUBROUTINE CompInterPol(NMode,NRes,NPer,NSta,NNT,
     >           DatInx,SenInx,Period,SknDepth,StaPos,
     >           DatErr,CdBB,Tau,Work)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER LWork
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  CdBB(NN0MX,LL0MX),Tau(NN0MX),Work(NN0MX)
      REAL*8  SknDepth(NMODMX,NPERMX)
      REAL*8  StaPos(NMODMX,NSTAMX)
     
      INTEGER ndg,ndd,np1,np2,np3,np4  
      INTEGER idx1,idx2,ido,im,ir,ip,is
      INTEGER im1,ir1,ip1,is1,im2,ir2,ip2,is2,info
      INTEGER io(8),jo(8),npp,nss,ii
      REAL*8  ds(8),dp(8),wt1,wt2,wt(8),wtt,wf,sumds,sumdp

      INTEGER  IndexData
      EXTERNAL IndexData

      ndg = NNT(2)
      ndd = NNT(3)
      np1 = NMODMX
      np2 = NRESMX
      np3 = NPERMX
      np4 = NSTAMX

      CALL ConstantMatrixR8(CdBB,NN0MX,LL0MX,ndg,ndd,D0)

C     interpolate in periods (for now)
      idx1 = 0
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im)
              IF (DatInx(im,ir,ip,is).EQ.1) THEN
                idx1 = idx1 + 1
                IF (SenInx(im,ir,ip,is).EQ.1) THEN
                  idx2 = IndexData(SenInx,np1,np2,np3,np4,
     >                             im,ir,ip,is,NRes,NPer,NSta)
                  CdBB(idx1,idx2) = D1
                ENDIF ! IF (SenInx(im,ir,ip,is).EQ.1) 

                IF (SenInx(im,ir,ip,is).EQ.0) THEN
                  npp = NPer(im)
                  nss = NSta(im)
                  CALL FindNearest(im,ir,ip,is,npp,nss,
     >                 DatInx,SenInx,SknDepth,Period,StaPos,
     >                 np1,np2,np3,np4,io,jo,ds,dp)
                  wtt = D0
                  sumds = D0
                  sumdp = D0
                  DO ii = 1,8
                    IF (io(ii).NE.0) THEN
                      sumds  = sumds + ds(ii)  
                      sumdp  = sumdp + dp(ii)  
                    ENDIF
                  ENDDO

                  DO ii = 1,8
                    wt(ii)   = D0
                    IF (io(ii).NE.0) THEN
                      wt1    = DEXP(-ds(ii)/sumds)
                      wt2    = DEXP(-dp(ii)/sumdp)
                      wt(ii) = wt1*wt2
                      wtt    = wtt + wt(ii)
                    ENDIF
                  ENDDO

                  DO ii = 1,8
                    IF (io(ii).NE.0) THEN
                      wf = wt(ii)/wtt
                      IF (wf.GE.1E-02) THEN
                        ido = IndexData(SenInx,np1,np2,np3,np4,
     >                        im,ir,io(ii),jo(ii),NRes,NPer,NSta)
                        CdBB(idx1,ido) = wf
                      ENDIF
                    ENDIF
                  ENDDO

                ENDIF ! IF (SenInx(im,ir,ip,is).EQ.0) 
              ENDIF ! IF (DatInx(im,ir,ip,is).EQ.1) 
            ENDDO ! is
          ENDDO ! ip
        ENDDO ! ir
      ENDDO ! im


C     Multiplying with the error
      idx1 = 0
      DO im1 = 1,NMode
        DO ir1 = 1,NRes(im1)
          DO ip1 = 1,NPer(im1)
            DO is1 = 1,NSta(im1)
              IF (DatInx(im1,ir1,ip1,is1).EQ.1) THEN
                idx1 = idx1 + 1
                idx2 = 1
                DO im2 = 1,NMode
                  DO ir2 = 1,NRes(im2)
                    DO ip2 = 1,NPer(im2)
                      DO is2 = 1,NSta(im2)
                        IF ((DatInx(im2,ir2,ip2,is2).EQ.1).AND.
     >                      (SenInx(im2,ir2,ip2,is2).EQ.1)) THEN
                          CdBB(idx1,idx2) = 
     >                    CdBB(idx1,idx2)/DatErr(im1,ir1,ip1,is1)
                          idx2 = idx2 + 1
                        ENDIF
                      ENDDO ! is2
                    ENDDO ! ip2
                  ENDDO ! ir2
                ENDDO ! im2
              ENDIF
            ENDDO ! is1
          ENDDO ! ip1
        ENDDO ! ir1
      ENDDO ! im1

C     Compute QR Factorization
      LWork = ndd
      CALL DGEQRF(ndg,ndd,CdBB,NN0MX,Tau,Work,LWork,info)
      IF (info.LT.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR QR FACTORIZATION !!!'
        WRITE(6,*)'!!! Please, correct your input file and restart!!!'
        STOP
      ENDIF


      RETURN
      END ! CompInterPol()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
