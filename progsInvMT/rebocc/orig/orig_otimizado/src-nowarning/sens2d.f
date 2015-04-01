
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SubSens2d(Nit,NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           Period,StaNod,NNT,DFStatus,SSPara,DatInx,SenInx,
     >           dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >           CRho,ETOL,NtMx,FWD_Only,
     >           CmHtIndx,CmHt,FOUT)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nzb,Nz,Ny,Nit
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*)
      REAL*8  CRho(NZ0MX,NY0MX)
      REAL*8  Period(NMODMX,NPERMX)
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*),NNT(*)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER DFStatus(*),NtMx,FWD_Only
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX),SSPara(NMODMX,NSTAMX),ETOL
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)
      REAL*8  FOUT(NMODMX,NRESMX,NPERMX,NSTAMX)

      INTEGER ipiv(MMIMX)
      COMPLEX*16 AII(NZ3MX,MMIMX)
      COMPLEX*16 HXI(MMIMX),HXB(MMBMX)
      COMPLEX*16 EXI(MMIMX),EXB(MMBMX)

      INTEGER im,ir,ip,is
      INTEGER np1,np2,igx
      INTEGER flagsens
      REAL*8  per,App(NSTAMX),Phs(NSTAMX)
      REAL*8  AA(MMIMX,4),BD(MMBMX)
      COMPLEX*16 Zyx(NSTAMX),Hxs(NSTAMX),Eys(NSTAMX),Tipper(NSTAMX)
      COMPLEX*16 Zxy(NSTAMX),Exs(NSTAMX),Hys(NSTAMX),Hzs(NSTAMX)

      INTEGER DoSens,fwdcond
     

      np1 = NMODMX
      np2 = NRESMX

      igx = 0
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            CmHtIndx(im,ir,ip) = igx
            DO is = 1,NSta(im)
              IF (SenInx(im,ir,ip,is).EQ.1) THEN
                igx = igx + 1 
              ENDIF
            ENDDO ! is
          ENDDO ! ip
        ENDDO ! ir
      ENDDO ! im

      DO im = 1,NMode
        IF (ModTyp(im).EQ.1)
     >     CALL SetupA_TM(Nzb,Ny,Dzb,Dy,Czb,Cy,CRho,AA,BD)
        IF ((ModTyp(im).EQ.2).OR.(ModTyp(im).EQ.3))
     >     CALL SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,CRho,AA,BD)
    
        DO ip = 1,NPer(im)
          per = Period(im,ip)
          flagsens = DoSens(im,ip,NRes,NSta,SenInx)
          IF (FWD_Only.EQ.1) flagsens = 0

          IF (flagsens.EQ.1) THEN
            IF (ModTyp(im).EQ.1) THEN
              CALL Fwd2DTM_LU(per,Nzb,Ny,Dzb,Dy,CRho,AA,BD,
     >        im,NSta,StaNod,AII,ipiv,HXI,HXB,App,Phs,Zyx,Hxs,Eys)

              CALL Sens2DTM(im,ip,per,Nzb,Ny,Dzb,Dy,CRho,
     >             NRes,ResTyp,NSta,StaNod,SenInx,DatInx,
     >             DFStatus,dATM,
     >             AII,ipiv,HXI,HXB,App,Phs,Zyx,CmHtIndx,CmHt)
            ENDIF
            IF (ModTyp(im).EQ.2) THEN
              CALL Fwd2DTE_LU(per,Nza,Nz,Ny,Dz,Dy,CRho,AA,BD,
     >        im,NSta,StaNod,AII,ipiv,EXI,EXB,App,Phs,Tipper,
     >        Zxy,Exs,Hys,Hzs)

              CALL Sens2DTE(im,ip,per,Nza,Nz,Ny,Dz,Dy,CRho,
     >             NRes,ResTyp,NSta,StaNod,SenInx,DatInx,
     >             DFStatus,dATE,AII,ipiv,
     >             EXI,EXB,App,Phs,Zxy,Exs,Hys,CmHtIndx,CmHt)
            ENDIF
            IF (ModTyp(im).EQ.3) THEN
              CALL Fwd2DTE_LU(per,Nza,Nz,Ny,Dz,Dy,CRho,AA,BD,
     >        im,NSta,StaNod,AII,ipiv,EXI,EXB,App,Phs,Tipper,
     >        Zxy,Exs,Hys,Hzs)

              CALL Sens2DTipper(im,ip,per,Nza,Nz,Ny,Dz,Dy,CRho,
     >             NRes,ResTyp,NSta,StaNod,SenInx,DatInx,
     >             DFStatus,dATE,
     >             AII,ipiv,EXI,EXB,Hys,Hzs,CmHtIndx,CmHt)
            ENDIF

          ELSE
            IF ((ModTyp(im).EQ.1).AND.(Nit.EQ.1)) THEN
              CALL Fwd2DTM_PCG(per,Nzb,Ny,Dzb,Dy,CRho,AA,BD,
     >             im,NSta,StaNod,ETOL,NtMx,fwdcond,HXI,HXB,
     >             App,Phs,Zyx,Hxs,Eys)
            ENDIF
            IF ((ModTyp(im).EQ.2).AND.(Nit.EQ.1)) THEN
              CALL Fwd2DTE_PCG(per,Nza,Nz,Ny,Dz,Dy,CRho,AA,BD,
     >             im,NSta,StaNod,ETOL,NtMx,fwdcond,
     >             EXI,EXB,App,Phs,Tipper,Zxy,Exs,Hys,Hzs)
            ENDIF 
            IF ((ModTyp(im).EQ.3).AND.(Nit.EQ.1)) THEN
              CALL Fwd2DTE_PCG(per,Nza,Nz,Ny,Dz,Dy,CRho,AA,BD,
     >             im,NSta,StaNod,ETOL,NtMx,fwdcond,
     >             EXI,EXB,App,Phs,Tipper,Zxy,Exs,Hys,Hzs)
            ENDIF 
          ENDIF ! flagsens

          IF (Nit.EQ.1) THEN
            DO ir = 1,NRes(im)
              IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
                DO is = 1,NSta(im) 
                  FOUT(im,ir,ip,is) = App(is) + SSPara(im,is)
                ENDDO
              ENDIF ! IF 1
              IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN 
                DO is = 1,NSta(im) 
                  FOUT(im,ir,ip,is) = Phs(is)
                ENDDO
              ENDIF ! IF 2
              IF (ResTyp(im,ir).EQ.5) THEN
                DO is = 1,NSta(im) 
                  FOUT(im,ir,ip,is) = DREAL(Tipper(is))
                ENDDO
              ENDIF ! IF 5
              IF (ResTyp(im,ir).EQ.6) THEN
                DO is = 1,NSta(im) 
                  FOUT(im,ir,ip,is) = IMAG(Tipper(is))
                ENDDO
              ENDIF ! IF 5
            ENDDO
          ENDIF ! IF Nit
       
        ENDDO ! ip
      ENDDO ! im


      RETURN
      END ! SubSens2d()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      INTEGER FUNCTION DoSens(im,ip,NRes,NSta,SenInx)
      INCLUDE 'parameter.h'
      INTEGER im,ip,NRes(*),NSta(*)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)

      INTEGER ir,is

      DoSens = 0
      DO ir = 1,NRes(im)
        DO is = 1,NSta(im)
          IF (SenInx(im,ir,ip,is).EQ.1) THEN
            DoSens = 1
            GOTO 100
          ENDIF
        ENDDO ! is
      ENDDO ! ir
100   CONTINUE

      RETURN
      END ! DoSens

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 
      SUBROUTINE Sens2DTM(im,ip,per,Nzb,Ny,Dzb,Dy,CRho,
     >           NRes,ResTyp,NSta,StaNod,SenInx,DatInx,DFStatus,
     >           dATM,AII,ipiv,HXI,HXB,App,Phs,Zyx,CmHtIndx,CmHt)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny,im,ip,NSta(*),NRes(*)
      REAL*8 per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER DFStatus(*),ipiv(*)
      COMPLEX*16 AII(NZ3MX,MMIMX),HXI(*),HXB(*),Zyx(*)
      REAL*8  App(*),Phs(*)
      REAL*8  dATM(MM0MX,3)
      INTEGER StaNod(NMODMX,NSTAMX)
    
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX),EA(MM0MX),EB(MM0MX)
      COMPLEX*16 dLR(MM0MX),dLL(MM0MX),dUL(MM0MX),dUR(MM0MX)
      COMPLEX*16 dQdR1(4),dQdR2(MM0MX)

      INTEGER jj,is,ir,igx,iss(NRESMX)
      INTEGER docompute,js,jss,jsr,jsl,iz,kl,ku,info,mmi,mmt
      REAL*8  omega,omue,dz1,dz2,dzz,aptm,phtm,cos2,ryy,rzr,rzl,dyy
      COMPLEX*16 ai(MMIMX),zztm

      INTEGER  DoSens

      omega = (D2*Pi)/per
      omue   = omega*Mue

      mmt = Ny*Nzb
      mmi = (Ny-1)*(Nzb-1)

      dz1 = Dzb(1)
      dz2 = Dzb(2)
      dzz = dz1 + dz2
      CALL CompZI1_TM(Nzb,Ny,DFStatus,HXI,HXB,dATM,dLR,dLL,dUL,dUR)

      kl = Nzb-1
      ku = Nzb-1

      DO ir = 1,NRes(im)
        iss(ir) = 0
      ENDDO ! ir
      DO is = 1,NSta(im)
        docompute = DoSens(im,ip,NRes,NSta,SenInx)

        IF (docompute.EQ.1) THEN
          CALL ConstantVectorC16(ai,mmi,D0)

          aptm = D10**App(is)         
          phtm = Phs(is)
          zztm = Zyx(is)
          cos2 = (DCOS(phtm)**D2)/DREAL(Zyx(is))
          iz   = 2
          js   = StaNod(im,is)
          dyy  = Dy(js) + Dy(js-1)
          ryy  = (CRho(1,js)*Dy(js) + CRho(1,js-1)*Dy(js-1))/dyy
          rzr  = (CRho(1,js)*dz1    + Crho(2,js)*dz2)/dzz
          rzl  = (CRho(1,js-1)*dz1  + Crho(2,js-1)*dz2)/dzz
         
C         one node beneath the surface node ...
          jss  = (js-2)*(Nzb-1) + (iz-1)
          jsr  = ((js+1)-2)*(Nzb-1) + (iz-1)
          jsl  = ((js-1)-2)*(Nzb-1) + (iz-1)
          ai(jss)   = ryy/dz1 + DCMPLX(D0,omue*dz1/D8) 
     >              - rzr*dz1/(D4*Dy(js)*dyy) 
     >              - rzl*dz1/(D4*Dy(js-1)*dyy)
          ai(jsr) = rzr*dz1/(D4*Dy(js)*dyy)
          ai(jsl) = rzl*dz1/(D4*Dy(js-1)*dyy)
          IF ((js.LT.1).or.(js.GT.Ny)) THEN 
            WRITE(6,*) '!!! ATTENTION, ERROR ai TM CASE !!! '
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         solve A^{T}*ui = ai wher ai indicate location node
          CALL ZGBTRS('T',mmi,kl,ku,1,AII,NZ3MX,ipiv,ai,MMIMX,info)
          IF (INFO.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR SOLVING SENS2D_TM ',info
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         dQdR0  =  (dai/dR)*XI+(dab/dR)*XB
C                   +------ dQdR1 --------+
C                   + -ui^{T}*[(dAi/dR)*XI+(dAb/dR)*Xb + (AIB)*(dXB/dR)]
C                             +---------- ZI1 -------+   +---- ZI2 ----+
C                    +--------------- dQdR2 ---------------------------+
C         assumes ZI2 = 0;

          CALL CompdQdR1_TM(im,is,Nzb,Ny,Dzb,Dy,StaNod,HXI,dQdR1)
          CALL CompdQdR2_TM(Nzb,Ny,dLR,dLL,dUL,dUR,ai,dQdR2)
          CALL CompdQdR0_TM(im,is,per,Nzb,Ny,StaNod,CRho,
     >         aptm,phtm,zztm,cos2,DFStatus,dQdR1,dQdR2,EA,EB)

          DO ir = 1,NRes(im)
            igx = CmHtIndx(im,ir,ip)
            IF (SenInx(im,ir,ip,is).EQ.1) THEN
              iss(ir) = iss(ir) + 1
              DO jj = 1,mmt
                IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
                  CmHt(jj,igx+iss(ir)) = EA(jj)
                ENDIF
                IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN 
                  CmHt(jj,igx+iss(ir)) = EB(jj)
                ENDIF
              ENDDO
            ENDIF
          ENDDO ! ir

        ENDIF ! docompute
      ENDDO ! is

      RETURN
      END ! Sens2DTM

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ZI1 = (dAii/dR)*Xi + (dAib/dR)*Xb

      SUBROUTINE CompZI1_TM(Nzb,Ny,DFStatus,HXI,HXB,dATM,
     >           dLR,dLL,dUL,dUR)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
  
      INTEGER Ny,Nzb,DFStatus(*)
      REAL*8  dATM(MM0MX,3)
      COMPLEX*16 HXI(*),HXB(*)
      COMPLEX*16  dLR(MM0MX),dLL(MM0MX),dUL(MM0MX),dUR(MM0MX)
      
      INTEGER iz,iy,ilr,ill,iul,iur,iz1,iz2,iz3,jj,mm
      REAL*8  dA1,dA2,dA3

      iz1 = Nzb+1
      iz2 = iz1 + (Ny-1)
      iz3 = iz2 + (Ny-1)

      mm = Ny*Nzb
      CALL ConstantVectorC16(dLR,mm,D0)
      CALL ConstantVectorC16(dLL,mm,D0)
      CALL ConstantVectorC16(dUL,mm,D0)
      CALL ConstantVectorC16(dUR,mm,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj  = jj + 1

          IF (DFStatus(jj).EQ.0) THEN
            GOTO 100
          ENDIF

          iul = (iy-2)*(Nzb-1) + (iz-1)
          ill = (iy-2)*(Nzb-1) +  iz
          iur = (iy-1)*(Nzb-1) + (iz-1)
          ilr = (iy-1)*(Nzb-1) +  iz
          dA1 = dATM(jj,1)
          dA2 = dATM(jj,2)
          dA3 = dATM(jj,3)
        
          IF (iy.EQ.1) THEN
            IF (iz.EQ.1) THEN
              dLR(jj) = dA1*HXI(ilr) + dA2*HXB(2) + dA3*HXB(iz1+1)
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dUR(jj) = dA1*HXI(iur) + dA2*HXB(Nzb) + dA3*HXB(iz2+1)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dUR(jj) = dA1*HXI(iur) + dA2*HXB(iz)   + dA3*HXI(ilr)  
              dLR(jj) = dA1*HXI(ilr) + dA2*HXB(iz+1) + dA3*HXI(iur)  
            ENDIF
          ENDIF

          IF (iy.EQ.Ny) THEN
            IF (iz.EQ.1) THEN
              dLL(jj) = dA1*HXI(ill)+dA2*HXB(iz3+iz+1)+dA3*HXB(iz2)
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dUL(jj) = dA1*HXI(iul)+dA2*HXB(iz3+iz) + dA3*HXB(iz3)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dUL(jj) = dA1*HXI(iul)+dA2*HXB(iz3+iz)+dA3*HXI(ill)
              dLL(jj) = dA1*HXI(ill)+dA2*HXB(iz3+iz+1)+dA3*HXI(iul)
            ENDIF
          ENDIF

          IF ((iy.GT.1).AND.(iy.LT.Ny)) THEN
            IF (iz.EQ.1) THEN
              dLR(jj) = dA1*HXI(ilr) + dA2*HXI(ill)+dA3*HXB(iz1+iy)  
              dLL(jj) = dA1*HXI(ill) + dA2*HXI(ilr)+dA3*HXB(iz1+iy-1)
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dUR(jj) = dA1*HXI(iur) + dA2*HXI(iul)+dA3*HXB(iz2+iy)
              dUL(jj) = dA1*HXI(iul) + dA2*HXI(iur)+dA3*HXB(iz2+iy-1)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dUL(jj) = dA1*HXI(iul) + dA2*HXI(iur) + dA3*HXI(ill)
              dLL(jj) = dA1*HXI(ill) + dA2*HXI(ilr) + dA3*HXI(iul)
              dUR(jj) = dA1*HXI(iur) + dA2*HXI(iul) + dA3*HXI(ilr)  
              dLR(jj) = dA1*HXI(ilr) + dA2*HXI(ill) + dA3*HXI(iur)  
            ENDIF
          ENDIF

100       CONTINUE
        ENDDO ! iz
      ENDDO ! iy

      RETURN
      END ! CompZI1_TM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Compute dQdR2 TM case (no boundary derivative)
C     zi = (dAii/dR)*Xi + (dAib/dR)*Xb + Aib*(dXb/dR)
C     zi = ZI1 - ZI2   where we assume ZI2 = 0.
C     dQdR2 = -ui^{T}*(ZI1 - ZI2); where   Aii*ui = ai

      SUBROUTINE CompdQdR2_TM(Nzb,Ny,dLR,dLL,dUL,dUR,ui,dQdR2)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb
      COMPLEX*16 dUR(*),dUL(*),dLR(*),dLL(*),dQdR2(*),ui(*)

      INTEGER iur,iul,ilr,ill,iz,iy,jj,mmt

      mmt = Ny*Nzb
      CALL ConstantVectorC16(dQdR2,mmt,D0)
     
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = (iy-1)*Nzb + iz

          iul = (iy-2)*(Nzb-1) + (iz-1)
          ill = (iy-2)*(Nzb-1) +  iz
          iur = (iy-1)*(Nzb-1) + (iz-1)
          ilr = (iy-1)*(Nzb-1) +  iz

          IF (iy.EQ.1) THEN
            IF (iz.EQ.1) THEN
              dQdR2(jj) = -ui(ilr)*dLR(jj) 
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dQdR2(jj) = -ui(iur)*dUR(jj)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dQdR2(jj) = -ui(iur)*dUR(jj) -ui(ilr)*dLR(jj) 
            ENDIF
          ENDIF

          IF (iy.EQ.Ny) THEN
            IF (iz.EQ.1) THEN
              dQdR2(jj) = -ui(ill)*dLL(jj)
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dQdR2(jj) = -ui(iul)*dUL(jj)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dQdR2(jj) = -ui(iul)*dUL(jj) -ui(ill)*dLL(jj)
            ENDIF
          ENDIF

          IF ((iy.GT.1).AND.(iy.LT.Ny)) THEN
            IF (iz.EQ.1) THEN
              dQdR2(jj) = -ui(ill)*dLL(jj) -ui(ilr)*dLR(jj) 
            ENDIF
            IF (iz.EQ.Nzb) THEN
              dQdR2(jj) = -ui(iul)*dUL(jj) -ui(iur)*dUR(jj)
            ENDIF
            IF ((iz.GT.1).AND.(iz.LT.Nzb)) THEN
              dQdR2(jj) = -ui(iul)*dUL(jj) -ui(ill)*dLL(jj)
     >                    -ui(iur)*dUR(jj) -ui(ilr)*dLR(jj) 
            ENDIF
          ENDIF


        ENDDO ! iz
      ENDDO !iy 

   
      RETURN
      END ! CompdQdR2_TM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     dQdR1  =  (dai/dR)*XI+(dab/dR)*XB

      SUBROUTINE CompdQdR1_TM(im,is,Nzb,Ny,Dzb,Dy,StaNod,HXI,dQdR1)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      REAL*8  Dzb(*),Dy(*)
      COMPLEX*16 HXI(*),dQdR1(*)

      INTEGER iz,iy,jss,jsr,jsl,js
      REAL*8  dz1,dz2,dzz,dyy,c1,c2,c3


      dz1 = Dzb(1)
      dz2 = Dzb(2)
      dzz = dz1 + dz2
      js  = StaNod(im,is)
      dyy = Dy(js) + Dy(js-1)

C     at R(1,js)
      iy = js
      c2 = (dz1*dz1)/(D4*Dy(iy)*dyy*dzz)
      c3 = -Dy(js)/(dyy*dz1) 
      c1 = -c3 - c2
      iz = 2
      jss  = (js-2)*(Nzb-1)     + (iz-1)
      jsr  = ((js+1)-2)*(Nzb-1) + (iz-1)
      dQdR1(1) = c1*HXI(jss) + c2*HXI(jsr) + c3*D1

C     at R(1,js-1)
      iy = js-1
      c2 = (dz1*dz1)/(D4*Dy(iy)*dyy*dzz)
      c3 = -Dy(iy)/(dyy*dz1) 
      c1 = -c3 - c2
      iz = 2
      jss  = (js-2)*(Nzb-1)     + (iz-1)
      jsl  = ((js-1)-2)*(Nzb-1) + (iz-1)
      dQdR1(2) = c1*HXI(jss) + c2*HXI(jsl) + c3*D1

C     at R(2,js)
      iy = js
      c1 = -(dz1*dz2)/(D4*Dy(iy)*dyy*dzz)
      c2 = -c1
      iz = 2
      jss  = (js-2)*(Nzb-1)     + (iz-1)
      jsr  = ((js+1)-2)*(Nzb-1) + (iz-1)
      dQdR1(3) = c1*HXI(jss) + c2*HXI(jsr)

C     at R(2,js-1)
      iy = js-1
      c1 = -(dz1*dz2)/(D4*Dy(iy)*dyy*dzz)
      c2 = -c1
      iz = 2
      jss  = (js-2)*(Nzb-1)     + (iz-1)
      jsl  = ((js-1)-2)*(Nzb-1) + (iz-1)
      dQdR1(4) = c1*HXI(jss) + c2*HXI(jsl)

      RETURN
      END ! CompdQdR1_TM()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Compute only in the sensitive zone

      SUBROUTINE CompdQdR0_TM(im,is,per,Nzb,Ny,StaNod,CRho,
     >                 aptm,phtm,zztm,cos2,DFStatus,dQdR1,dQdR2,EA,EB)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER DFStatus(MM0MX)
      REAL*8  CRho(NZ0MX,NY0MX)
      REAL*8  aptm,phtm,cos2,per,EA(*),EB(*)
      COMPLEX*16 dQdR1(*),dQdR2(*),zztm

      INTEGER jj,mmt,js,iy,iz
      REAL*8  omue,omega,dRdR,dPdR
      COMPLEX*16 dEdR(MM0MX),dZdR

      omega = (D2*Pi)/per
      omue  = omega*Mue
      mmt = Ny*Nzb
      js  = StaNod(im,is)

      CALL CopyVectorC16(1,mmt,dQdR2,1,mmt,dEdR)
      iz = 1
      iy = js
      jj = (iy-1)*Nzb + iz
      dEdR(jj) = dEdR(jj) + dQdR1(1)

      iz = 1
      iy = js-1
      jj = (iy-1)*Nzb + iz
      dEdR(jj) = dEdR(jj) + dQdR1(2)

      iz = 2
      iy = js
      jj = (iy-1)*Nzb + iz
      dEdR(jj) = dEdR(jj) + dQdR1(3)

      iz = 2
      iy = js-1
      jj = (iy-1)*Nzb + iz
      dEdR(jj) = dEdR(jj) + dQdR1(4)

C     TM : zyx = -Ey/Hx 
C     dZyx/dR  = -[(Hx*dEy-Ey*dHx)/(Hx^2)] = -dEy

      CALL ConstantVectorR8(EA,mmt,D0)
      CALL ConstantVectorR8(EB,mmt,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = jj+1

          dZdR = -dEdR(jj)
          dRdR = (D2/omue)*DREAL(CONJG(zztm)*dZdR)
          dPdR = cos2*(IMAG(dZdR)-DTAN(phtm)*DREAL(dZdR))

          EA(jj) = (CRho(iz,iy)*dRdR/aptm)*DFStatus(jj)
          EB(jj) = (CRho(iz,iy)*dPdR*DLOG(D10))*DFStatus(jj)
        ENDDO
      ENDDO ! iy


      RETURN
      END ! CompdQdR0_TM()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Sens2DTE(im,ip,per,Nza,Nz,Ny,Dz,Dy,CRho,
     >           NRes,ResTyp,NSta,StaNod,SenInx,DatInx,DFStatus,
     >           dATE,AII,ipiv,EXI,EXB,App,Phs,Zxy,Exs,Hys,
     >           CmHtIndx,CmHt)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz,Ny,im,ip,NSta(*),NRes(*)
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER ipiv(MMIMX)
      INTEGER DFStatus(MM0MX)
      COMPLEX*16 AII(NZ3MX,MMIMX)
      COMPLEX*16 EXI(MMIMX),EXB(MMBMX)
      REAL*8  App(NSTAMX),Phs(NSTAMX)
      COMPLEX*16 Zxy(NSTAMX),Exs(NSTAMX),Hys(NSTAMX)
      REAL*8  dATE(MM0MX,4)
      INTEGER StaNod(NMODMX,NSTAMX)
    
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX),EA(MM0MX),EB(MM0MX)
      COMPLEX*16 dLR(MM0MX),dLL(MM0MX),dUL(MM0MX),dUR(MM0MX)
      COMPLEX*16 dQdR1(2),dQdR2(MM0MX)

      
      INTEGER is,ir,Nzb,jj,j1s,j1r,j1l,igx,iss(NRESMX)
      INTEGER docompute,js,jss,jsr,jsl,iz,kl,ku,info,mmi,mmt
      REAL*8  omega,omue,dz1,dz2,dzz,apte,phte,cos2,dyy,sav
      COMPLEX*16 bi(MMIMX),zzte,dyz,zhte

      INTEGER  DoSens

      omega = (D2*Pi)/per
      omue   = omega*Mue

      Nzb = Nz-Nza
      mmt = Ny*Nzb
      mmi = (Ny-1)*(Nz-1)

      dz1 = Dz(Nza+1)
      dz2 = Dz(Nza+2)
      dzz = dz1 + dz2
      CALL CompZI1_TE(Nza,Nz,Ny,CRho,DFStatus,EXI,dATE,
     >     dLR,dLL,dUL,dUR)

      kl = Nz-1
      ku = Nz-1

      DO ir = 1,NRes(im)
        iss(ir) = 0
      ENDDO ! ir
      DO is = 1,NSta(im)
        docompute = DoSens(im,ip,NRes,NSta,SenInx)

        IF (docompute.EQ.1) THEN
          CALL ConstantVectorC16(bi,mmi,D0)

          apte = D10**App(is)         
c         phte = -D1*Phs(is)
          phte = Phs(is)
          zzte = Zxy(is)
          cos2 = (DCOS(phte)**D2)/DREAL(zzte)
          zhte = -zzte/Hys(is)

          js   = StaNod(im,is)
          dyy  = Dy(js) + Dy(js-1)
          sav  = (Dy(js)/CRho(1,js) + Dy(js-1)/CRho(1,js-1))/dyy
          dyz  = DCMPLX(D0,dz1/(dyy*omue))

          iz   = Nza+1
C         one node beneath the surface node ...
          jss  = (js-2)*(Nz-1)     + (iz-1)
          jsr  = ((js+1)-2)*(Nz-1) + (iz-1)
          jsl  = ((js-1)-2)*(Nz-1) + (iz-1)
          j1s  = (js-2)*(Nz-1)     + (iz)
          j1r  = ((js+1)-2)*(Nz-1) + (iz)
          j1l  = ((js-1)-2)*(Nz-1) + (iz)

          bi(jss) = D1/Hys(is) + 
     >              zhte*(DCMPLX(D0,D1/(omue*dz1)) + 
     >                    (D3/D8)*sav*dz1 + 
     >                    (D3/D4)*dyz*(D1/Dy(js) + D1/Dy(js-1)))
          bi(jsr) = zhte*(-D3/D4)*dyz*(D1/Dy(js))
          bi(jsl) = zhte*(-D3/D4)*dyz*(D1/Dy(js-1))
          bi(j1s) = zhte*(DCMPLX(D0,-D1/(omue*dz1)) +
     >                    (D1/D8)*sav*dz1 + 
     >                    (D1/D4)*dyz*(D1/Dy(js) + D1/Dy(js-1)))
          bi(j1r) = zhte*(-D1/D4)*dyz*(D1/Dy(js))
          bi(j1l) = zhte*(-D1/D4)*dyz*(D1/Dy(js-1))

          IF ((js.LT.1).or.(js.GT.Ny)) THEN 
            WRITE(6,*) '!!! ATTENTION, ERROR bi TE CASE !!! '
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         solve A^{T}*ui = bi wher bi indicate location node
          CALL ZGBTRS('T',mmi,kl,ku,1,AII,NZ3MX,ipiv,bi,MMIMX,info)
          IF (INFO.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR SOLVING SENS2D_TE ',info
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         dQdR0  =  (dai/dR)*XI+(dab/dR)*XB
C                   +------ dQdR1 --------+
C                   + -ui^{T}*[(dAi/dR)*XI+(dAb/dR)*Xb + (AIB)*(dXB/dR)]
C                             +---------- ZI1 -------+   +---- ZI2 ----+
C                    +--------------- dQdR2 ---------------------------+
C         assumes ZI2 = 0;

          CALL CompdQdR1_TE(im,is,Nza,Nz,Ny,Dz,Dy,CRho,Zxy,Hys,
     >                      StaNod,EXI,dQdR1)
          CALL CompdQdR2_TE(per,Nza,Nz,Ny,dLR,dLL,dUL,dUR,bi,dQdR2)
          CALL CompdQdR0_TE(im,is,per,Nzb,Ny,StaNod,CRho,
     >         apte,phte,zzte,cos2,DFStatus,dQdR1,dQdR2,EA,EB)

          DO ir = 1,NRes(im)
            igx = CmHtIndx(im,ir,ip)
            IF (SenInx(im,ir,ip,is).EQ.1) THEN
              iss(ir) = iss(ir) + 1
              DO jj = 1,mmt
                IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
                  CmHt(jj,igx+iss(ir)) = EA(jj)
                ENDIF
                IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN 
                  CmHt(jj,igx+iss(ir)) = EB(jj)
                ENDIF
              ENDDO
            ENDIF
          ENDDO ! ir
        ENDIF ! docompute
      ENDDO ! is

      RETURN
      END ! Sens2DTE

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ZI1 = (dAii/dR)*Xi + (dAib/dR)*Xb

      SUBROUTINE CompZI1_TE(Nza,Nz,Ny,CRho,DFStatus,EXI,dATE,
     >           dLR,dLL,dUL,dUR)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
  
      INTEGER Ny,Nza,Nz,DFStatus(*)
      REAL*8  dATE(MM0MX,4),CRho(NZ0MX,NY0MX)
      COMPLEX*16 EXI(*)
      COMPLEX*16  dLR(MM0MX),dLL(MM0MX),dUL(MM0MX),dUR(MM0MX)
      
      INTEGER iz,iy,ilr,ill,iul,iur,jj,Nzb,mm
      REAL*8  dA1,dA2,dA3,dA4,s2

      Nzb = Nz-Nza

      mm = Ny*Nzb
      CALL ConstantVectorC16(dLR,mm,D0)
      CALL ConstantVectorC16(dLL,mm,D0)
      CALL ConstantVectorC16(dUL,mm,D0)
      CALL ConstantVectorC16(dUR,mm,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = Nza+1,Nz
          jj = jj+1

          IF (DFStatus(jj).EQ.0) THEN
            GOTO 100
          ENDIF

          iul = (iy-2)*(Nz-1) + (iz-1)
          ill = (iy-2)*(Nz-1) +  iz
          iur = (iy-1)*(Nz-1) + (iz-1)
          ilr = (iy-1)*(Nz-1) +  iz
          dA1 = dATE(jj,1)
          dA2 = dATE(jj,2)
          dA3 = dATE(jj,3)
          dA4 = dATE(jj,4)
          s2  = -(D1/CRho(iz-Nza,iy))**D2
        
          IF (iy.EQ.1) THEN
            IF (iz.EQ.Nz) THEN
              dUR(jj) = dA2*EXI(iur)*s2
            ENDIF
            IF (iz.LT.Nz) THEN
              dUR(jj) = dA2*EXI(iur)*s2
              dLR(jj) = dA4*EXI(ilr)*s2
            ENDIF
          ENDIF

          IF (iy.EQ.Ny) THEN
            IF (iz.EQ.Nz) THEN
              dUL(jj) = dA1*EXI(iul)*s2
            ENDIF
            IF (iz.LT.Nz) THEN
              dUL(jj) = dA1*EXI(iul)*s2
              dLL(jj) = dA3*EXI(ill)*s2
            ENDIF
          ENDIF

          IF ((iy.GT.1).AND.(iy.LT.Ny)) THEN
            IF (iz.EQ.Nz) THEN
              dUL(jj) = dA1*EXI(iul)*s2
              dUR(jj) = dA2*EXI(iur)*s2
            ENDIF
            IF (iz.LT.Nz) THEN
              dUL(jj) = dA1*EXI(iul)*s2
              dUR(jj) = dA2*EXI(iur)*s2
              dLL(jj) = dA3*EXI(ill)*s2
              dLR(jj) = dA4*EXI(ilr)*s2
            ENDIF
          ENDIF

100       CONTINUE
        ENDDO ! iz
      ENDDO ! iy

      RETURN
      END ! CompZI1_TE()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Compute dQdR2 TE case (no boundary derivative)
C     zi = (dAii/dR)*Xi + (dAib/dR)*Xb + Aib*(dXb/dR)
C     zi = ZI1 - ZI2   where we assume ZI2 = 0.
C     dQdR2 = -vi^{T}*(ZI1 - ZI2); where   Aii*vi = bi

      SUBROUTINE CompdQdR2_TE(per,Nza,Nz,Ny,dLR,dLL,dUL,dUR,vi,dQdR2)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nza,Nz
      REAL*8  per
      COMPLEX*16 dUR(*),dUL(*),dLR(*),dLL(*),dQdR2(*),vi(*)

      INTEGER iur,iul,ilr,ill,iz,iy,jj,mmt,Nzb
      REAL*8  omega,omue
      COMPLEX*16 comue

      Nzb = Nz-Nza
      omega = (D2*Pi)/per
      omue  = omega*Mue
      comue = DCMPLX(D0,omue)

      mmt = Ny*Nzb
      CALL ConstantVectorC16(dQdR2,mmt,D0)
     
      DO iy = 1,Ny
        DO iz = Nza+1,Nz
          jj  = (iy-1)*Nzb + (iz-Nza)

          iul = (iy-2)*(Nz-1) + (iz-1)
          ill = (iy-2)*(Nz-1) +  iz
          iur = (iy-1)*(Nz-1) + (iz-1)
          ilr = (iy-1)*(Nz-1) +  iz

          IF (iy.EQ.1) THEN
            IF (iz.EQ.Nz) THEN
              dQdR2(jj) = (-vi(iur)*dUR(jj))*comue
            ENDIF
            IF (iz.LT.Nz) THEN
              dQdR2(jj) = (-vi(iur)*dUR(jj) -vi(ilr)*dLR(jj))*comue 
            ENDIF
          ENDIF

          IF (iy.EQ.Ny) THEN
            IF (iz.EQ.Nz) THEN
              dQdR2(jj) = (-vi(iul)*dUL(jj))*comue
            ENDIF
            IF (iz.LT.Nz) THEN
              dQdR2(jj) = (-vi(iul)*dUL(jj) -vi(ill)*dLL(jj))*comue
            ENDIF
          ENDIF

          IF ((iy.GT.1).AND.(iy.LT.Ny)) THEN
            IF (iz.EQ.Nz) THEN
              dQdR2(jj) = (-vi(iul)*dUL(jj) -vi(iur)*dUR(jj))*comue
            ENDIF
            IF (iz.LT.Nz) THEN
              dQdR2(jj) = (-vi(iul)*dUL(jj) -vi(ill)*dLL(jj)
     >                     -vi(iur)*dUR(jj) -vi(ilr)*dLR(jj))*comue 
            ENDIF
          ENDIF
        ENDDO ! iz
      ENDDO !iy 

      RETURN
      END ! CompdQdR2_TE()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     dQdR1 = (-zxy(js)/hys(js))*[(dai/dR)*XI + (dab/dR)*XB]


      SUBROUTINE CompdQdR1_TE(im,is,Nza,Nz,Ny,Dz,Dy,CRho,Zxy,Hys,
     >                        StaNod,EXI,dQdR1)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nza,Nz,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      REAL*8  Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 EXI(*),dQdR1(*),Zxy(*),Hys(*)

      INTEGER iz,jss,js
      REAL*8  dz1,dz2,dzz,dyy
      COMPLEX*16 s2,c00,c10


      dz1 = Dz(Nza+1)
      dz2 = Dz(Nza+2)
      dzz = dz1 + dz2
      js  = StaNod(im,is)
      dyy = Dy(js) + Dy(js-1)
      s2  = (-Zxy(is)/Hys(is))*(-D1/(CRho(1,js)**D2))

      iz  = Nza + 1
      jss = (js-2)*(Nz-1) + (iz-1)

C     at R(Nza+1,js)
      c00 = s2*(D3/D8)*dz1*Dy(js)/dyy
      c10 = s2*(D1/D8)*dz1*Dy(js)/dyy
      dQdR1(1) = c00*EXI(jss) + c10*EXI(jss+1)

C     at R(Nza+1,js-1)
      c00 = s2*(D3/D8)*dz1*Dy(js-1)/dyy
      c10 = s2*(D1/D8)*dz1*Dy(js-1)/dyy
      dQdR1(2) = c00*EXI(jss) + c10*EXI(jss+1)

      RETURN
      END ! CompdQdR1_TE()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TE : Zxy = Ex/Hy
C          dZxy/dR  = [(Hy*dEx/dR-Ex*dHy/dR)/(Hy^2)]
C          dQdR = dQ1dR + dQ2dR


      SUBROUTINE CompdQdR0_TE(im,is,per,Nzb,Ny,StaNod,CRho,
     >               apte,phte,zzte,cos2,DFStatus,dQdR1,dQdR2,EA,EB)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER DFStatus(MM0MX)
      REAL*8  CRho(NZ0MX,NY0MX)
      REAL*8  apte,phte,cos2,per,EA(*),EB(*)
      COMPLEX*16 dQdR1(*),dQdR2(*),zzte

      INTEGER jj,mmt,js,iy,iz
      REAL*8  omue,omega,dRdR,dPdR
      COMPLEX*16 dQdR(MM0MX),dZdR

      omega = (D2*Pi)/per
      omue  = omega*Mue
      mmt = Ny*Nzb
      js  = StaNod(im,is)

      CALL CopyVectorC16(1,mmt,dQdR2,1,mmt,dQdR)
      iz = 1
      iy = js
      jj = (iy-1)*Nzb + iz
      dQdR(jj) = dQdR(jj) + dQdR1(1)

      iz = 1
      iy = js-1
      jj = (iy-1)*Nzb + iz
      dQdR(jj) = dQdR(jj) + dQdR1(2)

      CALL ConstantVectorR8(EA,mmt,D0)
      CALL ConstantVectorR8(EB,mmt,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj   = jj+1
          dZdR = dQdR(jj)
          dRdR = (D2/omue)*DREAL(CONJG(zzte)*dZdR)
          dPdR = cos2*(IMAG(dZdR)-DTAN(phte)*DREAL(dZdR))

          EA(jj) = (CRho(iz,iy)*dRdR/apte)*DFStatus(jj)
          EB(jj) = (CRho(iz,iy)*dPdR*DLOG(D10))*DFStatus(jj)
        ENDDO ! iz
      ENDDO ! iy

      RETURN
      END ! CompdQdR0_TE()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE Sens2DTipper(im,ip,per,Nza,Nz,Ny,Dz,Dy,CRho,
     >           NRes,ResTyp,NSta,StaNod,SenInx,DatInx,DFStatus,
     >           dATE,AII,ipiv,EXI,EXB,Hys,Hzs,
     >           CmHtIndx,CmHt)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz,Ny,im,ip,NSta(*),NRes(*)
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER ipiv(MMIMX)
      INTEGER DFStatus(MM0MX)
      COMPLEX*16 AII(NZ3MX,MMIMX)
      COMPLEX*16 EXI(*),EXB(*),Hzs(*),Hys(*)
      REAL*8  dATE(MM0MX,4)
      INTEGER StaNod(NMODMX,NSTAMX)
    
      INTEGER CmHtIndx(NMODMX,NRESMX,NSTAMX)
      REAL*8  CmHt(MM0MX,LL0MX),EA(MM0MX),EB(MM0MX)
      COMPLEX*16 dLR(MM0MX),dLL(MM0MX),dUL(MM0MX),dUR(MM0MX)
      COMPLEX*16 dQdR1(2),dQdR2(MM0MX)

      
      INTEGER is,ir,mmt,iss(NRESMX)
      INTEGER docompute,js,jss,jsr,jsl,iz,kl,ku,info,mmi
      INTEGER j1s,j1r,j1l,Nzb,jj,igx
      REAL*8  omega,omue,dz1,dz2,dzz,dyy,sav
      COMPLEX*16 ci(MMIMX),dyz,hzhy

      INTEGER  DoSens

      omega = (D2*Pi)/per
      omue   = omega*Mue

      Nzb = Nz-Nza
      mmi = (Ny-1)*(Nz-1)
      mmt = Ny*Nzb


      dz1 = Dz(Nza+1)
      dz2 = Dz(Nza+2)
      dzz = dz1 + dz2
      CALL CompZI1_TE(Nza,Nz,Ny,CRho,DFStatus,EXI,dATE,
     >     dLR,dLL,dUL,dUR)

      kl = Nz-1
      ku = Nz-1

      DO ir = 1,NRes(im)
        iss(ir) = 0
      ENDDO ! ir
      DO is = 1,NSta(im)
        docompute = DoSens(im,ip,NRes,NSta,SenInx)

        IF (docompute.EQ.1) THEN
          CALL ConstantVectorC16(ci,mmi,D0)

          hzhy = -Hzs(is)/(Hys(is)**D2)

          js   = StaNod(im,is)
          dyy  = Dy(js) + Dy(js-1)
          sav  = (Dy(js)/CRho(1,js) + Dy(js-1)/CRho(1,js-1))/dyy
          dyz  = DCMPLX(D0,dz1/(dyy*omue))

          iz   = Nza+1
C         one node beneath the surface node ...
          jss  = (js-2)*(Nz-1)     + (iz-1)
          jsr  = ((js+1)-2)*(Nz-1) + (iz-1)
          jsl  = ((js-1)-2)*(Nz-1) + (iz-1)
          j1s  = (js-2)*(Nz-1)     + (iz)
          j1r  = ((js+1)-2)*(Nz-1) + (iz)
          j1l  = ((js-1)-2)*(Nz-1) + (iz)

          ci(jss) = hzhy*(DCMPLX(D0,D1/(omue*dz1)) + 
     >      (D3/D8)*sav*dz1 + (D3/D4)*dyz*(D1/Dy(js) + D1/Dy(js-1)))
          ci(jsr) = (D1/Hys(is))*DCMPLX(D0,D1/(omue*dyy)) +
     >              hzhy*(-D3/D4)*dyz*(D1/Dy(js))
          ci(jsl) = (D1/Hys(is))*DCMPLX(D0,-D1/(omue*dyy)) +
     >              hzhy*(-D3/D4)*dyz*(D1/Dy(js-1))
          ci(j1s) = hzhy*(DCMPLX(D0,-D1/(omue*dz1)) +
     >      (D1/D8)*sav*dz1 + (D1/D4)*dyz*(D1/Dy(js) + D1/Dy(js-1)))
          ci(j1r) = hzhy*(-D1/D4)*dyz*(D1/Dy(js))
          ci(j1l) = hzhy*(-D1/D4)*dyz*(D1/Dy(js-1))

          IF ((js.LT.1).or.(js.GT.Ny)) THEN 
            WRITE(6,*) '!!! ATTENTION, ERROR ci TP CASE !!! '
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         solve A^{T}*ui = ci wher ci indicate location node
          CALL ZGBTRS('T',mmi,kl,ku,1,AII,NZ3MX,ipiv,ci,MMIMX,info)
          IF (INFO.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR SOLVING SENS2D_TP ',info
            WRITE(6,*) '!!! Please, correct input file and restart!!!'
            STOP
          ENDIF

C         dQdR0  =  (dai/dR)*XI+(dab/dR)*XB
C                   +------ dQdR1 --------+
C                   + -ui^{T}*[(dAi/dR)*XI+(dAb/dR)*Xb + (AIB)*(dXB/dR)]
C                             +---------- ZI1 -------+   +---- ZI2 ----+
C                    +--------------- dQdR2 ---------------------------+
C         assumes ZI2 = 0;

          CALL CompdQdR1_Tipper(im,is,Nza,Nz,Ny,Dz,Dy,CRho,Hys,Hzs,
     >         StaNod,EXI,dQdR1)
          CALL CompdQdR2_TE(per,Nza,Nz,Ny,dLR,dLL,dUL,dUR,ci,dQdR2)
          CALL CompdQdR0_Tipper(im,is,Nzb,Ny,StaNod,CRho,
     >                          DFStatus,dQdR1,dQdR2,EA,EB)

          DO ir = 1,NRes(im)
            igx = CmHtIndx(im,ir,ip)
            IF (SenInx(im,ir,ip,is).EQ.1) THEN
              iss(ir) = iss(ir) + 1
              DO jj = 1,mmt
                IF (ResTyp(im,ir).EQ.5) THEN 
                  CmHt(jj,igx+iss(ir)) = EA(jj)
                ENDIF
                IF (ResTyp(im,ir).EQ.6) THEN 
                  CmHt(jj,igx+iss(ir)) = EB(jj)
                ENDIF
              ENDDO
            ENDIF
          ENDDO ! ir

        ENDIF ! docompute
      ENDDO ! is


      RETURN
      END ! Sens2DTipper

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     dQdR1 = (-hzs(js)/hys(js)**D2)*[(dai/dR)*XI + (dab/dR)*XB]


      SUBROUTINE CompdQdR1_Tipper(im,is,Nza,Nz,Ny,Dz,Dy,CRho,Hys,Hzs,
     >                            StaNod,EXI,dQdR1)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nza,Nz,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      REAL*8  Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 EXI(*),dQdR1(*),Hzs(*),Hys(*)

      INTEGER iz,jss,js
      REAL*8  dz1,dz2,dzz,dyy
      COMPLEX*16 s2,c00,c10


      dz1 = Dz(Nza+1)
      dz2 = Dz(Nza+2)
      dzz = dz1 + dz2
      js  = StaNod(im,is)
      dyy = Dy(js) + Dy(js-1)
      s2  = (-Hzs(is)/(Hys(is)**D2))*(-D1/(CRho(1,js)**D2))

      iz  = Nza + 1
      jss = (js-2)*(Nz-1) + (iz-1)

C     at R(Nza+1,js)
      c00 = s2*(D3/D8)*dz1*Dy(js)/dyy
      c10 = s2*(D1/D8)*dz1*Dy(js)/dyy
      dQdR1(1) = c00*EXI(jss) + c10*EXI(jss+1)

C     at R(Nza+1,js-1)
      c00 = s2*(D3/D8)*dz1*Dy(js-1)/dyy
      c10 = s2*(D1/D8)*dz1*Dy(js-1)/dyy
      dQdR1(2) = c00*EXI(jss) + c10*EXI(jss+1)

      RETURN
      END ! CompdQdR1_Tipper()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     HZ/HY : Tipper = Hz/Hy
C          dTip/dR  = [(Hy*dHz/dR-Hz*dHy/dR)/(Hy^2)]
C          dQdR = dQ1dR + dQ2dR


      SUBROUTINE CompdQdR0_Tipper(im,is,Nzb,Ny,StaNod,CRho,
     >           DFStatus,dQdR1,dQdR2,EA,EB)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,im,is
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER DFStatus(MM0MX)
      REAL*8  CRho(NZ0MX,NY0MX)
      REAL*8  EA(*),EB(*)
      COMPLEX*16 dQdR1(*),dQdR2(*)

      INTEGER jj,mmt,js,iy,iz
      COMPLEX*16 dQdR(MM0MX)

      mmt = Ny*Nzb
      js  = StaNod(im,is)

      CALL CopyVectorC16(1,mmt,dQdR2,1,mmt,dQdR)
      iz = 1
      iy = js
      jj = (iy-1)*Nzb + iz
      dQdR(jj) = dQdR(jj) + dQdR1(1)

      iz = 1
      iy = js-1
      jj = (iy-1)*Nzb + iz
      dQdR(jj) = dQdR(jj) + dQdR1(2)

      CALL ConstantVectorR8(EA,mmt,D0)
      CALL ConstantVectorR8(EB,mmt,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = jj+1
          EA(jj) = (CRho(iz,iy)*DLOG(D10)*DREAL(dQdR(jj)))*DFStatus(jj)
          EB(jj) = (CRho(iz,iy)*DLOG(D10)*IMAG(dQdR(jj)))*DFStatus(jj)

        ENDDO
      ENDDO ! iy

      RETURN
      END ! CompdQdR0_Tipper()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
