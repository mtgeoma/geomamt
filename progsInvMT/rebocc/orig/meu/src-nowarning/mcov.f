C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.


      SUBROUTINE Setup1DCM(Nzb,Ny,Dzb,Dy,Czb,Cy,ModGrd,HGamma,VGamma,
     >           HDiff,VDiff)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      INTEGER ModGrd(NZ0MX,NY0MX)
      REAL*8  Dzb(*),Dy(*),Czb(*),Cy(*)
      REAL*8  HGamma(NZ0MX,NY0MX),VGamma(*)
      REAL*8  HDiff(2,NZ0MX,NY0MX),VDiff(2,NZ0MX,NY0MX)

      INTEGER mgd(MM0MX),iy,iz
      REAL*8  hdif(2,NY0MX),vdif(2,NZ0MX),hgam(NY0MX),vgam(NZ0MX)


C     Horizontal Diffusion
      DO iz = 1,Nzb
        DO iy = 1,Ny
          mgd(iy)  = ModGrd(iz,iy)
          hgam(iy) = HGamma(iz,iy)
        ENDDO ! iy
        CALL SetupHorCM(Ny,Dy,Cy,mgd,hgam,hdif)
        DO iy = 1,Ny
          HDiff(1,iz,iy) = hdif(1,iy)
          HDiff(2,iz,iy) = hdif(2,iy)
        ENDDO ! iy
      ENDDO ! iz

C     Vertical Diffusion
      DO iy = 1,Ny
        DO iz = 1,Nzb
          mgd(iz)  = ModGrd(iz,iy)
          vgam(iz) = VGamma(iz)
        ENDDO ! iz
        CALL SetupVerCM(Nzb,Dzb,Czb,mgd,vgam,vdif)
        DO iz = 1,Nzb
          VDiff(1,iz,iy) = vdif(1,iz)
          VDiff(2,iz,iy) = vdif(2,iz)
        ENDDO ! iz
      ENDDO ! iy

      RETURN
      END ! SUBROUTINE Setup1DCM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupHorCM(Ny,Dy,Cy,mgd,hgam,hdif)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,mgd(*)
      REAL*8  Dy(*),Cy(*),hgam(*),hdif(2,NY0MX)

      INTEGER iy,info
      REAL*8  cr,cl,ygam,DHC(NY0MX),DHR(NY0MX),DHL(NY0MX)

C     Horizontal Diffusion
      CALL ConstantVectorR8(DHC,Ny,D0)
      CALL ConstantVectorR8(DHR,Ny,D0)
      CALL ConstantVectorR8(DHL,Ny,D0)

      DO iy = 1,Ny
        cr = D0
        cl = D0
        IF (iy.LT.Ny) THEN
          ygam = (hgam(iy)*Dy(iy)+hgam(iy+1)*Dy(iy+1))/
     >           (Dy(iy)+Dy(iy+1))
          cr = (D2*ygam)/Cy(iy+1)

          IF (((mgd(iy).EQ.2).OR.(mgd(iy).EQ.6)).OR.
     >         (mgd(iy).EQ.8)) THEN
            cr = D0
          ENDIF
        ENDIF
        IF (iy.GT.1) THEN
          ygam = (hgam(iy)*Dy(iy)+hgam(iy-1)*Dy(iy-1))/
     >           (Dy(iy)+Dy(iy-1))
          cl = (D2*ygam)/Cy(iy)

          IF (((mgd(iy).EQ.3).OR.(mgd(iy).EQ.7)).OR.
     >         (mgd(iy).EQ.9)) THEN
            cl = D0
          ENDIF
        ENDIF

        DHC(iy) = Dy(iy)
        DHC(iy) = DHC(iy) + cr + cl
        DHR(iy) = -cr
        DHL(iy) = -cl

c100     CONTINUE
      ENDDO ! iy

      CALL ConstantMatrixR8(hdif,2,NY0MX,2,Ny,D0)
      DO iy = 2,Ny
        hdif(1,iy) = DHR(iy-1)
      ENDDO
      DO iy = 1,Ny
        hdif(2,iy) = DHC(iy)
      ENDDO

      CALL DPBTRF('U',Ny,1,hdif,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : HDiff',info
        STOP
      ENDIF

      RETURN
      END ! SUBROUTINE SetupHorCM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupVerCM(Nzb,Dzb,Czb,mgd,vgam,vdif)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,mgd(*)
      REAL*8  Dzb(*),Czb(*),vgam(*),vdif(2,NZ0MX)

      INTEGER iz,info
      REAL*8  cd,cu,zgam,DVC(NZ0MX),DVD(NZ0MX),DVU(NZ0MX)


C     Vertical Diffusion
      CALL ConstantVectorR8(DVC,Nzb,D0)
      CALL ConstantVectorR8(DVD,Nzb,D0)
      CALL ConstantVectorR8(DVU,Nzb,D0)

      DO iz = 1,Nzb
        cd = D0
        cu = D0
        IF (iz.LT.Nzb) THEN
          zgam = (vgam(iz)*Dzb(iz)+vgam(iz+1)*Dzb(iz+1))/
     >           (Dzb(iz)+Dzb(iz+1))
          cd = (D2*zgam)/Czb(iz+1)

          IF (((mgd(iz).EQ.5).OR.(mgd(iz).EQ.8)).OR.
     >         (mgd(iz).EQ.9)) THEN
            cd = D0
          ENDIF
        ENDIF
        IF (iz.GT.1) THEN
          zgam = (vgam(iz)*Dzb(iz)+vgam(iz-1)*Dzb(iz-1))/
     >           (Dzb(iz)+Dzb(iz-1))
          cu = (D2*zgam)/Czb(iz)
          IF (((mgd(iz).EQ.4).OR.(mgd(iz).EQ.6)).OR.
     >         (mgd(iz).EQ.7)) THEN
            cu = D0
          ENDIF
        ENDIF

        DVC(iz) = Dzb(iz)
        DVC(iz) = DVC(iz) + cd + cu
        DVD(iz) = -cd
        DVU(iz) = -cu

c100     CONTINUE
      ENDDO ! iz

      CALL ConstantMatrixR8(vdif,2,NZ0MX,2,Nzb,D0)
      DO iz = 2,Nzb
        vdif(1,iz) = DVD(iz-1)
      ENDDO
      DO iz = 1,Nzb
        vdif(2,iz) = DVC(iz)
      ENDDO

      CALL DPBTRF('U',Nzb,1,vdif,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : VDiff',info
        STOP
      ENDIF

      RETURN
      END ! SUBROUTINE SetupVerCM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SolveDiff(order,PDifTime,HDiff,VDiff,
     >           Nzb,Ny,Dzb,Dy,DFStatus,uu)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER PDifTime,Nzb,Ny,DFStatus(*),order
      REAL*8  uu(*),Dzb(*),Dy(*)
      REAL*8  HDiff(2,NZ0MX,NY0MX),VDiff(2,NZ0MX,NY0MX)

      INTEGER iy,iz,jj,dif2(NZ0MX,NY0MX),it
      REAL*8  u2(NZ0MX,NY0MX)

      jj = 1
      DO iy = 1,Ny
       DO iz = 1,Nzb
        u2(iz,iy)   = uu(jj)
        dif2(iz,iy) = DFStatus(jj)
        jj = jj + 1
       ENDDO
      ENDDO

      IF (order.EQ.1) THEN
        DO it = 1,PDifTime
C         Horizontal Diffusion
          CALL HorDiff(1,VDiff,Ny,Nzb,Dzb,dif2,u2)
C         Vertical Diffusion
          CALL VerDiff(1,HDiff,Ny,Nzb,Dy,dif2,u2)
        ENDDO 
      ENDIF
      IF (order.EQ.2) THEN
        DO it = 1,PDifTime
C         Vertical Diffusion
          CALL VerDiff(2,HDiff,Ny,Nzb,Dy,dif2,u2)
C         Horizontal Diffusion
          CALL HorDiff(2,VDiff,Ny,Nzb,Dzb,dif2,u2)
        ENDDO
      ENDIF

      jj = 1
      DO iy = 1,Ny
       DO iz = 1,Nzb
        uu(jj) = u2(iz,iy)
        jj = jj + 1
       ENDDO
      ENDDO

      RETURN
      END ! SolveDiff()
     
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE HorDiff(order,VDiff,Ny,Nzb,Dzb,dif2,u2)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,dif2(NZ0MX,NY0MX)
      REAL*8  VDiff(2,NZ0MX,NY0MX),Dzb(*)
      REAL*8  u2(NZ0MX,NY0MX)

      INTEGER iy,iz,info,order
      REAL*8  uz(NZ0MX),vdif(2,NZ0MX)

C     Horizontal Diffusion

      IF (order.EQ.1) THEN
        DO iy = 1,Ny
          DO iz = 1,Nzb
            uz(iz) = u2(iz,iy)
            vdif(1,iz) = VDiff(1,iz,iy)
            vdif(2,iz) = VDiff(2,iz,iy)
          ENDDO
          DO iz = 1,Nzb
            uz(iz) = uz(iz)*Dzb(iz)*dif2(iz,iy)
          ENDDO
          CALL DPBTRS('U',Nzb,1,1,vdif,2,uz,NZ0MX,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : VDiff',info
            STOP
          ENDIF
          DO iz = 1,Nzb
            u2(iz,iy) = uz(iz)*dif2(iz,iy)
          ENDDO
        ENDDO ! iy
      ENDIF

      IF (order.EQ.2) THEN
        DO iy = 1,Ny
          DO iz = 1,Nzb
            uz(iz) = u2(iz,iy)
            vdif(1,iz) = VDiff(1,iz,iy)
            vdif(2,iz) = VDiff(2,iz,iy)
          ENDDO
          CALL DPBTRS('U',Nzb,1,1,vdif,2,uz,NZ0MX,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : VDiff',info
            STOP
          ENDIF
          DO iz = 1,Nzb
            uz(iz) = uz(iz)*Dzb(iz)*dif2(iz,iy)
          ENDDO
          DO iz = 1,Nzb
            u2(iz,iy) = uz(iz)
          ENDDO
        ENDDO ! iy
      ENDIF

      RETURN
      END ! HorDiff

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE VerDiff(order,HDiff,Ny,Nzb,Dy,dif2,u2)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nzb,dif2(NZ0MX,NY0MX)
      REAL*8  HDiff(2,NZ0MX,NY0MX),Dy(*)
      REAL*8  u2(NZ0MX,NY0MX),hdif(2,NY0MX)

      INTEGER iy,iz,info,order
      REAL*8  uy(NY0MX)

C     Vertical Diffusion

      IF (order.EQ.1) THEN
        DO iz = 1,Nzb
          DO iy = 1,Ny
            uy(iy) = u2(iz,iy)
            hdif(1,iy) = HDiff(1,iz,iy)
            hdif(2,iy) = HDiff(2,iz,iy)
          ENDDO
          DO iy = 1,Ny
            uy(iy) = uy(iy)*Dy(iy)*dif2(iz,iy)
          ENDDO
          CALL DPBTRS('U',Ny,1,1,hdif,2,uy,NY0MX,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : HDiff',info
            STOP
          ENDIF
          DO iy = 1,Ny
            u2(iz,iy) = uy(iy)*dif2(iz,iy)
          ENDDO
        ENDDO ! iz
      ENDIF

      IF (order.EQ.2) THEN
        DO iz = 1,Nzb
          DO iy = 1,Ny
            uy(iy) = u2(iz,iy)
            hdif(1,iy) = HDiff(1,iz,iy)
            hdif(2,iy) = HDiff(2,iz,iy)
          ENDDO
          CALL DPBTRS('U',Ny,1,1,hdif,2,uy,NY0MX,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : HDiff',info
            STOP
          ENDIF
          DO iy = 1,Ny
            uy(iy) = uy(iy)*Dy(iy)*dif2(iz,iy)
          ENDDO
          DO iy = 1,Ny
            u2(iz,iy) = uy(iy)
          ENDDO
        ENDDO ! iz
      ENDIF

      RETURN
      END ! VerDiff

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupNormCM(Nzb,Ny,Dzb,Dy,HDiff,VDiff,PDifTime,
     >           CMNorm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER PDifTime,Nzb,Ny
      REAL*8  Dzb(*),Dy(*)
      REAL*8  HDiff(2,NZ0MX,NY0MX),VDiff(2,NZ0MX,NY0MX)
      REAL*8  CMNorm(*)

      INTEGER izz,iz,iyy,iy,it,info,jj
      REAL*8  uy(NY0MX),uz(NZ0MX)
      REAL*8  ynorm(NZ0MX,NY0MX),znorm(NZ0MX,NY0MX)
      REAL*8  hdif(2,NY0MX),vdif(2,NZ0MX)


C     Horizontal Diffusion
      DO iyy = 1,Ny
        DO izz = 1,Nzb
          vdif(1,izz) = VDiff(1,izz,iyy)
          vdif(2,izz) = VDiff(2,izz,iyy)
        ENDDO ! izz
        DO izz = 1,Nzb
          CALL ConstantVectorR8(uz,Nzb,D0)
          uz(izz) = D1
          DO it = 1,PDifTime
            DO iz = 1,Nzb
              uz(iz) = uz(iz)*Dzb(iz)
            ENDDO
            CALL DPBTRS('U',Nzb,1,1,vdif,2,uz,NZ0MX,info)
            IF (info.NE.0) THEN
              WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
              WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : VDiff',
     >                   info
              STOP
            ENDIF
          ENDDO ! it
          DO it = 1,PDifTime
            CALL DPBTRS('U',Nzb,1,1,vdif,2,uz,NZ0MX,info)
            IF (info.NE.0) THEN
              WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
              WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : VDiff',
     >                   info
              STOP
            ENDIF
            DO iz = 1,Nzb
              uz(iz) = uz(iz)*Dzb(iz)
            ENDDO
          ENDDO ! it
          znorm(izz,iyy) = uz(izz)
        ENDDO ! izz
      ENDDO ! iyy

C     Vertical Diffusion
      DO izz = 1,Nzb
        DO iyy = 1,Ny
          hdif(1,iyy) = HDiff(1,izz,iyy)
          hdif(2,iyy) = HDiff(2,izz,iyy)
        ENDDO
        DO iyy = 1,Ny
          CALL ConstantVectorR8(uy,Ny,D0)
          uy(iyy) = D1
          DO it = 1,PDifTime
            DO iy = 1,Ny
              uy(iy) = uy(iy)*Dy(iy)
            ENDDO
            CALL DPBTRS('U',Ny,1,1,hdif,2,uy,NY0MX,info)
            IF (info.NE.0) THEN
              WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
              WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : HDiff',
     >                   info
              STOP
            ENDIF
          ENDDO ! it
          DO it = 1,PDifTime
            CALL DPBTRS('U',Ny,1,1,hdif,2,uy,NY0MX,info)
            IF (info.NE.0) THEN
              WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
              WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : HDiff',
     >                   info
              STOP
            ENDIF
            DO iy = 1,Ny
              uy(iy) = uy(iy)*Dy(iy)
            ENDDO
          ENDDO ! it
          ynorm(izz,iyy) = uy(iyy)
        ENDDO ! iyy
      ENDDO ! izz
 
      jj = 1
      DO iy = 1,Ny
       DO iz = 1,Nzb
         CMNorm(jj) = ynorm(iz,iy)*znorm(iz,iy)
         jj = jj+1
       ENDDO ! iz
      ENDDO ! iy

      RETURN
      END ! SetupNormCM
   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
