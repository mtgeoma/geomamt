
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupA_TM(Nzb,Ny,Dzb,Dy,Czb,Cy,CRho,ATM,BTM)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      REAL*8  Dzb(*),Dy(*),Czb(*),Cy(*),CRho(NZ0MX,NY0MX)
      REAL*8  ATM(MMIMX,4),BTM(MMBMX)
      
C     ATM(:,1) : Real diagonal term
C     ATM(:,2) : Real first upper strip
C     ATM(:,3) : Real second upper strip
C     ATM(:,4) : Imaginary diagonal term
C     BTM      : Boundary

      INTEGER jj,iz,iy,nblx,iblx
      REAL*8  r00,r10,r01,r11
      REAL*8  ar,al,ad,au

      nblx = Ny*Nzb
      iblx = 2*Ny + 2*Nzb

      CALL ConstantMatrixR8(ATM,MMIMX,4,nblx,4,D0)
      CALL ConstantVectorR8(BTM,iblx,D0)

      jj = 1
      DO iy = 2,Ny
        DO iz = 2,Nzb
          r00 = CRho(iz,iy)
          r10 = CRho(iz-1,iy)
          r01 = CRho(iz,iy-1)
          r11 = CRho(iz-1,iy-1)

          ar  = D2*(Dzb(iz)*r00 + Dzb(iz-1)*r10)/Dy(iy)
          al  = D2*(Dzb(iz)*r01 + Dzb(iz-1)*r11)/Dy(iy-1)
          ad  = D2*(Dy(iy)*r00  + Dy(iy-1)*r01)/Dzb(iz)
          au  = D2*(Dy(iy)*r10  + Dy(iy-1)*r11)/Dzb(iz-1)
          ATM(jj,1) = - ar - al - ad - au
          ATM(jj,4) = Cy(iY)*Czb(iz)
          ATM(jj,2) = ad
          IF (iz.EQ.Nzb) THEN
            ATM(jj,2) = D0
          ENDIF
          IF (iy.LT.NY) THEN 
            ATM(jj,3) = ar
          ENDIF
          jj = jj + 1
        ENDDO ! iz
      ENDDO ! iy

C     boundary fields

C     left side fields
      iy = 2
      DO iz = 2,Nzb
        r01 = CRho(iz,iy-1)
        r11 = CRho(iz-1,iy-1)
        BTM(iz) = D2*(Dzb(iz)*r01 + Dzb(iz-1)*r11)/Dy(iy-1)
      ENDDO ! iz

C     surface fields
      iz = 2
      DO iy = 2,Ny
        r10 = CRho(iz-1,iy)
        r11 = CRho(iz-1,iy-1)
        BTM(Nzb+iy) = D2*(Dy(iy)*r10  + Dy(iy-1)*r11)/Dzb(iz-1)
      ENDDO ! iy

C     Bottom fields
      iz = Nzb
      DO iy = 2,Ny
        r00 = CRho(iz,iy)
        r01 = CRho(iz,iy-1)
        BTM(Nzb+Ny+iy-1) = D2*(Dy(iy)*r00  + Dy(iy-1)*r01)/Dzb(iz)
      ENDDO ! iy

C     Right fields
      iy = Ny
      DO iz = 2,Nzb
        r00 = CRho(iz,iy)
        r10 = CRho(iz-1,iy)
        BTM(Nzb+2*Ny-1+iz) = D2*(Dzb(iz)*r00 + Dzb(iz-1)*r10)/Dy(iy)
      ENDDO ! iz


      RETURN
      END ! SetupA_TM()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,CRho,ATE,BTE)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz,Ny
      REAL*8  Dz(*),Dy(*),Cz(*),Cy(*),CRho(NZ0MX,NY0MX)
      REAL*8  ATE(MMIMX,4),BTE(MMBMX)
      
C     ATE(:,1) : Real diagonal term
C     ATE(:,2) : Real second strip
C     ATE(:,3) : Real thrid strip
C     ATE(:,4) : Imaginary diagonal term
C     BTE      : Boundary

      INTEGER jj,iz,iy,nblx,iblx
      REAL*8  s00,s10,s01,s11,s2d
      REAL*8  ar,al,ad,au,CCon(NZ0MX,NY0MX)

      nblx = Ny*Nz
      iblx = 2*Ny + 2*Nz

      CALL ConstantMatrixR8(CCon,NZ0MX,NY0MX,Nza,Ny,CondAir)
      DO iy = 1,Ny
        DO iz = Nza+1,Nz
          CCon(iz,iy) = 1./CRho(iz-Nza,iy)
        ENDDO
      ENDDO

      CALL ConstantMatrixR8(ATE,MMIMX,4,nblx,4,D0)
      CALL ConstantVectorR8(BTE,iblx,D0)

      jj = 1
      DO iy = 2,Ny
        DO iz = 2,Nz
          s00 = CCon(iz,iy)
          s10 = CCon(iz-1,iy)
          s01 = CCon(iz,iy-1)
          s11 = CCon(iz-1,iy-1)
          s2d = (s00*Dz(iz)*Dy(iy)   + s01*Dz(iz)*Dy(iy-1) +
     >           s10*Dz(iz-1)*Dy(iy) + s11*Dz(iz-1)*Dy(iy-1))/
     >          (Dz(iz)*Dy(iy)   + Dz(iz)*Dy(iy-1) + 
     >           Dz(iz-1)*Dy(iy) + Dz(iz-1)*Dy(iy-1)) 

          ar  = D2*Cz(iz)/Dy(iy)
          al  = D2*Cz(iz)/Dy(iy-1)
          ad  = D2*Cy(iy)/Dz(iz)
          au  = D2*Cy(iy)/Dz(iz-1)
          ATE(jj,1) = - ar - al - ad - au
          ATE(jj,4) = Cy(iy)*Cz(iz)*s2d

          ATE(jj,2) = ad
          IF (iz.EQ.Nz) THEN
            ATE(jj,2) = D0
          ENDIF
          IF (iy.LT.NY) THEN 
            ATE(jj,3) = ar
          ENDIF
          jj = jj + 1
        ENDDO ! iz
      ENDDO ! iy

C     boundary fields

C     left side fields
      iy = 2
      DO iz = 2,Nz
        BTE(iz) = D2*Cz(iz)/Dy(iy-1)
      ENDDO ! iz

C     surface fields
      iz = 2
      DO iy = 2,Ny
        BTE(Nz+iy) = D2*Cy(iy)/Dz(iz-1)
      ENDDO ! iy

C     Bottom fields
      iz = Nz
      DO iy = 2,Ny
        BTE(Nz+Ny+iy-1) = D2*Cy(iy)/Dz(iz)
      ENDDO ! iy

C     Right fields
      iy = Ny
      DO iz = 2,Nz
        BTE(Nz+2*Ny-1+iz) = D2*Cz(iz)/Dy(iy)
      ENDDO ! iz


      RETURN
      END ! SetupA_TE()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupDiffOper(NMode,ModTyp,Nza,Nzb,Nz,Ny,
     >                         Dzb,Dz,Dy,Cz,Cy,
     >                         dATM,dATE)
      INCLUDE 'parameter.h'

      INTEGER NMode,ModTyp(*),Nzb,Ny,Nza,Nz
      REAL*8  Dy(*),Dzb(*),Cz(*),Cy(*),Dz(*)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)

      INTEGER im,done_te

      done_te = 0
      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) 
     >     CALL SetupDiffOper_TM(Nzb,Ny,Dzb,Dy,dATM)
        IF ((ModTyp(im).EQ.2).OR.(ModTyp(im).EQ.3)) THEN
          IF (done_te.EQ.0) 
     >      CALL SetupDiffOper_TE(Nza,Nzb,Nz,Ny,Dz,Dy,Cz,Cy,dATE)
          done_te = 1
        ENDIF
      ENDDO

      RETURN
      END ! SetupDiffOper()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupDiffOper_TM(Nzb,Ny,Dzb,Dy,dATM)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      REAL*8  Dzb(*),Dy(*)
      REAL*8  dATM(MM0MX,3)

      INTEGER jj,iz,iy
 
      CALL ConstantMatrixR8(dATM,MM0MX,3,Ny*Nzb,3,D0)
      jj = 1
      DO iy = 1,Ny
        DO iz = 1,Nzb
C         left or right
          dATM(jj,2) = D2*Dzb(iz)/Dy(iy)
C         up or down
          dATM(jj,3) = D2*Dy(iy)/Dzb(iz)
          dATM(jj,1) = - dATM(jj,2) - dATM(jj,3)
          jj = jj + 1
        ENDDO ! izb
      ENDDO ! iy

      RETURN
      END ! SetupDiffOper_TM()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetupDiffOper_TE(Nza,Nzb,Nz,Ny,Dz,Dy,Cz,Cy,
     >                            dATE)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nzb,Nz,Ny
      REAL*8  Dz(*),Dy(*),Cz(*),Cy(*)
      REAL*8  dATE(MM0MX,4)

      INTEGER jj,iz,iy

C     Upper Left  : dATE(jj,1)
C     Upper Right : dATE(jj,2)
C     Lower Left  : dATE(jj,3)
C     Lower Right : dATE(jj,4)

      CALL ConstantMatrixR8(dATE,MM0MX,4,Ny*Nzb,4,D0)

      DO iy = 2,Ny-1
        DO iz = Nza+1,Nz-1
          jj = (iy-1)*Nzb + iz-Nza
          CALL DiffTE(iz,iy,iz,iy,Dz,Dy,Cz,Cy,dATE(jj,1))
          CALL DiffTE(iz,iy,iz,iy+1,Dz,Dy,Cz,Cy,dATE(jj,2))
          CALL DiffTE(iz,iy,iz+1,iy,Dz,Dy,Cz,Cy,dATE(jj,3))
          CALL DiffTE(iz,iy,iz+1,iy+1,Dz,Dy,Cz,Cy,dATE(jj,4))
        ENDDO
      ENDDO

      iy = 1
      DO iz = Nza+1,Nz-1
        jj = (iy-1)*Nzb + iz-Nza
        CALL DiffTE(iz,iy,iz,iy+1,Dz,Dy,Cz,Cy,dATE(jj,2))
        CALL DiffTE(iz,iy,iz+1,iy+1,Dz,Dy,Cz,Cy,dATE(jj,4))
      ENDDO

      iy = Ny
      DO iz = Nza+1,Nz-1
        jj = (iy-1)*Nzb + iz-Nza
        CALL DiffTE(iz,iy,iz,iy,Dz,Dy,Cz,Cy,dATE(jj,1))
        CALL DiffTE(iz,iy,iz+1,iy,Dz,Dy,Cz,Cy,dATE(jj,3))
      ENDDO

      iz = Nz
      iy = 1
      jj = (iy-1)*Nzb + iz-Nza
      CALL DiffTE(iz,iy,iz,iy+1,Dz,Dy,Cz,Cy,dATE(jj,2))

      iz = Nz
      iy = Ny
      jj = (iy-1)*Nzb + iz-Nza
      CALL DiffTE(iz,iy,iz,iy,Dz,Dy,Cz,Cy,dATE(jj,1))

      iz = Nz
      DO iy = 2,Ny-1
        jj = (iy-1)*Nzb + iz-Nza
        CALL DiffTE(iz,iy,iz,iy,Dz,Dy,Cz,Cy,dATE(jj,1))
        CALL DiffTE(iz,iy,iz,iy+1,Dz,Dy,Cz,Cy,dATE(jj,2))
      ENDDO

      RETURN
      END ! SetupDiffOper_TE()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE DiffTE(iz,iy,kz,ky,Dz,Dy,Cz,Cy,dA)
      INTEGER ky,kz,iy,iz
      REAL*8  Dz(*),Dy(*),Cz(*),Cy(*),dA

      REAL*8  dv,dv00,dv01,dv10,dv11

      dv00 = Dz(kz)*Dy(ky)
      dv01 = Dz(kz)*Dy(ky-1)
      dv10 = Dz(kz-1)*Dy(ky)
      dv11 = Dz(kz-1)*Dy(ky-1)
      dv   = dv00 + dv01 + dv10 + dv11
      dA   = Cz(kz)*Cy(ky)*Dy(iy)*Dz(iz)/dv

      RETURN
      END ! DiffTE_UpperLeft()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
