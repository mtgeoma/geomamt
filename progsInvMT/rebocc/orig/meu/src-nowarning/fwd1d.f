C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb
      REAL*8  per,Dzb(*),r1d(*)
      COMPLEX*16 X1D(*)
      
      INTEGER iz,jj,nzb1,nzb2
      REAL*8  Omega,Omue,skd
      REAL*8  dz1(NZ1MX),cz1(NZ2MX)
      REAL*8  au1d(NZ0MX),ad1d(NZ0MX)
      COMPLEX*16 ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

      nzb1 = Nzb + 1
      nzb2 = Nzb + 2
      Omega = (D2*Pi)/Per
      Omue  = Omega*Mue
 
      r1d(nzb1) = r1d(Nzb)
      skd  = DSQRT(D2*r1d(nzb1)/Omue)
      CALL CopyVectorR8(1,Nzb,Dzb,1,Nzb,dz1)
      dz1(nzb1) = skd
      CALL DistanceBetweenBlocks(nzb1,dz1,cz1)

C     Assign operators
      CALL ConstantVectorR8(au1d,nzb1,D0)
      CALL ConstantVectorR8(ad1d,nzb1,D0)
      CALL ConstantVectorC16(ac1d,nzb1,D0)
      jj = 1
      DO iz = 2,nzb1
        au1d(jj) = D2*r1d(iz-1)/dz1(iz-1)
        ad1d(jj) = D2*r1d(iz)/dz1(iz)
        ac1d(jj) = DCMPLX(D0,Omue)*cz1(iz) - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

C     Boundary condition for 1D
      xb1d(1) = D1
      xb1d(2) = D0

C     solve Aii*Xi = -Aib*Xb
      CALL ConstantVectorC16(xb,nzb,D0)
      xb(1) = -au1d(1)*xb1d(1)

      CALL Solve1D(nzb,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      CALL CopyVectorC16(1,nzb,xb,2,nzb1,X1D)

c100   FORMAT(7e11.3)

      Return
      END ! Fwd1d_TM

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd1D_TE(per,Nza,Nz,Dz,s1,X1D)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz
      REAL*8  per,Dz(*),s1(*)
      COMPLEX*16 X1D(*)
      
      INTEGER iz,jj,nz1,nz2
      REAL*8  Omega,Omue,skd,s1d(NZ1MX)
      REAL*8  dz1(NZ1MX),cz1(NZ2MX)
      REAL*8  au1d(NZ0MX),ad1d(NZ0MX),ss
      COMPLEX*16 ac1d(NZ0MX),xb1d(2),xb(NZ0MX)

      nz1 = Nz + 1
      nz2 = Nz + 2
      Omega = (D2*Pi)/Per
      Omue  = Omega*Mue

      DO iz = 1,Nza
        s1d(iz) = CondAir
      ENDDO
      DO iz = Nza+1,Nz
        s1d(iz) = s1(iz-Nza)
      ENDDO
      s1d(nz1) = s1d(Nz)

      skd  = DSQRT(D2/(s1d(nz1)*Omue))
      DO iz = 1,Nz
        dz1(iz) = Dz(iz)
      ENDDO ! iz

c     CALL CopyVectorR8(1,Nz,Dz,1,Nz,dz1)
      dz1(nz1) = skd
      CALL DistanceBetweenBlocks(nz1,dz1,cz1)

C     Assign operators
      CALL ConstantVectorR8(au1d,nz1,D0)
      CALL ConstantVectorR8(ad1d,nz1,D0)
      CALL ConstantVectorC16(ac1d,nz1,D0)
      jj = 1
      DO iz = 2,nz1
        au1d(jj) = D2/dz1(iz-1)
        ad1d(jj) = D2/dz1(iz)
        ss       = (s1d(iz)*dz1(iz)+s1d(iz-1)*dz1(iz-1))
        ac1d(jj) = DCMPLX(D0,Omue*ss) - au1d(jj) - ad1d(jj)
        jj = jj + 1
      ENDDO

C     Boundary condition for 1D
      xb1d(1) = D1
      xb1d(2) = D0

C     solve Aii*Xi = -Aib*Xb
      CALL ConstantVectorC16(xb,Nz,D0)
      xb(1) = -au1d(1)*xb1d(1)

      CALL Solve1D(Nz,ad1d,ac1d,xb)
      X1D(1) = xb1d(1)
      CALL CopyVectorC16(1,nz,xb,2,nz1,X1D)


c100   FORMAT(7e11.3)

      Return
      END ! Fwd1d_TE

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Solve1D(nz0,ad1d,ac1d,xb)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER nz0
      REAL*8  ad1d(*)
      COMPLEX*16 ac1d(*),xb(*)

      INTEGER jj,ipiv(NZ0MX),info
      COMPLEX*16 a1d(4,NZ0MX)
     
C     make matrix in LAPACK's format
      CALL ConstantMatrixC16(a1d,4,NZ0MX,4,nz0,D0)

      DO jj = 2,nz0
        a1d(2,jj) = ad1d(jj-1)
      ENDDO
      DO jj = 1,nz0
        a1d(3,jj) = ac1d(jj)
      ENDDO
      DO jj = 1,nz0-1
        a1d(4,jj) = ad1d(jj)
      ENDDO

      CALL ZGBTRF(nz0,nz0,1,1,a1d,4,ipiv,info)
      CALL ZGBTRS('N',nz0,1,1,1,a1d,4,ipiv,xb,NZ0MX,info)

      END ! Solve1D

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
