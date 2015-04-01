
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd2DTM_LU(per,Nzb,Ny,Dzb,Dy,CRho,ATM,BTM,
     >                      im,NSta,StaNod,
     >                      AII,ipiv,HXI,HXB,App,Phs,Zyx,Hxs,Eys)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER Nzb,Ny,StaNod(NMODMX,NSTAMX),im,NSta(*)
      REAL*8  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      REAL*8  ATM(MMIMX,4),BTM(MMBMX),App(*),Phs(*)
      COMPLEX*16 HXI(*),HXB(*),Zyx(*),Hxs(*),Eys(*)

      INTEGER ipiv(*),mmi,mmb,kl,ku,info,is
      REAL*8  AppTM(NY1MX),PhsTM(NY1MX)
      COMPLEX*16 AII(NZ3MX,MMIMX),HXI0(MMIMX),Ey0(NY1MX),
     >           Ztm(NY1MX),Hx0(NY1MX)

      mmi = (Ny-1)*(Nzb-1)
      mmb = 2*Ny + 2*Nzb
      kl = Nzb-1
      ku = Nzb-1
      
      CALL SetBound2D_TM(per,Nzb,Ny,Dzb,Dy,CRho,HXI0,HXB)
      CALL FormAII(per,Nzb,Ny,ATM,AII,ipiv)
      CALL MulAibWithXb(Nzb,Ny,BTM,HXB,HXI)
      CALL ZGBTRS('N',mmi,kl,ku,1,AII,NZ3MX,ipiv,HXI,MMIMX,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR WHILE SOLVING FWD !!!'
        WRITE(6,*) '!!! Please, check your model and restart !!!'
        WRITE(6,*) 
     >  '!!! If problem persists, contact wsiripun@oce.orst.edu !!'
        STOP
      ENDIF
      CALL CompEyAtSurface(per,Nzb,Ny,Dzb,Dy,HXI,HXB,CRho,Ey0,Hx0)
      CALL CompRespTM(per,Ny,Ey0,Ztm,AppTM,PhsTM)

      DO is = 1,NSta(im)
        App(is) = DLOG10(AppTM(StaNod(im,is)))
        Phs(is) = PhsTM(StaNod(im,is))*(Pi/D180)
        Zyx(is) = Ztm(StaNod(im,is))
        Hxs(is) = Hx0(StaNod(im,is))
        Eys(is) = Ey0(StaNod(im,is))
      ENDDO
      
c100   FORMAT(i5,6e15.7)
c200   FORMAT(i5,6e11.3)
c300   FORMAT(2i5,8e11.3)

      END ! Fwd2DTM_LU

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd2DTM_PCG(per,Nzb,Ny,Dzb,Dy,CRho,ATM,BTM,
     >                      im,NSta,StaNod,ETOL,NtMx,fwdcond,
     >                      HXI,HXB,App,Phs,Zyx,Hxs,Eys)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER Nzb,Ny,StaNod(NMODMX,NSTAMX),im,NSta(*),NtMx,
     >        fwdcond
      REAL*8  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX),ETOL
      REAL*8  ATM(MMIMX,4),BTM(MMBMX),App(*),Phs(*)
      COMPLEX*16 HXI(*),HXB(*),Zyx(*),Hxs(*),Eys(*)

      INTEGER mmi,is,iter
      REAL*8  RERR
      REAL*8  AppTM(NY1MX),PhsTM(NY1MX)
      COMPLEX*16 HXI0(MMIMX),Ey0(NY1MX),BB0(MMIMX),
     >           Ztm(NY1MX),Hx0(NY1MX)

      mmi = (Ny-1)*(Nzb-1)
      
      CALL SetBound2D_TM(per,Nzb,Ny,Dzb,Dy,CRho,HXI0,HXB)
      CALL MulAibWithXb(Nzb,Ny,BTM,HXB,BB0)
      CALL CopyVectorC16(1,mmi,HXI0,1,mmi,HXI)
      CALL PCG(per,Ny,Nzb,ATM,BB0,HXI,ETol,NtMx,iter,RERR,fwdcond)
      CALL CompEyAtSurface(per,Nzb,Ny,Dzb,Dy,HXI,HXB,CRho,Ey0,Hx0)
      CALL CompRespTM(per,Ny,Ey0,Ztm,AppTM,PhsTM)

      DO is = 1,NSta(im)
        App(is) = DLOG10(AppTM(StaNod(im,is)))
        Phs(is) = PhsTM(StaNod(im,is))*(Pi/D180)
        Zyx(is) = Ztm(StaNod(im,is))
        Hxs(is) = Hx0(StaNod(im,is))
        Eys(is) = Ey0(StaNod(im,is))
      ENDDO


      END ! Fwd2DTM_PCG
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd2DTE_LU(per,Nza,Nz,Ny,Dz,Dy,CRho,ATE,BTE,
     >           im,NSta,StaNod,AII,ipiv,EXI,EXB,
     >           App,Phs,Tipper,Zxy,Exs,Hys,Hzs)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER Nza,Nz,Ny,StaNod(NMODMX,NSTAMX),im,NSta(*)
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      REAL*8  ATE(MMIMX,4),BTE(MMBMX),App(*),Phs(*)
      COMPLEX*16 EXI(*),EXB(*),Zxy(*),Exs(*),Hys(*),Hzs(*)
      COMPLEX*16 Tipper(*)

      INTEGER ipiv(*),mmi,mmb,kl,ku,info,is
      REAL*8  AppTE(NY1MX),PhsTE(NY1MX)
      COMPLEX*16 AII(NZ3MX,MMIMX),EXI0(MMIMX),Hy0(NY1MX),
     >           Zte(NY1MX),Ex0(NY1MX),Hz0(NY1MX)

      mmi = (Ny-1)*(Nz-1)
      mmb = 2*Ny+2*Nz
      kl = Nz-1
      ku = Nz-1

      CALL SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,CRho,EXI0,EXB)
      CALL FormAII(per,Nz,Ny,ATE,AII,ipiv)
      CALL MulAibWithXb(Nz,Ny,BTE,EXB,EXI)
      CALL ZGBTRS('N',mmi,kl,ku,1,AII,NZ3MX,ipiv,EXI,MMIMX,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR WHILE SOLVING FWD !!!'
        WRITE(6,*) '!!! Please, check your model and restart !!!'
        WRITE(6,*) 
     >  '!!! If problem persists, contact wsiripun@oce.orst.edu !!'
        STOP
      ENDIF

      CALL CompHyHzAtSurface(per,Nza,Nz,Ny,Dz,Dy,EXI,EXB,CRho,
     >     Hy0,Ex0,Hz0)
      CALL CompRespTE(per,Ny,Ex0,Hy0,Zte,AppTE,PhsTE)


      DO is = 1,NSta(im)
        App(is) = DLOG10(AppTE(StaNod(im,is)))
        Phs(is) = PhsTE(StaNod(im,is))*(Pi/D180)
        Zxy(is) = Zte(StaNod(im,is))
        Exs(is) = Ex0(StaNod(im,is))
        Hys(is) = Hy0(StaNod(im,is))
        Hzs(is) = Hz0(StaNod(im,is))
        Tipper(is) = Hzs(is)/Hys(is)
      ENDDO

      END ! Fwd2DTE_LU

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Fwd2DTE_PCG(per,Nza,Nz,Ny,Dz,Dy,CRho,ATE,BTE,
     >                       im,NSta,StaNod,ETOL,NtMx,fwdcond,
     >                    EXI,EXB,App,Phs,Tipper,Zxy,Exs,Hys,Hzs)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER Nza,Nz,Ny,StaNod(NMODMX,NSTAMX),im,NSta(*),NtMx,
     >        fwdcond
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX),ETOL
      REAL*8  ATE(MMIMX,4),BTE(MMBMX),App(*),Phs(*)
      COMPLEX*16 EXI(*),EXB(*),Zxy(*),Exs(*),Hys(*),Hzs(*)
      COMPLEX*16 Tipper(*)

      INTEGER mmi,is,iter
      REAL*8  AppTE(NY1MX),PhsTE(NY1MX),RERR
      COMPLEX*16 EXI0(MMIMX),Hy0(NY1MX),BB0(MMIMX),
     >           Zte(NY1MX),Ex0(NY1MX),Hz0(NY1MX)

      mmi = (Ny-1)*(Nz-1)
      CALL SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,CRho,EXI0,EXB)
      CALL MulAibWithXb(Nz,Ny,BTE,EXB,BB0)
      CALL CopyVectorC16(1,mmi,EXI0,1,mmi,EXI)
      CALL PCG(per,Ny,Nz,ATE,BB0,EXI,ETol,NtMx,iter,RERR,fwdcond)
      CALL CompHyHzAtSurface(per,Nza,Nz,Ny,Dz,Dy,EXI,EXB,CRho,
     >     Hy0,Ex0,Hz0)
      CALL CompRespTE(per,Ny,Ex0,Hy0,Zte,AppTE,PhsTE)

      DO is = 1,NSta(im)
        App(is) = DLOG10(AppTE(StaNod(im,is)))
        Phs(is) = PhsTE(StaNod(im,is))*(Pi/D180)
        Zxy(is) = Zte(StaNod(im,is))
        Exs(is) = Ex0(StaNod(im,is))
        Hys(is) = Hy0(StaNod(im,is))
        Hzs(is) = Hz0(StaNod(im,is))
        Tipper(is) = Hzs(is)/Hys(is)
      ENDDO

      END ! Fwd2DTE_PCG

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE FormAII(per,Nz0,Ny,AA,AII,ipiv)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER Nz0,Ny,ipiv(*)
      REAL*8  per,AA(MMIMX,4)
      COMPLEX*16 AII(NZ3MX,MMIMX)

      INTEGER mmi,jj,nz3,kl,ku,kc,info
      REAL*8  Omega,Omue

      Omega = (D2*Pi)/per
      Omue  = Omega*Mue

C     Forming matrix Aii in LAPACK's format
      nz3 = 3*(Nz0-1) + 1
      mmi = (Ny-1)*(Nz0-1)
      CALL ConstantMatrixC16(AII,NZ3MX,MMIMX,nz3,mmi,D0)

      kl = Nz0-1
      ku = Nz0-1
      kc = kl + ku + 1

C     diagonal term
      DO jj = 1,mmi
        AII(kc,jj) = DCMPLX(AA(jj,1),Omue*AA(jj,4))
      ENDDO

C     first lower and upper diagonal terms
      DO jj = 1,mmi-1
        AII(kc-1,jj+1) = AA(jj,2)
        AII(kc+1,jj)   = AA(jj,2)
      ENDDO

C     second lower and upper diagonal terms
      DO jj = 1,mmi-(Nz0-1)
        AII(kc-(Nz0-1),jj+(Nz0-1)) = AA(jj,3)
        AII(kc+(Nz0-1),jj)         = AA(jj,3)
      ENDDO
  
C     LU Decompose 
      CALL ZGBTRF(mmi,mmi,kl,ku,AII,NZ3MX,ipiv,info)
      IF (info.NE.0) THEN
        WRITE(6,*) 
     >  '!!! ATTENTION, ERROR WHILE DECOMPOSING MATRIX !!!'
        WRITE(6,*) '!!! Please, check your model and restart !!!'
        WRITE(6,*) 
     >  '!!! If problem persists, contact wsiripun@oce.orst.edu !!'
        STOP
      ENDIF

c100   FORMAT(7E11.3)
      
      RETURN
      END  ! FormAII()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE MulAibWithXb(Nz0,Ny,Aib,Xb,AXb)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nz0,Ny
      REAL*8  Aib(*)
      COMPLEX*16 Xb(*),AXb(*)

      INTEGER mmi,mmb,iy,iz,jj,izz

      mmi = (Ny-1)*(Nz0-1)
      mmb = 2*Ny + 2*Nz0

      CALL ConstantVectorC16(AXb,mmi,D0)

C     left side
      iy = 2
      DO iz = 2,Nz0
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = - Aib(iz)*Xb(iz)
      ENDDO

C     right side
      iy = Ny
      izz = 2*Ny + Nz0 + 1
      DO iz = 2,Nz0
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

C     top side
      izz = Nz0 + 2
      iz  = 2
      DO iy = 2,Ny
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = AXb(jj) - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

C     bottom side
      izz = Nz0 + Ny + 1
      iz  = Nz0
      DO iy = 2,Ny
        jj = (iy-2)*(Nz0-1) + iz-1
        AXb(jj) = AXb(jj) - Aib(izz)*Xb(izz)
        izz = izz + 1
      ENDDO

      END ! MulAibWithXb

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE PCG(per,Ny,Nz0,AI0,B0,XX,ETol,NtMx,iter,RERR,
     >           fwdcond)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nz0,NtMx,iter,fwdcond
      REAL*8  AI0(MMIMX,4),ETol,per,RERR
      COMPLEX*16 XX(*),B0(*)

      INTEGER jj,mmi,it
      REAL*8  norm,dnorm
      COMPLEX*16 R1,R2,R3,alpha,beta
      COMPLEX*16 AX(MMIMX),RR(MMIMX),DD(MMIMX),AD(MMIMX)
      COMPLEX*16 ZR(MMIMX),ID(MMIMX)
      COMPLEX*16 ilu(9,MMIMX)

      COMPLEX*16 ZDOTU,ZDOTC

      fwdcond = 0
      mmi  = (Ny-1)*(Nz0-1)
      CALL ConstantVectorC16(ID,mmi,D0)
      CALL ConstantVectorC16(RR,mmi,D0)
      CALL ConstantVectorC16(AX,mmi,D0)
      CALL ConstantVectorC16(DD,mmi,D0)
      CALL ConstantVectorC16(AD,mmi,D0)
      CALL ConstantVectorC16(ZR,mmi,D0)

C     compute ILU decomposition
      CALL ILU3(per,Nz0,Ny,AI0,ilu)

C     compute norm
      dnorm = CDSQRT(ZDOTC(mmi,B0,1,B0,1))

      CALL MulAX(per,Ny,Nz0,AI0,XX,AX)
      DO jj = 1,mmi
        RR(jj) = B0(jj)-AX(jj)
      ENDDO
      norm = CDSQRT(ZDOTC(mmi,RR,1,RR,1))

      it = 0
      RERR = norm/dnorm

      IF (RERR.LE.ETol) THEN
        GOTO 2000
      ELSE
        GOTO 1000
      ENDIF
1000  CONTINUE

      it = it+1
      CALL Solve_ILU3(Ny,Nz0,ilu,RR,ZR)
      R1 = ZDOTU(mmi,RR,1,ZR,1)
      IF (it.EQ.1) THEN
        CALL ZCOPY(mmi,ZR,1,DD,1)
      ELSE
        beta = R1/R3
        CALL ZAXPY(mmi,D1/beta,ZR,1,DD,1)
        CALL ZSCAL(mmi,beta,DD,1)
      ENDIF
     
      CALL MulAX(per,Ny,Nz0,AI0,DD,AD)
      R2 = ZDOTU(mmi,DD,1,AD,1)
      alpha = R1/R2
      CALL ZAXPY(mmi,alpha,DD,1,XX,1)
      CALL ZAXPY(mmi,-alpha,AD,1,RR,1)
      R3 = R1

      norm = CDSQRT(ZDOTC(mmi,RR,1,RR,1))
      RERR = norm/dnorm
    
      IF (RERR.LE.ETol) THEN
        fwdcond = 0
        GOTO 2000
      ELSE
        IF (it.GE.NtMx) THEN
          fwdcond = 1
          WRITE(6,3050)
          WRITE(6,3060) 
          WRITE(99,3050)
          WRITE(99,3060)
          GOTO 2000
        ENDIF
        GOTO 1000
      ENDIF
2000  CONTINUE
      iter = it

3050  FORMAT('MAX NUMBER OF <PCG> ITERATION IS REACHED ')
3060  FORMAT('HIGHER VALUE OF <LGM> WILL BE REPLACED ')

      RETURN
      END ! PCG()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     IL(i,i)    = 1.0
C     IL(i+1,i)  = ilu2(i+1)
C     IL(i+n-2,i)= ilu9(i+n-1)
C     IL(i+n-1,i)= ilu7(i+n-1)
C     IL(i+n,i)  = ilu1(i+n)
C     IL(i,j)    = 0  otherwise

C     IU(i,i)    = ilu3(i)
C     IU(i,i+1)  = ilu4(i)
C     IU(i,i+n-2)= ilu8(i)
C     IU(i,i+n-1)= ilu6(i)
C     IU(i,i+n)  = ilu5(i)
C     IU(i,j)    = 0  otherwise


      SUBROUTINE ILU3(per,Nz0,Ny,AI,ilu)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nz0,Ny
      REAL*8  AI(MMIMX,4),per
      COMPLEX*16 ilu(9,MMIMX)

      INTEGER ii,jj,kk,pp,qq,nzz,mmi
      REAL*8  omega,omue

      mmi = (Ny-1)*(Nz0-1)
      nzz = (Nz0-1)
      omega = (D2*Pi)/per
      omue  = omega*Mue

      CALL ConstantMatrixC16(ilu,9,MMIMX,9,mmi,D0)

      DO jj = 1,mmi
        ilu(3,jj)  = DCMPLX(AI(jj,1),omue*AI(jj,4))
      ENDDO ! jj
      DO jj = 1,mmi-1
        ilu(4,jj)  = AI(jj,2)
      ENDDO ! jj
      DO jj = 1,mmi-nzz
        ilu(5,jj)  = AI(jj,3)
      ENDDO ! jj

      DO jj = 2,mmi
        ilu(2,jj) = ilu(4,jj-1)
      ENDDO ! jj
      DO jj = nzz+1,mmi
        ilu(1,jj) = ilu(5,jj-nzz)
      ENDDO ! jj

      ilu(9,nzz-1) = D0
      ilu(7,nzz)   = D0
      ilu(8,1)     = D0
      ilu(6,1)     = D0

      DO ii = 1,mmi-nzz
        jj = ii+1
        qq = ii+nzz-2
        kk = ii+nzz-1
        pp = ii+nzz

        ilu(2,jj) = ilu(2,jj)/ilu(3,ii)
        ilu(9,qq) = ilu(9,qq)/ilu(3,ii)
        ilu(7,kk) = ilu(7,kk)/ilu(3,ii)
        ilu(1,pp) = ilu(1,pp)/ilu(3,ii)

        ilu(3,jj) = ilu(3,jj) - ilu(4,ii)*ilu(2,jj)
        ilu(8,jj) =           - ilu(6,ii)*ilu(2,jj)
        ilu(6,jj) =           - ilu(5,ii)*ilu(2,jj)

        ilu(3,qq) = ilu(3,qq) - ilu(8,ii)*ilu(9,qq)
        ilu(4,qq) = ilu(4,qq) - ilu(6,ii)*ilu(9,qq)

        ilu(9,kk) =           - ilu(4,ii)*ilu(7,kk)
        ilu(2,kk) = ilu(2,kk) - ilu(8,ii)*ilu(7,kk)
        ilu(3,kk) = ilu(3,kk) - ilu(6,ii)*ilu(7,kk)
        ilu(4,kk) = ilu(4,kk) - ilu(5,ii)*ilu(7,kk)

        ilu(7,pp) =           - ilu(4,ii)*ilu(1,pp)
        ilu(2,pp) = ilu(2,pp) - ilu(6,ii)*ilu(1,pp)
        ilu(3,pp) = ilu(3,pp) - ilu(5,ii)*ilu(1,pp)
      ENDDO

      ii = mmi-nzz+1
      jj = ii+1
      qq = ii+nzz-2
      kk = ii+nzz-1

      ilu(2,jj) = ilu(2,jj)/ilu(3,ii)
      ilu(9,qq) = ilu(9,qq)/ilu(3,ii)
      ilu(7,kk) = ilu(7,kk)/ilu(3,ii)
 
      ilu(3,jj) = ilu(3,jj) - ilu(4,ii)*ilu(2,jj)
      ilu(8,jj) =           - ilu(6,ii)*ilu(2,jj)
 
      ilu(3,qq) = ilu(3,qq) - ilu(8,ii)*ilu(9,qq)
      ilu(4,qq) = ilu(4,qq) - ilu(6,ii)*ilu(9,qq)

      ilu(9,kk) =           - ilu(4,ii)*ilu(7,kk)
      ilu(2,kk) = ilu(2,kk) - ilu(8,ii)*ilu(7,kk)
      ilu(3,kk) = ilu(3,kk) - ilu(6,ii)*ilu(7,kk)

      ii = mmi-nzz+2
      jj = ii+1
      qq = ii+nzz-2

      ilu(2,jj) = ilu(2,jj)/ilu(3,ii)
      ilu(9,qq) = ilu(9,qq)/ilu(3,ii)

      ilu(3,jj) = ilu(3,jj) - ilu(4,ii)*ilu(2,jj)

      ilu(3,qq) = ilu(3,qq) - ilu(8,ii)*ilu(9,qq)

      DO ii = mmi-nzz+3,mmi-1
        jj = ii+1
        ilu(2,jj) = ilu(2,jj)/ilu(3,ii)

        ilu(3,jj) = ilu(3,jj) - ilu(4,ii)*ilu(2,jj)
      ENDDO

      RETURN
      END ! ILU

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
     
      SUBROUTINE MulAX(per,Ny,Nz0,AI,XI,YI)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nz0
      REAL*8  AI(MMIMX,4),per
      COMPLEX*16 XI(*),YI(*)

      INTEGER mmi,jj
      REAL*8  omega,omue

      omega = (D2*Pi)/per
      omue  = omega*Mue
      mmi = (Ny-1)*(Nz0-1)
      CALL ConstantVectorC16(YI,mmi,D0)

      CALL ZCOPY(mmi,XI,1,YI,1)
      DO jj = 1,mmi
        YI(jj) = DCMPLX(AI(jj,1),omue*AI(jj,4))*XI(jj)
      ENDDO
      DO jj = 1,mmi-1
        YI(jj) = YI(jj) + AI(jj,2)*XI(jj+1)
      ENDDO
      DO jj = 1,mmi-(Nz0-1)
        YI(jj) = YI(jj) + AI(jj,3)*XI(jj+(Nz0-1))
      ENDDO

      DO jj = 2,mmi
        YI(jj) = YI(jj) + AI(jj-1,2)*XI(jj-1)
      ENDDO
      DO jj = (Nz0-1)+1,mmi
        YI(jj) = YI(jj) + AI(jj-(Nz0-1),3)*XI(jj-(Nz0-1))
      ENDDO

      RETURN
      END ! MulAX

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Solve_ILU3(Ny,Nz0,ilu,RR,ZR)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny,Nz0
      COMPLEX*16 ilu(9,MMIMX),RR(*),ZR(*)

      INTEGER mmi,nzz,jj
      COMPLEX*16 YR(MMIMX)

C     IL(i,i)    = 1.0
C     IL(i+1,i)  = ilu2(i+1)
C     IL(i+n-2,i)= ilu9(i+n-1)
C     IL(i+n-1,i)= ilu7(i+n-1)
C     IL(i+n,i)  = ilu1(i+n)
C     IL(i,j)    = 0  otherwise
C
C     IU(i,i)    = ilu3(i)
C     IU(i,i+1)  = ilu4(i)
C     IU(i,i+n-2)= ilu8(i)
C     IU(i,i+n-1)= ilu6(i)
C     IU(i,i+n)  = ilu5(i)
C     IU(i,j)    = 0  otherwise


      mmi = (Ny-1)*(Nz0-1)
      nzz = Nz0-1

C     Forward substitution
      CALL ConstantVectorC16(YR,mmi,D0)

      YR(1) = RR(1)
      DO jj = 2,nzz-3
        YR(jj) = RR(jj) - ilu(2,jj)*YR(jj-1)
      ENDDO
      jj = nzz-2
      YR(jj) = RR(jj) - ilu(2,jj)*YR(jj-1)-ilu(9,jj)*YR(jj-(nzz-2))
      jj = nzz-1
      YR(jj) = RR(jj) - ilu(2,jj)*YR(jj-1)-ilu(9,jj)*YR(jj-(nzz-2))
     >                - ilu(7,jj)*YR(jj-(nzz-1))
      DO jj = nzz,mmi
        YR(jj) = RR(jj)-ilu(2,jj)*YR(jj-1)-ilu(9,jj)*YR(jj-(nzz-2))
     >         - ilu(7,jj)*YR(jj-(nzz-1)) -ilu(1,jj)*YR(jj-(nzz))
      ENDDO


C     Backward subsitution
      ZR(mmi) = YR(mmi)/ilu(3,mmi)
      DO jj = mmi-1,mmi-(nzz-3),-1
        ZR(jj) = (YR(jj) - ilu(4,jj)*ZR(jj+1))/ilu(3,jj) 
      ENDDO
      jj = mmi-(nzz-2)
      ZR(jj) = (YR(jj) - ilu(4,jj)*ZR(jj+1) - ilu(8,jj)*ZR(jj+(nzz-2)))/
     >          ilu(3,jj) 
      jj = mmi-(nzz-1)
      ZR(jj) = (YR(jj) - ilu(4,jj)*ZR(jj+1) - ilu(8,jj)*ZR(jj+(nzz-2))
     >         -ilu(6,jj)*ZR(jj+(nzz-1)))/ilu(3,jj) 
      DO jj = mmi-nzz,1,-1
      ZR(jj) = (YR(jj) - ilu(4,jj)*ZR(jj+1) - ilu(8,jj)*ZR(jj+(nzz-2))
     >        -ilu(6,jj)*ZR(jj+(nzz-1)) -ilu(5,jj)*ZR(jj+nzz))/ilu(3,jj)
      ENDDO


      END ! Solve_ILU3()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetBound2D_TM(per,Nzb,Ny,Dzb,Dy,CRho,HXI0,HXB)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      REAL*8  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 HXI0(*),HXB(*)

      INTEGER nzb1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      REAL*8  r1d(NZ1MX)
      COMPLEX*16 X1D(NZ2MX)

      nzb1 = Nzb + 1
      ny1  = Ny  + 1

      mmb = 2*Nzb + 2*Ny
      mmi = (Nzb-1)*(Ny-1)
   
      CALL ConstantVectorC16(HXB,mmb,D0)
      CALL ConstantVectorC16(HXI0,mmi,D0)

      np1 = NZ0MX
      np2 = NY0MX

C     left fields
      CALL TransMatrixToVectorR8(CRho,np1,np2,1,1,Nzb,1,r1d,1,Nzb)
      CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      CALL CopyVectorC16(1,Nzb+1,X1D,1,Nzb+1,HXB)

C     surface fields
      DO ii = Nzb+2,Nzb+Ny
        HXB(ii) = D1
      ENDDO


C     Interior and bottom fields
      DO iy = 2,Ny
        DO iz = 1,Nzb
          r1d(iz) = (CRho(iz,iy-1)*Dy(iy-1)+CRho(iz,iy)*Dy(iy))/
     >              (Dy(iy-1) + Dy(iy))
        ENDDO
        CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
        HXB(Nzb+Ny-1+iy) = X1D(Nzb+1)
        DO iz = 2,Nzb
          jj = (iy-2)*(Nzb-1) + iz-1
          HXI0(jj) = X1D(iz)
        ENDDO
      ENDDO 

C     right fields
      CALL TransMatrixToVectorR8(CRho,np1,np2,1,Ny,Nzb,Ny,r1d,1,Nzb)
      CALL Fwd1D_TM(per,Nzb,Dzb,r1d,X1D)
      CALL CopyVectorC16(1,Nzb+1,X1D,Nzb+2*Ny,2*Nzb+2*Ny,HXB)
      

c100   FORMAT(i5,2e11.3)

      RETURN
      END ! SetBound2D_TM

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SetBound2D_TE(per,Nza,Nz,Ny,Dz,Dy,CRho,EXI0,EXB)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz,Ny
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 EXI0(*),EXB(*)

      INTEGER nz1,ny1,mmi,mmb,np1,np2,iz,iy,ii,jj
      REAL*8  s1d(NZ1MX),CCon(NZ0MX,NY0MX)
      COMPLEX*16 X1D(NZ2MX)


      nz1 = Nz + 1
      ny1 = Ny + 1

      mmb = 2*Nz + 2*Ny
      mmi = (Nz-1)*(Ny-1)

      DO iy = 1,Ny
        DO iz = 1,Nz-Nza
          CCon(iz,iy) = 1./CRho(iz,iy)
        ENDDO ! iz
      ENDDO ! iy
   
      CALL ConstantVectorC16(EXB,mmb,D0)
      CALL ConstantVectorC16(EXI0,mmi,D0)

      np1 = NZ0MX
      np2 = NY0MX

C     left fields
      CALL TransMatrixToVectorR8(CCon,np1,np2,1,1,Nz-Nza,1,
     >                           s1d,1,Nz-Nza)
      CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      CALL CopyVectorC16(1,Nz+1,X1D,1,Nz+1,EXB)


C     surface fields
      DO ii = Nz+2,Nz+Ny
        EXB(ii) = D1
      ENDDO


C     Interior and bottom fields
      DO iy = 2,Ny
        DO iz = 1,Nz-Nza
          s1d(iz) = (CCon(iz,iy-1)*Dy(iy-1)+CCon(iz,iy)*Dy(iy))/
     >              (Dy(iy-1) + Dy(iy))
        ENDDO
        CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
        EXB(Nz+Ny-1+iy) = X1D(Nz+1)
        DO iz = 2,Nz
          jj = (iy-2)*(Nz-1) + iz-1
          EXI0(jj) = X1D(iz)
        ENDDO
      ENDDO 

C     right fields
      CALL TransMatrixToVectorR8(CCon,np1,np2,1,Ny,Nz-Nza,Ny,
     >                           s1d,1,Nz-Nza)
      CALL Fwd1D_TE(per,Nza,Nz,Dz,s1d,X1D)
      CALL CopyVectorC16(1,Nz+1,X1D,Nz+2*Ny,2*Nz+2*Ny,EXB)

c100   FORMAT(i5,2e11.3)

      RETURN
      END ! SetBound2D_TE


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompEyAtSurface(per,Nzb,Ny,Dzb,Dy,HXI,HXB,CRho,
     >                           Ey0,Hx0)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      REAL*8  per,Dzb(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 HXI(*),HXB(*),Ey0(*),Hx0(*)

      INTEGER iz,iy,jj,jj1,jj2
      REAL*8  omega,omue,dz1,dz2,dzz,dyy,ryy,rzr,rzl
      COMPLEX*16 ai00,ai10,ai01,ab00

      omega = (D2*Pi)/per
      omue  = omega*Mue

      iz  = 2
      dz1 = Dzb(1)
      dz2 = Dzb(2)
      dzz = dz1 + dz2
      DO iy = 2,Ny
        dyy  = Dy(iy) + Dy(iy-1)
        ryy  = (CRho(1,iy)*Dy(iy)+CRho(1,iy-1)*Dy(iy-1))/dyy
        rzr  = (CRho(1,iy)*dz1   +CRho(2,iy)*dz2)/dzz
        rzl  = (CRho(1,iy-1)*dz1 +CRho(2,iy-1)*dz2)/dzz
        ai00 = ryy/dz1 + DCMPLX(D0,omue*dz1/D8) -
     >         rzr*dz1/(D4*Dy(iy)*dyy) + rzl*dz1/(D4*Dy(iy-1)*dyy)
        ai10 = rzr*dz1/(D4*Dy(iy)*dyy)
        ai01 = -rzl*dz1/(D4*Dy(iy-1)*dyy)
        ab00 = -ryy/dz1 + DCMPLX(D0,D3*omue*dz1/D8)
        IF (iy.EQ.2) THEN
          jj  = (iy-2)*(Nzb-1) + iz - 1
          jj1 = (iy+1-2)*(Nzb-1) + iz - 1
          Ey0(iy) = ai00*HXI(jj) + ai10*HXI(jj1) +
     >              ai01*HXB(2)  + ab00*HXB(Nzb+iy)
        ENDIF
        IF (iy.EQ.Ny) THEN
          jj  = (iy-2)*(Nzb-1) + iz - 1
          jj2 = (iy-1-2)*(Nzb-1) + iz - 1
          Ey0(iy) = ai00*HXI(jj)  + ai10*HXB(Nzb+2*Ny+1) +
     >              ai01*HXI(jj2) + ab00*HXB(Nzb+iy)
        ENDIF
        IF ((iy.GT.2).AND.(iy.LT.Ny)) THEN
          jj  = (iy-2)*(Nzb-1) + iz - 1
          jj1 = (iy+1-2)*(Nzb-1) + iz - 1
          jj2 = (iy-1-2)*(Nzb-1) + iz - 1
          Ey0(iy) = ai00*HXI(jj)  + ai10*HXI(jj1) + 
     >              ai01*HXI(jj2) + ab00*HXB(Nzb+iy)
        ENDIF
      ENDDO ! iy
      ryy = CRho(1,1)
      ai00 = ryy/dz1 + DCMPLX(D0,omue*dz1/D8) 
      ab00 = -ryy/dz1 + DCMPLX(D0,D3*omue*dz1/D8)
      Ey0(1) = ai00*HXB(2) + ab00*HXB(1)

      ryy = CRho(1,Ny)
      ai00 = ryy/dz1 + DCMPLX(D0,omue*dz1/D8)
      ab00 = -ryy/dz1 + DCMPLX(D0,D3*omue*dz1/D8)
      Ey0(Ny+1) = ai00*HXB(2*Ny+Nzb+1) + ab00*HXB(2*Ny+Nzb)

      DO iy = 1,Ny+1
        Hx0(iy) = D1
      ENDDO ! iy

      RETURN
      END ! CompEyAtSurface
   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
     
      SUBROUTINE CompRespTM(per,Ny,Eys,Zyx,AppTM,PhsTM)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny
      REAL*8  per,AppTM(*),PhsTM(*)
      COMPLEX*16 Eys(*),Zyx(*)
      
      INTEGER iy
      REAL*8  omega,omue
  
      omega = (D2*Pi)/per
      omue  = omega*Mue

      DO iy = 1,Ny+1
        Zyx(iy)   = -Eys(iy)
        AppTM(iy) = CDABS(Zyx(iy)*Zyx(iy))/omue 
        PhsTM(iy) = DATAN2(DIMAG(Zyx(iy)),DREAL(Zyx(iy)))*(D180/Pi)
      ENDDO
      
      RETURN
      END ! CompRespTM

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompHyHzAtSurface(per,Nza,Nz,Ny,Dz,Dy,EXI,EXB,CRho,
     >                             Hy0,Ex0,Hz0)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nz,Ny
      REAL*8  per,Dz(*),Dy(*),CRho(NZ0MX,NY0MX)
      COMPLEX*16 EXI(*),EXB(*),Hy0(*),Ex0(*),Hz0(*)

      INTEGER iz,iy,j00,j0r,j0l,j10,j1r,j1l,izz
      REAL*8  omega,omue,dz1,dz2,dzz,dyy,sav
      COMPLEX*16 a00,a0r,a0l,a10,a1r,a1l,dyz

      omega = (D2*Pi)/per
      omue  = omega*Mue

      dz1 = Dz(Nza+1)
      dz2 = Dz(Nza+2)
      dzz = dz1 + dz2
      iz  = Nza + 1
      DO iy = 2,Ny
        dyy  = Dy(iy) + Dy(iy-1)
        sav  = (Dy(iy)/CRho(1,iy) + Dy(iy-1)/CRho(1,iy-1))/dyy
        dyz  = DCMPLX(D0,dz1/(dyy*omue))

        a00 = DCMPLX(D0,D1/(omue*dz1)) + (D3/D8)*sav*dz1
     >      + (D3/D4)*dyz*(D1/Dy(iy) + D1/Dy(iy-1))
        a0r = -(D3/D4)*dyz*(D1/Dy(iy))
        a0l = -(D3/D4)*dyz*(D1/Dy(iy-1))
        a10 = -DCMPLX(D0,D1/(omue*dz1)) + (D1/D8)*sav*dz1
     >      + (D1/D4)*dyz*(D1/Dy(iy) + D1/Dy(iy-1))
        a1r = -(D1/D4)*dyz*(D1/Dy(iy))
        a1l = -(D1/D4)*dyz*(D1/Dy(iy-1))
    
        IF (iy.EQ.2) THEN
          j00 = (iy-2)*(Nz-1)   + iz - 1
          j0r = (iy+1-2)*(Nz-1) + iz - 1
          j10 = (iy-2)*(Nz-1)   + (iz+1) - 1
          j1r = (iy+1-2)*(Nz-1) + (iz+1) - 1
          izz = Nza
          Hy0(iy) = a00*EXI(j00) + a0l*EXB(izz+1) + a0r*EXI(j0r) + 
     >              a10*EXI(j10) + a1l*EXB(izz+2) + a1r*EXI(j1r)
        ENDIF
        IF (iy.EQ.Ny) THEN
          j00 = (iy-2)*(Nz-1)   + iz - 1
          j0l = (iy-1-2)*(Nz-1) + iz - 1
          j10 = (iy-2)*(Nz-1)   + (iz+1) - 1
          j1l = (iy-1-2)*(Nz-1) + (iz+1) - 1
          izz = (Nz+1)+2*(Ny-1)+Nza
          Hy0(iy) = a00*EXI(j00) + a0l*EXI(j0l) + a0r*EXB(izz+1) + 
     >              a10*EXI(j10) + a1l*EXI(j1l) + a1r*EXI(izz+2)
        ENDIF
        IF ((iy.GT.2).AND.(iy.LT.Ny)) THEN
          j00 = (iy-2)*(Nz-1)   + iz - 1
          j0l = (iy-1-2)*(Nz-1) + iz - 1
          j0r = (iy+1-2)*(Nz-1) + iz - 1
          j10 = (iy-2)*(Nz-1)   + (iz+1) - 1
          j1r = (iy+1-2)*(Nz-1) + (iz+1) - 1
          j1l = (iy-1-2)*(Nz-1) + (iz+1) - 1
          Hy0(iy) = a00*EXI(j00) + a0l*EXI(j0l) + a0r*EXI(j0r) + 
     >              a10*EXI(j10) + a1l*EXI(j1l) + a1r*EXI(j1r)
        ENDIF
      ENDDO ! iy

      izz = Nza
      sav = D1/CRho(1,1)
      a00 =  DCMPLX(D0,D1/(omue*dz1)) + (D3/D8)*sav*dz1
      a10 = -DCMPLX(D0,D1/(omue*dz1)) + (D1/D8)*sav*dz1
      Hy0(1) = a00*EXB(izz+1) + a10*EXB(izz+2) 

      izz = (Nz+1)+2*(Ny-1)+Nza
      sav = D1/CRho(1,Ny)
      a00 =  DCMPLX(D0,D1/(omue*dz1)) + (D3/D8)*sav*dz1
      a10 = -DCMPLX(D0,D1/(omue*dz1)) + (D1/D8)*sav*dz1
      Hy0(Ny+1) = a00*EXB(izz+1) + a10*EXB(izz+2) 

      iz  = Nza + 1
      DO iy = 2,Ny
        j00 = (iy-2)*(Nz-1)   + iz - 1
        Ex0(iy) = EXI(j00)
      ENDDO ! iy
      Ex0(1)    = EXB(Nza+1)
      Ex0(Ny+1) = EXB((Nz+1)+2*(Ny-1)+(Nza+1)) 

      DO iy = 2,Ny-1
        Hz0(iy) = DCMPLX(D0,D1/omue)*(Ex0(iy+1)-Ex0(iy-1))/
     >            (Dy(iy)+Dy(iy-1))
      ENDDO ! iy
      Hz0(1)    = DCMPLX(D0,D1/omue)*(Ex0(2)-Ex0(1))/Dy(1)
      Hz0(Ny+1) = DCMPLX(D0,D1/omue)*(Ex0(Ny+1)-Ex0(Ny))/Dy(Ny)

      RETURN
      END ! CompHyAtSurface
   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
     
      SUBROUTINE CompRespTE(per,Ny,Exs,Hys,Zxy,AppTE,PhsTE)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Ny
      REAL*8  per,AppTE(*),PhsTE(*)
      COMPLEX*16 Exs(*),Hys(*),Zxy(*)
      
      INTEGER iy
      REAL*8  omega,omue
  
      omega = (D2*Pi)/per
      omue  = omega*Mue

      DO iy = 1,Ny+1
        Zxy(iy)   = Exs(iy)/Hys(iy)
        AppTE(iy) = CDABS(Zxy(iy)*Zxy(iy))/omue 
        PhsTE(iy) = DATAN2(DIMAG(Zxy(iy)),DREAL(Zxy(iy)))*(D180/Pi)
      ENDDO

      RETURN
      END ! CompRespTE

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

