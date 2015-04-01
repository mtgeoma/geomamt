      SUBROUTINE G03ADZ(N,NX,NY,CVX,LDCVX,CVY,LDCVY,QX,QY,LQY,RDF,TOL,
     *                  IRANKX,IRANKY,NCV,E,WK,IERROR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RDF, TOL
      INTEGER           IERROR, IRANKX, IRANKY, LDCVX, LDCVY, LQY, N,
     *                  NCV, NX, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  CVX(LDCVX,*), CVY(LDCVY,*), E(NY), QX(N,*),
     *                  QY(LQY), WK(NX+NY)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, TEMP
      INTEGER           I, IFAULT, IWK, J
      LOGICAL           SVDX, SVDY
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1,1), WKSP2(1,1)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ
      INTEGER           F06KLF
      EXTERNAL          F02WDZ, F06KLF
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F01QDF, F01QEF, F02WEF, F02WUF, F06FCF,
     *                  F06QFF, DCOPY, DGEMV, DSCAL, DTRSV
C     .. Executable Statements ..
C
C       Asssume NX GE NY
C
      IERROR = 0
      IWK = NX*NX + 1
      IFAULT = 1
      CALL F01QCF(N,NX,QX,N,WK,IFAULT)
      IFAULT = 1
      CALL F01QCF(N,NY,QY,N,WK(NX+1),IFAULT)
C       Copy Ry into  CVY
      CALL F06QFF('U',NY,NY,QY,N,CVY,LDCVY)
C       Put Qx'Qy in  CVX
      IFAULT = 1
      CALL F01QEF('S',N,NY,NY,QY,N,WK(NX+1),CVX,IFAULT)
      IFAULT = 1
      CALL F01QDF('T','S',N,NX,QX,N,WK,NY,QY,N,CVX,IFAULT)
      CALL F06QFF('G',NX,NY,QY,N,CVX,LDCVX)
C       Copy Ry back into lower part of X workspace
      CALL F06QFF('U',NY,NY,CVY,LDCVY,QX(NX+1,1),N)
      COND = F02WDZ(NX,QX,N,WK)
      IF (COND*TOL.GT.1.0D0) THEN
         SVDX = .TRUE.
         IFAULT = 1
         CALL F02WUF(NX,QX,N,0,WKSP1,1,.TRUE.,QY,NX,WK,.TRUE.,QY(IWK),
     *               IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 1
            RETURN
         END IF
         IRANKX = F06KLF(NX,WK,1,TOL)
         IF (IRANKX.LE.0) THEN
            IERROR = 2
            RETURN
         END IF
         DO 20 I = 1, IRANKX
            WK(I) = 1.0D0/WK(I)
   20    CONTINUE
         DO 40 I = 1, NX
            CALL F06FCF(IRANKX,WK,1,QX(1,I),1)
   40    CONTINUE
         DO 60 I = 1, NY
            CALL DCOPY(NX,CVX(1,I),1,WK,1)
            CALL DGEMV('T',NX,IRANKX,1.0D0,QY,NX,WK,1,0.0D0,CVX(1,I),1)
   60    CONTINUE
      ELSE
         SVDX = .FALSE.
         IRANKX = NX
      END IF
      COND = F02WDZ(NY,QX(NX+1,1),N,WK)
      IF (COND*TOL.GT.1.0D0) THEN
         SVDY = .TRUE.
         IFAULT = 1
         CALL F02WUF(NY,QX(NX+1,1),N,0,WKSP1,1,.TRUE.,QY,NY,WK,.TRUE.,
     *               QY(IWK),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 1
            RETURN
         END IF
         IRANKY = F06KLF(NY,WK,1,TOL)
         IF (IRANKY.LE.0) THEN
            IERROR = 3
            RETURN
         END IF
         DO 80 I = 1, IRANKY
            WK(I) = 1.0D0/WK(I)
   80    CONTINUE
         DO 100 I = 1, NY
            CALL F06FCF(IRANKY,WK,1,QX(NX+1,I),1)
  100    CONTINUE
         DO 120 I = 1, IRANKY
            CALL DGEMV('N',IRANKX,NY,1.0D0,CVX,LDCVX,QY((I-1)*NY+1),1,
     *                 0.0D0,QY(NY*NY+1+(I-1)*IRANKX),1)
  120    CONTINUE
         CALL F06QFF('G',IRANKX,IRANKY,QY(NY*NY+1),IRANKX,CVX,LDCVX)
      ELSE
         SVDY = .FALSE.
         IRANKY = NY
      END IF
      IF (IRANKX.GE.IRANKY) THEN
         IFAULT = 1
         CALL F02WEF(IRANKX,IRANKY,CVX,LDCVX,0,WKSP1,1,.TRUE.,WKSP2,1,E,
     *               .TRUE.,CVY,LDCVY,QY,IFAULT)
         IF (IFAULT.GT.0) THEN
            IERROR = 1
            RETURN
         END IF
         NCV = F06KLF(IRANKY,E,1,TOL)
      ELSE
         CALL F06QFF('G',IRANKX,IRANKY,CVX,LDCVX,CVY,LDCVY)
         IFAULT = 1
         CALL F02WEF(IRANKX,IRANKY,CVY,LDCVY,0,WKSP1,1,.TRUE.,CVX,LDCVX,
     *               E,.TRUE.,WKSP2,1,QY,IFAULT)
         IF (IFAULT.GT.0) THEN
            IERROR = 1
            RETURN
         END IF
         NCV = F06KLF(IRANKX,E,1,TOL)
      END IF
      IF ( .NOT. SVDX) THEN
         DO 140 I = 1, NCV
            CALL DTRSV('U','N','N',NX,QX,N,CVX(1,I),1)
            CALL DSCAL(NX,RDF,CVX(1,I),1)
  140    CONTINUE
      ELSE
         DO 160 I = 1, NCV
            CALL DCOPY(IRANKX,CVX(1,I),1,WK,1)
            CALL DGEMV('T',IRANKX,NX,RDF,QX,N,WK,1,0.0D0,CVX(1,I),1)
  160    CONTINUE
      END IF
      DO 200 I = 1, NY
         DO 180 J = 1, I - 1
            TEMP = CVY(J,I)
            CVY(J,I) = CVY(I,J)
            CVY(I,J) = TEMP
  180    CONTINUE
  200 CONTINUE
      IF ( .NOT. SVDY) THEN
         DO 220 I = 1, NCV
            CALL DTRSV('U','N','N',NY,QX(NX+1,1),N,CVY(1,I),1)
            CALL DSCAL(NY,RDF,CVY(1,I),1)
  220    CONTINUE
      ELSE
         DO 240 I = 1, NCV
            CALL DCOPY(IRANKY,CVY(1,I),1,WK,1)
            CALL DGEMV('T',IRANKY,NY,RDF,QX(NX+1,1),N,WK,1,0.0D0,
     *                 CVY(1,I),1)
  240    CONTINUE
      END IF
      END
