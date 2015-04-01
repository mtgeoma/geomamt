      SUBROUTINE G02HKF(N,M,X,LDX,EPS,COV,THETA,MAXIT,NITMON,TOL,NIT,WK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-921 (APR 1991).
C
C     MASTER ROUTINE FOR ROBUST COVARIANCE
C     BASED ON ROUTINES IN ROBETH BY A. MARAZZI
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HKF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, TOL
      INTEGER           IFAIL, LDX, M, MAXIT, N, NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  COV(M*(M+1)/2), THETA(M), WK(N+M*(M+5)/2),
     *                  X(LDX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  BD, BL, SW, T, XMD, XSD
      INTEGER           I, IERROR, IFAULT, IJ, J, MM, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  B(3)
      CHARACTER*80      EREC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02BUF, G02HKU, G02HKW, G02HKX, G02HKZ, G02HLF,
     *                  G07DAF, DSCAL
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (EREC(1),FMT=99999) M
      ELSE IF (N.LT.2) THEN
         WRITE (EREC(1),FMT=99998) N
      ELSE IF (N.LT.M) THEN
         WRITE (EREC(1),FMT=99997) N, M
      ELSE IF (LDX.LT.N) THEN
         WRITE (EREC(1),FMT=99996) LDX, N
      ELSE IF (EPS.LT.0.0D0 .OR. EPS.GE.1.0D0) THEN
         WRITE (EREC(1),FMT=99995) EPS
      ELSE IF (TOL.LE.0.0D0) THEN
         WRITE (EREC(1),FMT=99994) TOL
      ELSE IF (MAXIT.LE.0) THEN
         WRITE (EREC(1),FMT=99993) MAXIT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         MM = (M*(M+1))/2
         IF (EPS.EQ.0.0D0) THEN
            IFAULT = 1
            CALL G02BUF('M','U',N,M,X,LDX,WK,SW,THETA,COV,IFAULT)
            T = 1.0D0/N
            CALL DSCAL(MM,T,COV,1)
         ELSE
C
C           FIND INITIAL ESTIMATE OF A AND CHECK FOR CONSTANT X
C
            IJ = 0
            DO 20 I = 1, M
               DO 10 J = IJ + 1, IJ + I - 1
                  WK(J) = 0.0D0
   10          CONTINUE
               IJ = IJ + I
               IFAULT = 1
               CALL G07DAF(N,X(1,I),WK(MM+1),THETA(I),XMD,XSD,IFAULT)
               IF (XMD.LE.0.0D0) THEN
                  XMD = WK(MM+N) - WK(MM+1)
                  IF (XMD.EQ.0.0D0) GO TO 40
                  XSD = XMD/2.0D0
               END IF
               WK(IJ) = 1.0D0/XSD
   20       CONTINUE
            GO TO 60
   40       IERROR = 2
            WRITE (EREC(1),FMT=99990) I
            GO TO 80
   60       CONTINUE
C
C           CALCULATE CONSTANTS
C
            CALL G02HKZ(EPS,M,B(1),B(2))
            CALL G02HKX(EPS,B(3))
            BL = 0.9D0
            BD = 0.9D0
C
C           CALL ROBUST COVARIANCE ROUTINE
C
            IFAULT = 1
            CALL G02HLF(G02HKU,B,0,N,M,X,LDX,COV,WK,WK(MM+1),THETA,BL,
     *                  BD,MAXIT,NITMON,TOL,NIT,WK(N+MM+1),IFAULT)
            IF (IFAULT.EQ.6) THEN
               IERROR = 4
               WRITE (EREC(1),FMT=99991)
               GO TO 80
            ELSE IF (IFAULT.EQ.5) THEN
               IERROR = 3
               WRITE (EREC(1),FMT=99992)
            END IF
C
C           SCALE RESULT
C
            CALL G02HKW(B(1),B(2),M,T)
            CALL DSCAL(MM,T,COV,1)
         END IF
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,EREC)
C
99999 FORMAT (' ** On entry, M.lt.1 : M =',I16)
99998 FORMAT (' ** On entry, N.lt.2 : N =',I16)
99997 FORMAT (' ** On entry, N.lt.M : N =',I16,' M =',I16)
99996 FORMAT (' ** On entry, LDX.lt.N : LDX =',I16,' N =',I16)
99995 FORMAT (' ** On entry, EPS.lt.0.0 or EPS.ge.1.0 : EPS =',D13.5)
99994 FORMAT (' ** On entry, TOL.le.0.0 : TOL =',D13.5)
99993 FORMAT (' ** On entry, MAXIT.le.0 : MAXIT =',I16)
99992 FORMAT (' ** Iterations have failed to converge')
99991 FORMAT (' ** Iterations have become unstable')
99990 FORMAT (' ** ',I16,' TH column of X has constant value')
      END
