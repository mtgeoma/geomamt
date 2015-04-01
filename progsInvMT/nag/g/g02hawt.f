      SUBROUTINE G02HAW(X,A,CUCV,N,M,MM,MDX,MAXIT,NITMON,TOL,NIT,SZ,SC2,
     *                  IFUN,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-691 (DEC 1989).
C
C     ITERATIVE ALGORITHM FOR THE COMPUTATION OF THE MATRIX A
C     (STANDARDIZED CASE, A LOWER TRIANGULAR)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CUCV, TOL
      INTEGER           IFAIL, IFUN, M, MAXIT, MDX, MM, N, NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MM), SC2(MM), SZ(N), X(MDX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  BLII, BLIJ, BUII, BUIJ, CX, SMAX, SV, U, XN, ZNR
      INTEGER           I, I1, IFLAG, IJ, IL, J, J1, J2, JJ, JJ1, JJ2,
     *                  JL, KK, L, MM4, NN, NOUT
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  G02HAV, DNRM2
      EXTERNAL          G02HAV, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTPMV, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, MOD, DBLE
C     .. Data statements ..
      DATA              BLIJ, BUIJ, BLII, BUII/-0.9D0, 0.9D0, -0.9D0,
     *                  0.9D0/
C     .. Executable Statements ..
C
C     PARAMETER CHECK AND INITIALIZATION
C
      IFAIL = 0
      NN = (M+1)*M/2
      NIT = 0
      XN = DBLE(N)
C
C     SET UP FOR ITERATION MONITORING
C
      IF (NITMON.GT.0) THEN
         IFLAG = 0
         CALL X04ABF(IFLAG,NOUT)
         WRITE (REC,FMT=99999)
         CALL X04BAF(NOUT,REC)
      END IF
C
   20 CONTINUE
C
C     ITERATIONS
C
      NIT = NIT + 1
C
C     COMPUTE AVERAGES
C
      DO 40 I = 1, MM
         SC2(I) = 0.0D0
   40 CONTINUE
      SV = 0.0D0
      DO 120 L = 1, N
         DO 60 J = 1, M
            SZ(J) = X(L,J)
   60    CONTINUE
         CALL DTPMV('U','T','N',M,A,SZ,1)
         ZNR = DNRM2(M,SZ,1)
         U = G02HAV(ZNR,IFUN,CUCV)
         IJ = 0
         DO 100 I = 1, M
            DO 80 J = 1, I
               IJ = IJ + 1
               SC2(IJ) = (SZ(I)*U)*SZ(J) + SC2(IJ)
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
C
C     FIND IMPROVEMENT MATRIX (SC2) FOR A;
C     TRUNCATE IF NECESSARY; FIND MAXIMUM IMPROVEMENT
C
      SMAX = 0.0D0
      IJ = 0
      DO 160 I = 1, M
         IF (I.NE.1) THEN
            I1 = I - 1
            DO 140 J = 1, I1
               IJ = IJ + 1
               SC2(IJ) = -MIN(MAX(SC2(IJ)/XN,BLIJ),BUIJ)
               SMAX = MAX(SMAX,ABS(SC2(IJ)))
  140       CONTINUE
         END IF
C
         IJ = IJ + 1
         CX = -MIN(MAX((SC2(IJ)/XN-1.0D0)/2.0D0,BLII),BUII)
         SMAX = MAX(SMAX,ABS(CX))
         SC2(IJ) = CX + 1.0D0
  160 CONTINUE
C
C     FIND NEW TRANSFORMATION MATRIX A=SC2*A
C
      IL = 1
      DO 180 I = 1, M
         CALL DTPMV('U','N','N',I,A,SC2(IL),1)
         IL = IL + I
  180 CONTINUE
      CALL DCOPY(MM,SC2,1,A,1)
C
C     ITERATION MONITORING
C
      IF (NITMON.GT.0) THEN
         IF ((MOD(NIT,NITMON).EQ.0) .OR. (NIT.EQ.1)) THEN
            WRITE (REC,FMT=99998) NIT, SMAX
            CALL X04BAF(NOUT,REC)
            WRITE (REC,FMT=99997)
            CALL X04BAF(NOUT,REC)
            WRITE (REC,FMT=99996)
            CALL X04BAF(NOUT,REC)
            MM4 = MIN(7,M)
            J2 = 0
            DO 200 J = 1, MM4
               J1 = J2 + 1
               J2 = J2 + J
               WRITE (REC,FMT=99995) J, (A(KK),KK=J1,J2)
               CALL X04BAF(NOUT,REC)
  200       CONTINUE
            DO 240 J = 8, M
               J1 = J2 + 1
               J2 = J2 + J
               JL = (J2-J1)/7 + 1
               DO 220 JJ = 1, JL
                  JJ1 = (JJ-1)*7 + J1
                  JJ2 = MIN(J2,JJ1+6)
                  IF (JJ.EQ.1) THEN
                     WRITE (REC,FMT=99995) J, (A(KK),KK=JJ1,JJ2)
                  ELSE
                     WRITE (REC,FMT=99994) (A(KK),KK=JJ1,JJ2)
                  END IF
                  CALL X04BAF(NOUT,REC)
  220          CONTINUE
  240       CONTINUE
         END IF
C
      END IF
C
C     STOP ITERATIONS IF GIVEN TOLERANCE OR MAXIMUM
C     NUMBER OF ITERATIONS IS REACHED
C
      IF (SMAX.GE.TOL .AND. NIT.LT.MAXIT) GO TO 20
      IF (SMAX.GE.TOL .AND. NIT.GE.MAXIT) THEN
         IFAIL = 1
C
      ELSE
C
C        CALCULATE Z
C
         DO 280 I = 1, N
            DO 260 J = 1, M
               SC2(J) = X(I,J)
  260       CONTINUE
            CALL DTPMV('U','T','N',M,A,SC2,1)
            SZ(I) = DNRM2(M,SC2,1)
  280    CONTINUE
      END IF
      RETURN
C
99999 FORMAT ('                ** ITERATION MONITORING FOR WEIGHTS **')
99998 FORMAT (' Iteration',I5,' MAX(ABS(S(I,J)))=',D13.5)
99997 FORMAT ('       A')
99996 FORMAT (' Row')
99995 FORMAT (1X,I2,7D10.2)
99994 FORMAT (3X,7D10.2)
      END
