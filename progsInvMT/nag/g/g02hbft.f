      SUBROUTINE G02HBF(UCV,N,M,X,IX,A,Z,BL,BD,TOL,MAXIT,NITMON,NIT,WK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ITERATIVE ALGORITHM FOR THE COMPUTATION OF THE MATRIX A
C     (STANDARDIZED CASE, A LOWER TRIANGULAR)
C     WHERE INV(A'A) IS A ROBUST EXTIMATE OF X'X
C     BASED ON ROUTINES IN ROBETH BY A. MARAZZI
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BD, BL, TOL
      INTEGER           IFAIL, IX, M, MAXIT, N, NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  A(M*(M+1)/2), WK(M*(M+1)/2), X(IX,M), Z(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  UCV
      EXTERNAL          UCV
C     .. Local Scalars ..
      DOUBLE PRECISION  CX, SMAX, U, XN, ZNR
      INTEGER           I, I1, IERROR, IFLAG, IJ, IL, J, J1, J2, JJ,
     *                  JJ1, JJ2, JL, KK, L, MM, MM4, NOUT, NREC
      LOGICAL           AZERO
      CHARACTER*80      REC
C     .. Local Arrays ..
      CHARACTER*80      EREC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           P01ABF
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTPMV, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, MOD, DBLE
C     .. Executable Statements ..
C
C     PARAMETER CHECK AND INITIALIZATION
C
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (EREC(1),FMT=99999) M
      ELSE IF (N.LT.2) THEN
         WRITE (EREC(1),FMT=99998) N
      ELSE IF (N.LT.M) THEN
         WRITE (EREC(1),FMT=99997) N, M
      ELSE IF (IX.LT.N) THEN
         WRITE (EREC(1),FMT=99996) IX, N
      ELSE
         IERROR = 2
      END IF
      IF (IERROR.NE.1) THEN
         AZERO = .FALSE.
         IJ = 0
         DO 20 I = 1, M
            IJ = IJ + I
            IF (A(IJ).EQ.0.0D0) GO TO 40
   20    CONTINUE
         GO TO 60
C
   40    AZERO = .TRUE.
   60    IF (AZERO) THEN
            WRITE (EREC(1),FMT=99995) I
         ELSE IF (TOL.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99994) TOL
         ELSE IF (MAXIT.LE.0) THEN
            WRITE (EREC(1),FMT=99993) MAXIT
         ELSE IF (BL.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99992) BL
         ELSE IF (BD.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99991) BD
         ELSE
            IERROR = 0
         END IF
         IF (IERROR.NE.2) THEN
C
C           SET UP CONSTANTS ETC
C
            NIT = 0
            XN = DBLE(N)
            MM = (M+1)*M/2
C
C           SET UP FOR ITERATION MONITORING
C
            IF (NITMON.GT.0) THEN
               IFLAG = 0
               CALL X04ABF(IFLAG,NOUT)
               WRITE (REC,FMT=99988)
               CALL X04BAF(NOUT,REC)
            END IF
   80       CONTINUE
C
C           ITERATIONS
C
            NIT = NIT + 1
C
C           COMPUTE AVERAGES
C
            DO 100 I = 1, MM
               WK(I) = 0.0D0
  100       CONTINUE
            DO 180 L = 1, N
               DO 120 J = 1, M
                  Z(J) = X(L,J)
  120          CONTINUE
               CALL DTPMV('U','T','N',M,A,Z,1)
               ZNR = DNRM2(M,Z,1)
               U = UCV(ZNR)
               IF (U.LT.0.0D0) THEN
                  GO TO 360
C
               ELSE
                  IJ = 0
                  DO 160 I = 1, M
                     DO 140 J = 1, I
                        IJ = IJ + 1
                        WK(IJ) = (Z(I)*U)*Z(J) + WK(IJ)
  140                CONTINUE
  160             CONTINUE
               END IF
  180       CONTINUE
C
C           FIND IMPROVEMENT MATRIX (WK) FOR A;
C           TRUNCATE IF NECESSARY; FIND MAXIMUM IMPROVEMENT
C
            SMAX = 0.0D0
            IJ = 0
            DO 220 I = 1, M
               IF (I.NE.1) THEN
                  I1 = I - 1
                  DO 200 J = 1, I1
                     IJ = IJ + 1
                     WK(IJ) = -MIN(MAX(WK(IJ)/XN,-BL),BL)
                     SMAX = MAX(SMAX,ABS(WK(IJ)))
  200             CONTINUE
               END IF
               IJ = IJ + 1
               CX = -MIN(MAX((WK(IJ)/XN-1.0D0)/2.0D0,-BD),BD)
               SMAX = MAX(SMAX,ABS(CX))
               WK(IJ) = CX + 1.0D0
  220       CONTINUE
C
C           FIND NEW TRANSFORMATION MATRIX A=WK*A
C
            IL = 1
            DO 240 I = 1, M
               CALL DTPMV('U','N','N',I,A,WK(IL),1)
               IL = IL + I
  240       CONTINUE
            CALL DCOPY(MM,WK,1,A,1)
C
C           ITERATION MONITORING
C
            IF (NITMON.GT.0) THEN
               IF ((MOD(NIT,NITMON).EQ.0) .OR. (NIT.EQ.1)) THEN
                  WRITE (REC,FMT=99987) NIT, SMAX
                  CALL X04BAF(NOUT,REC)
                  WRITE (REC,FMT=99986)
                  CALL X04BAF(NOUT,REC)
                  WRITE (REC,FMT=99985)
                  CALL X04BAF(NOUT,REC)
                  MM4 = MIN(7,M)
                  J2 = 0
                  DO 260 J = 1, MM4
                     J1 = J2 + 1
                     J2 = J2 + J
                     WRITE (REC,FMT=99984) J, (A(KK),KK=J1,J2)
                     CALL X04BAF(NOUT,REC)
  260             CONTINUE
                  DO 300 J = 8, M
                     J1 = J2 + 1
                     J2 = J2 + J
                     JL = (J2-J1)/7 + 1
                     DO 280 JJ = 1, JL
                        JJ1 = (JJ-1)*7 + J1
                        JJ2 = MIN(J2,JJ1+6)
                        IF (JJ.EQ.1) THEN
                           WRITE (REC,FMT=99984) J, (A(KK),KK=JJ1,JJ2)
                        ELSE
                           WRITE (REC,FMT=99983) (A(KK),KK=JJ1,JJ2)
                        END IF
                        CALL X04BAF(NOUT,REC)
  280                CONTINUE
  300             CONTINUE
               END IF
            END IF
C
C           STOP ITERATIONS IF GIVEN TOLERANCE OR MAXIMUM
C           NUMBER OF ITERATIONS IS REACHED
C
            IF (SMAX.GE.TOL .AND. NIT.LT.MAXIT) GO TO 80
            IF (SMAX.GE.TOL .AND. NIT.GE.MAXIT) THEN
               IERROR = 4
               NREC = 2
               WRITE (EREC,FMT=99989) MAXIT
            ELSE
C
C              CALCULATE Z
C
               DO 340 I = 1, N
                  DO 320 J = 1, M
                     WK(J) = X(I,J)
  320             CONTINUE
                  CALL DTPMV('U','T','N',M,A,WK,1)
                  Z(I) = DNRM2(M,WK,1)
  340          CONTINUE
            END IF
            GO TO 380
C
  360       IERROR = 3
            WRITE (EREC(1),FMT=99990) ZNR, U
         END IF
      END IF
  380 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,EREC)
      RETURN
C
99999 FORMAT (' ** On entry, M .lt. 1: M=',I16)
99998 FORMAT (' ** On entry, N .lt. 2: N=',I16)
99997 FORMAT (' ** On entry, N .lt. M: N=',I16,' M=',I16)
99996 FORMAT (' ** On entry, IX .lt. N: N=',I16,' IX=',I16)
99995 FORMAT (' ** On entry, diagonal element ',I4,' of A is 0')
99994 FORMAT (' ** On entry, TOL .le. 0: TOL=',D13.5)
99993 FORMAT (' ** On entry, MAXIT .le. 0: MAXIT=',I16)
99992 FORMAT (' ** On entry, BL .le. 0: BL=',D13.5)
99991 FORMAT (' ** On entry, BD .le. 0: BD=',D13.5)
99990 FORMAT (' ** Value returned by UCV function .lt. 0: U(',D13.5,
     *       ') =',D13.5)
99989 FORMAT (' ** Iterations to calculate weights failed to converge',
     *       /'    in MAXIT iterations: MAXIT = ',I16)
99988 FORMAT ('                    ** ITERATION MONITORING **')
99987 FORMAT (//' Iteration',I5,' MAX(ABS(S(I,J)))=',D13.5)
99986 FORMAT ('       A')
99985 FORMAT (' Row')
99984 FORMAT (1X,I2,7D10.2)
99983 FORMAT (3X,7D10.2)
      END
