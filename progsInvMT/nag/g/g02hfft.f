      SUBROUTINE G02HFF(PSI,PSP,INDW,INDC,SIGMA,N,M,X,IX,RS,WGT,C,IC,WK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1034 (JUN 1993).
C     MARK 17 REVISED. IER-1678 (JUN 1995).
C
C     CALCULATION OF ROBUST ASYMPTOTIC VARIANCE-COVARIANCE MATRIX
C     FOR ROBUST REGRESSION
C     FOR HUBER TYPE REGRESSION C=FACT*INV(X'X)
C     FOR MALLOWS OR SCHWEPPE TYPE REGRESSION
C           C=FACT*INV(S1)*S2*INV(S1)
C     WHERE S1=(1/N)*(X'DX)  AND S2=(1/N)*(X'EX)
C
C     BASED ON ROUTINES FROM ROBETH BY A. MARAZZI
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGMA
      INTEGER           IC, IFAIL, INDC, INDW, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(IC,M), RS(N), WGT(N), WK(M*(N+M+1)+2*N),
     *                  X(IX,M)
C     .. Function Arguments ..
      DOUBLE PRECISION  PSI, PSP
      EXTERNAL          PSI, PSP
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, FACT, PS, RN, S, S1, S2, SP, SUM2, TMP1,
     *                  TMP2, VAR, XKAPPA, XMU, XMU2, Z
      INTEGER           I, IERROR, IFAIL2, IJ, J, J1, M1, MM1, NREC,
     *                  NUMC
      LOGICAL           WTNEG
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AAZ, F01ACZ, F06FCF, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     CHECK FOR ERRORS IN INPUT PARAMETERS
C
      NREC = 1
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (REC(1),FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (REC(1),FMT=99998) M
      ELSE IF (N.LE.M) THEN
         WRITE (REC(1),FMT=99997) N, M
      ELSE IF (IX.LT.N) THEN
         WRITE (REC(1),FMT=99996) N, IX
      ELSE IF (IC.LT.M) THEN
         WRITE (REC(1),FMT=99995) M, IC
      ELSE
         IERROR = 2
      END IF
      IF (IERROR.NE.1) THEN
         WTNEG = .FALSE.
         IF (INDW.NE.0) THEN
            DO 20 I = 1, N
               IF (WGT(I).LT.0.0D0) GO TO 40
   20       CONTINUE
            GO TO 60
C
   40       WTNEG = .TRUE.
         END IF
   60    IF (WTNEG) THEN
            WRITE (REC(1),FMT=99990)
         ELSE IF (SIGMA.LE.0.0D0) THEN
            WRITE (REC(1),FMT=99994) SIGMA
         ELSE
            IERROR = 0
         END IF
         IF (IERROR.NE.2) THEN
            M1 = M + 1
            MM1 = M1*M + 1
C
C           FOR HUBER TYPE REGRESSION
C
            IF (INDW.EQ.0) THEN
C
C              CALCULATE INV(X'X)
C
               EPS = X02AJF()
               DO 100 J = 1, M
                  DO 80 I = 1, J
                     IJ = (J-1)*(M+1) + I
                     WK(IJ) = DDOT(N,X(1,I),1,X(1,J),1)
   80             CONTINUE
  100          CONTINUE
               IFAIL2 = 1
               CALL F01ACZ(M,EPS,WK,M1,WK(MM1),N,WK(2*MM1),NUMC,IFAIL2)
               IF (IFAIL2.NE.0) THEN
                  IERROR = 3
                  WRITE (REC(1),FMT=99993)
               ELSE
C                 -------
C                 COMPUTES CORRECTION FACTORS XKAPPA AND SUM2 FOR
C                 THE COVARIANCE MATRIX.
C
                  TMP1 = 0.0D0
                  TMP2 = 0.0D0
                  DO 120 J = 1, N
                     S = RS(J)/SIGMA
                     SP = PSP(S)
                     WK(J+MM1) = SP
                     TMP1 = SP + TMP1
                     PS = PSI(S)
                     TMP2 = PS*PS + TMP2
  120             CONTINUE
                  XMU = TMP1/DBLE(N)
                  SUM2 = TMP2
                  VAR = 0.0D0
                  DO 140 J = 1, N
                     VAR = (WK(J+MM1)-XMU)**2 + VAR
  140             CONTINUE
                  VAR = VAR/DBLE(N)
                  XKAPPA = 0.0D0
                  XMU2 = XMU*XMU
                  IF (XMU2.GT.0.0D0) XKAPPA = 1.0D0 + DBLE(M)/DBLE(N)
     *                *VAR/XMU2
                  IF (XMU2.GT.0.0D0) SUM2 = SUM2/XMU2/DBLE(N-M)
                  IF (XKAPPA.LE.0.0D0 .OR. SUM2.LE.0.0D0) THEN
                     IERROR = 4
                     WRITE (REC(1),FMT=99991)
                     FACT = 1.0D0
                  ELSE
                     FACT = (XKAPPA*XKAPPA)*SUM2
                     FACT = FACT*SIGMA*SIGMA
                  END IF
                  DO 180 J = 1, M
                     J1 = J + 1
                     IJ = (J-1)*(M+1) + J1
                     C(J,J) = WK(IJ)*FACT
                     DO 160 I = J1, M
                        IJ = (J-1)*(M+1) + I + 1
                        C(I,J) = WK(IJ)*FACT
                        C(J,I) = C(I,J)
  160                CONTINUE
  180             CONTINUE
               END IF
            ELSE
C
C              FOR MALLOWS AND SCHWEPPE CASES
C
C              CALCULATE VALUES OF D AND E
C              WK(I), I=1,N CONTAINS D
C              WK(N+I), I=1,N CONTAINS E
C
               IF (INDC.EQ.1) THEN
C
C                 CASE INDC = 1 : AVERAGING OVER THE WEIGHTS
C
                  IF (INDW.LT.0) THEN
C
C                    MALLOWS CASE
C
                     S1 = 0.0D0
                     S2 = 0.0D0
                     DO 200 J = 1, N
                        IF (WGT(J).GT.0.0D0) THEN
                           Z = RS(J)/SIGMA
                           S1 = PSP(Z) + S1
                           S2 = PSP(Z)**2 + S2
                        END IF
  200                CONTINUE
                     DO 220 I = 1, N
                        WK(I) = S1/DBLE(N)*WGT(I)
                        WK(N+I) = S2/DBLE(N)*WGT(I)*WGT(I)
  220                CONTINUE
                  ELSE
C
C                    SCHWEPPE CASE
C
                     DO 260 I = 1, N
                        S1 = 0.0D0
                        S2 = 0.0D0
                        IF (WGT(I).GT.0.0D0) THEN
                           DO 240 J = 1, N
                              Z = RS(J)/SIGMA/WGT(I)
                              S1 = PSP(Z) + S1
                              S2 = PSI(Z)**2 + S2
  240                      CONTINUE
                        END IF
                        WK(I) = S1/DBLE(N)
                        WK(N+I) = S2/DBLE(N)*WGT(I)*WGT(I)
  260                CONTINUE
                  END IF
C
C                 DEFAULT CASE EXPECTED EQUAL OBSERVED
C
               ELSE IF (INDW.LT.0) THEN
C
C                 MALLOWS CASE
C
                  DO 280 I = 1, N
                     Z = RS(I)/SIGMA
                     WK(I) = PSP(Z)*WGT(I)
                     WK(N+I) = (PSI(Z)*WGT(I))**2
  280             CONTINUE
               ELSE
C
C                 SCHWEPPE CASE
C
                  DO 300 I = 1, N
                     IF (WGT(I).LE.0.0D0) THEN
                        WK(I) = 0.0D0
                        WK(N+I) = 0.0D0
                     ELSE
                        Z = RS(I)/SIGMA/WGT(I)
                        WK(I) = PSP(Z)
                        WK(N+I) = (PSI(Z)*WGT(I))**2
                     END IF
  300             CONTINUE
               END IF
C
C              CALCULATE ESTIMATE OF THE FORM INV(S1)*S2*INV(S1)
C
               RN = DBLE(N)
               FACT = (SIGMA*SIGMA)/RN
C
C              CALCULATE S1
C
               DO 320 J = 1, M
                  IJ = (J-1)*N + 1
                  CALL DCOPY(N,X(1,J),1,WK(2*N+IJ),1)
                  CALL F06FCF(N,WK(1),1,WK(2*N+IJ),1)
  320          CONTINUE
               DO 360 J = 1, M
                  J1 = J - 1
                  IJ = (J-1)*N + 1
                  C(J,J) = DDOT(N,X(1,J),1,WK(2*N+IJ),1)/RN
                  DO 340 I = 1, J1
                     C(I,J) = DDOT(N,X(1,I),1,WK(2*N+IJ),1)/RN
                     C(J,I) = C(I,J)
  340             CONTINUE
  360          CONTINUE
C
C              CALCULATE INV(S1)
C
               IFAIL2 = 1
               CALL F01AAZ(C,IC,M,WK(2*N+1),M,WK(2*N+MM1),IFAIL2)
               IF (IFAIL2.NE.0) THEN
                  IERROR = 3
                  WRITE (REC(1),FMT=99992)
               ELSE
C
C                 CALCULATE S2
C
                  DO 380 J = 1, M
                     IJ = (J-1)*N + MM1
                     CALL DCOPY(N,X(1,J),1,WK(2*N+IJ),1)
                     CALL F06FCF(N,WK(N+1),1,WK(2*N+IJ),1)
  380             CONTINUE
                  DO 420 J = 1, M
                     DO 400 I = 1, J
                        IJ = (J-1)*N + MM1
                        C(I,J) = DDOT(N,X(1,I),1,WK(2*N+IJ),1)/RN
  400                CONTINUE
  420             CONTINUE
C
C                 CALCULATE INV(S1)*S2
C
                  DO 460 J = 1, M
                     DO 440 I = 1, M
                        IJ = (J-1)*M + I + MM1
                        WK(2*N+IJ) = DDOT(J,WK(2*N+I),M,C(1,J),1)
                        IF (J.NE.M) WK(2*N+IJ) = WK(2*N+IJ) + DDOT(M-J,
     *                      WK(2*N+I+J*M),M,C(J,J+1),IC)
  440                CONTINUE
  460             CONTINUE
C
C                 CALCULATE INV(S1)*S2*INV(S1)
C
                  DO 500 J = 1, M
                     J1 = J - 1
                     IJ = (J-1)*M + 1
                     C(J,J) = DDOT(M,WK(2*N+J+MM1),M,WK(2*N+IJ),1)*FACT
                     DO 480 I = 1, J1
                        C(I,J) = DDOT(M,WK(2*N+I+MM1),M,WK(2*N+IJ),1)
     *                           *FACT
                        C(J,I) = C(I,J)
  480                CONTINUE
  500             CONTINUE
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, N lt 2: N=',I16)
99998 FORMAT (' ** On entry, M lt 1: M=',I16)
99997 FORMAT (' ** On entry, N le M: N=',I16,' M=',I16)
99996 FORMAT (' ** On entry, IX lt N: N=',I16,' IX=',I16)
99995 FORMAT (' ** On entry, IC lt M: M=',I16,' IC=',I16)
99994 FORMAT (' ** On entry, SIGMA lt 0: SIGMA=',D13.5)
99993 FORMAT (' ** X''X Matrix not positive definite')
99992 FORMAT (' ** S1 Matrix is singular or almost singular')
99991 FORMAT (' ** Correction factor = 0 (Huber type regression)')
99990 FORMAT (' ** On entry, an element of WGT lt 0')
      END
