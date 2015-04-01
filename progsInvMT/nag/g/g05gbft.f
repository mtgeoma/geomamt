      SUBROUTINE G05GBF(N,D,C,LDC,EPS,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       This algorithm randomly selects a correlation matrix from the
C       class of all correlation matrices with specified eigenvalues
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, FOUR, ZP25, ZP5
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0,
     *                  ZP25=0.25D0,ZP5=0.5D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05GBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IFAIL, LDC, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,N), D(N), WK(2*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CA, CB, CC, CII, CIIMJJ, CIIPJJ, CIJ, CIJIJ,
     *                  CJJ, CS, CSCS, CSSN, DT, SN, SNSN, SUM, THETA,
     *                  TT
      INTEGER           I, I1, IERROR, IFAULT, IG, IGTEMP, II, IL,
     *                  ILTEMP, J, J1, JJ, K, L, NK, NM, NN
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05GBZ, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, ACOS, COS, DBLE, INT, MAX, MIN, SIN, SQRT
C     .. Executable Statements ..
C
      IERROR = 1
      IF (N.LE.0) THEN
         WRITE (REC,FMT=99999) N
      ELSE IF (N.GT.LDC) THEN
         WRITE (REC,FMT=99998) N, LDC
      ELSE IF (EPS.LT.DBLE(N)*X02AJF()) THEN
         WRITE (REC,FMT=99997) EPS
      ELSE
         IERROR = 0
         TT = ZERO
         DO 20 I = 1, N
            WK(I) = D(I)
   20    CONTINUE
         IFAULT = 0
         CALL M01CAF(WK,1,N,'A',IFAULT)
         IF (WK(1).LT.ZERO) THEN
            IERROR = 2
            WRITE (REC,FMT=99996)
            GO TO 220
         ELSE
            DO 40 I = 1, N
               TT = TT + WK(I)
   40       CONTINUE
         END IF
         IF (ABS(TT-DBLE(N)).GT.EPS) THEN
            IERROR = 2
            WRITE (REC,FMT=99995)
            GO TO 220
         ELSE
C
C           Generate covariance matrix
C
            CALL G05GBZ(N,D,C,LDC,WK)
C
C           Randomly select two diagonal elements of C, one greater
C           than one and the other less than one, then construct the
C           appropriate 2x2 rotational matrix which sets one diagonal
C           element equal to one while maintaining the structure of
C           eigenvalues
C
            NK = 0
            SUM = ZERO
            DO 60 I = 1, N
               WK(I) = -2.0D0
               IF (ABS(C(I,I)-ONE).LE.EPS) THEN
                  NK = NK + 1
                  WK(I) = 0.0D0
               END IF
   60       CONTINUE
            NM = N - NK - 1
            IF (NM.EQ.-1) GO TO 240
            IF (NM.NE.0) THEN
               DO 180 NN = 1, NM
                  IG = 0
                  IL = 0
                  DO 80 K = 1, N
C
C                    WK(1..N) stores indicators showing whether diagonal
C                    elements of C are .EQ. ONE, .GT. ONE, and .LT. ONE
C
                     IF (WK(K).NE.0.0D0) THEN
                        IF (C(K,K).GT.(ONE+EPS)) THEN
                           IG = IG + 1
                           WK(K) = 1.0D0
                        ELSE
                           IL = IL + 1
                           WK(K) = -1.0D0
                        END IF
                     END IF
   80             CONTINUE
                  IF (IL*IG.EQ.0) THEN
                     GO TO 200
                  ELSE
                     II = INT(G05CAF(0.0D0)*IG) + 1
                     JJ = INT(G05CAF(0.0D0)*IL) + 1
                     I1 = 0
                     J1 = 0
                     IGTEMP = 0
                     ILTEMP = 0
                     DO 100 K = 1, N
                        IF (I1.EQ.0) THEN
                           IF (WK(K).EQ.1.0D0) IGTEMP = IGTEMP + 1
                           IF (IGTEMP.EQ.II) I1 = K
                        END IF
                        IF (J1.EQ.0) THEN
                           IF (WK(K).EQ.-1.0D0) ILTEMP = ILTEMP + 1
                           IF (ILTEMP.EQ.JJ) J1 = K
                        END IF
  100                CONTINUE
C
                     I = MIN(I1,J1)
                     J = MAX(I1,J1)
                     IF (ABS(C(I,I)-ONE).GT.EPS) THEN
C
C                       Solve the quadratic equation and determine the
C                       rotation angle
C
                        CII = C(I,I)
                        CIJ = C(I,J)
                        CJJ = C(J,J)
                        CIIMJJ = CII - CJJ
                        CIIPJJ = CII + CJJ - TWO
                        CIJIJ = CIJ*CIJ
                        CA = CIJIJ + ZP25*CIIMJJ*CIIMJJ
                        CB = ZP5*CIIPJJ*CIIMJJ
                        CC = ZP25*CIIPJJ*CIIPJJ - CIJIJ
                        DT = CB*CB - FOUR*CA*CC
                        IF (DT.LT.0.0D0) THEN
                           DT = 0.0D0
                        ELSE
                           DT = SQRT(DT)
                        END IF
C
C                       Randomly select one of the two possible
C                       solutions
C
                        IF (G05CAF(0.0D0).GE.ZP5) DT = -DT
                        THETA = ACOS(ZP5*(-CB+DT)/CA)/TWO
                        CS = COS(THETA)
                        SN = SIN(THETA)
                        CSCS = CS*CS
                        CSSN = CS*SN
                        SNSN = SN*SN
C
C                       Choose the correct direction of rotation
C
                        IF (ABS(CSCS*CII+SNSN*CJJ+TWO*CSSN*CIJ-ONE)
     *                      .GT.EPS) THEN
                           SN = -SN
                           CSSN = -CSSN
                        END IF
C
C                       Compute C=PCP' using the fact that P is a
C                       rotational matrix with the rotation being
C                       operated by 2x2 matrix (CS SN)/(-SN CS). only
C                       elements of I and J-th rows (hence columns) of
C                       C are changed by the rotation.
C
                        C(I,I) = CSCS*CII + SNSN*CJJ + TWO*CSSN*CIJ
                        IF (ABS(C(I,I)-ONE).GT.EPS) THEN
                           GO TO 200
                        ELSE
                           WK(I) = 0.0D0
                           NK = NK + 1
                           C(J,J) = CSCS*CJJ + SNSN*CII - TWO*CSSN*CIJ
                           C(I,J) = (CSCS-SNSN)*CIJ + CSSN*(CJJ-CII)
                           C(J,I) = C(I,J)
                           DO 120 L = 1, N
                              IF (L.NE.I .AND. L.NE.J) THEN
                                 IF (L.LT.I) THEN
                                    C(I,L) = CS*C(L,I) + SN*C(L,J)
                                    C(J,L) = -SN*C(L,I) + CS*C(L,J)
                                 ELSE IF (L.GT.J) THEN
                                    C(L,I) = CS*C(I,L) + SN*C(J,L)
                                    C(L,J) = -SN*C(I,L) + CS*C(J,L)
                                 ELSE
                                    C(L,I) = CS*C(I,L) + SN*C(L,J)
                                    C(J,L) = -SN*C(I,L) + CS*C(L,J)
                                 END IF
                              END IF
  120                      CONTINUE
                           DO 160 II = 1, N
                              DO 140 JJ = II + 1, N
                                 C(II,JJ) = C(JJ,II)
  140                         CONTINUE
  160                      CONTINUE
                        END IF
                     END IF
                     SUM = SUM + C(I,I)
                  END IF
  180          CONTINUE
C
C              Set the remaining diagonal element to be (N-SUM) and
C              see if it is within the specified precision limit
C
               C(N,N) = DBLE(N) - SUM
               IF (ABS(C(N,N)-ONE).GT.EPS) GO TO 200
               GO TO 240
            END IF
  200       IERROR = 3
            WRITE (REC,FMT=99994)
         END IF
  220    CONTINUE
      END IF
  240 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
      RETURN
99999 FORMAT (' ** On entry, N.lt.1: N = ',I16)
99998 FORMAT (' ** On entry, N.gt.LDC: N = ',I16,' LDC = ',I16)
99997 FORMAT (' ** On entry, EPS.lt.N*machine precision: EPS = ',D13.5)
99996 FORMAT (' ** On entry, an eigenvalue is negative.')
99995 FORMAT (' ** On entry, the eigenvalues do not sum to N.')
99994 FORMAT (' ** Diagonals of returned matrix are not unity.')
      END
