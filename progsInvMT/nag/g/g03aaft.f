      SUBROUTINE G03AAF(MATRIX,STD,WEIGHT,N,M,X,LDX,ISX,S,WT,NVAR,E,LDE,
     *                  P,LDP,V,LDV,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1035 (JUN 1993).
C     MARK 17 REVISED. IER-1661 (JUN 1995).
C
C     PRINCIPAL COMPONENT ANALYSIS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03AAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDE, LDP, LDV, LDX, M, N, NVAR
      CHARACTER         MATRIX, STD, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  E(LDE,6), P(LDP,NVAR), S(M), V(LDV,NVAR),
     *                  WK(NVAR*NVAR+5*(NVAR-1)), WT(*), X(LDX,M)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  R, RN, SCALE, SRN, TEMP, WSUM
      INTEGER           I, IERROR, IFAULT, IND0, IVAR, J, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1,1), WKSP2(1,1)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, G01ECF
      INTEGER           P01ABF
      EXTERNAL          DDOT, G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02WEF, G03AAZ, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) M
      ELSE IF (NVAR.LT.1) THEN
         WRITE (P01REC(1),FMT=99983) NVAR
      ELSE IF (NVAR.GT.M) THEN
         WRITE (P01REC(1),FMT=99982) NVAR, M
      ELSE IF (N.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) N
      ELSE IF (NVAR.GE.N) THEN
         WRITE (P01REC(1),FMT=99985) N, NVAR
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDX, N
      ELSE IF (LDV.LT.N) THEN
         WRITE (P01REC(1),FMT=99996) LDV, N
      ELSE IF (LDP.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99995) LDP, NVAR
      ELSE IF (LDE.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99994) LDE, NVAR
      ELSE IF (MATRIX.NE.'S' .AND. MATRIX.NE.'C' .AND. MATRIX.NE.
     *         's' .AND. MATRIX.NE.'c' .AND. MATRIX.NE.'U' .AND.
     *         MATRIX.NE.'V' .AND. MATRIX.NE.'u' .AND. MATRIX.NE.'v')
     *         THEN
         WRITE (P01REC(1),FMT=99993) MATRIX
      ELSE IF (STD.NE.'S' .AND. STD.NE.'U' .AND. STD.NE.'s' .AND.
     *         STD.NE.'u' .AND. STD.NE.'Z' .AND. STD.NE.'z' .AND.
     *         STD.NE.'E' .AND. STD.NE.'e') THEN
         WRITE (P01REC(1),FMT=99991) STD
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99992) WEIGHT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
C
C        FIND  NO OF SELECTED VARIABLES
C
         IVAR = 0
         DO 20 I = 1, M
            IF (ISX(I).GT.0) IVAR = IVAR + 1
   20    CONTINUE
         IF (IVAR.NE.NVAR) THEN
            IERROR = 3
            NREC = 2
            WRITE (P01REC,FMT=99989) IVAR, NVAR
            GO TO 320
         END IF
C
C        CHECK WEIGHTS
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            WSUM = 0.0D0
            SRN = 0.0D0
            DO 40 I = 1, N
               IF (WT(I).LT.0.0D0) GO TO 60
               IF (WT(I).GT.0.0D0) THEN
                  WSUM = WSUM + WT(I)
                  V(I,IVAR) = SQRT(WT(I))
               ELSE
                  V(I,IVAR) = 0.0D0
               END IF
   40       CONTINUE
            IF (DBLE(IVAR).GT.WSUM-1.0D0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99988)
               GO TO 320
            END IF
         ELSE
            WSUM = DBLE(N)
         END IF
         GO TO 80
   60    IERROR = 2
         WRITE (P01REC(1),FMT=99990) I
         GO TO 320
   80    CONTINUE
         SRN = SQRT(WSUM-1.0D0)
C
C        CHECK INPUT SCALE FACTORS
C
         IF (MATRIX.EQ.'S' .OR. MATRIX.EQ.'s') THEN
            DO 100 J = 1, M
               IF (ISX(J).GT.0) THEN
                  IF (S(J).LE.0.0D0) GO TO 120
               END IF
  100       CONTINUE
            GO TO 140
  120       IERROR = 4
            WRITE (P01REC(1),FMT=99987) J
            GO TO 320
  140       CONTINUE
         END IF
C
C        CREATE STANDARDIZED DATA MATRIX
C
         CALL G03AAZ(MATRIX,WEIGHT,N,X,LDX,M,ISX,IVAR,WT,WSUM,V,LDV,S,E)
         IFAULT = 1
         CALL F02WEF(N,IVAR,V,LDV,0,WKSP1,1,.TRUE.,WKSP2,1,E,.TRUE.,P,
     *               LDP,WK,IFAULT)
         IF (IFAULT.GT.0) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99986)
            GO TO 320
         END IF
         IF (MATRIX.EQ.'V' .OR. MATRIX.EQ.'v') THEN
            SCALE = 1.0D0/SRN
            CALL DSCAL(IVAR,SCALE,E,1)
         END IF
         IF (STD.EQ.'U' .OR. STD.EQ.'u') THEN
            DO 160 I = 1, IVAR
               CALL DSCAL(N,E(I,1),V(1,I),1)
  160       CONTINUE
         ELSE IF (STD.EQ.'Z' .OR. STD.EQ.'z') THEN
            DO 170 I = 1, IVAR
               CALL DSCAL(N,SRN,V(1,I),1)
 170        CONTINUE
         ELSE IF (STD.EQ.'E' .OR. STD.EQ.'e') THEN
            DO 175 I = 1, IVAR
               CALL DSCAL(N,E(I,1)*SRN,V(1,I),1)
 175        CONTINUE
         END IF
         DO 200 I = 1, IVAR
            DO 180 J = I + 1, IVAR
               TEMP = P(I,J)
               P(I,J) = P(J,I)
               P(J,I) = TEMP
  180       CONTINUE
  200    CONTINUE
         IF (MATRIX.EQ.'C' .OR. MATRIX.EQ.'c') THEN
            SCALE = DBLE(IVAR)
         ELSE
            SCALE = DDOT(IVAR,E,1,E,1)
         END IF
         IF (SCALE.LE.0.0D0) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99984)
            GO TO 320
         END IF
         E(IVAR,4) = 0.0D0
         E(IVAR,5) = 0.0D0
         E(IVAR,6) = 0.0D0
         E(1,1) = E(1,1)*E(1,1)
         E(1,2) = E(1,1)/SCALE
         E(1,3) = E(1,2)
         IF (IVAR.EQ.1) GO TO 320
         IND0 = 1
         DO 220 I = 2, IVAR
            E(I,1) = E(I,1)*E(I,1)
            IF (E(I,1).GT.0.0D0) IND0 = I
            E(I,2) = E(I,1)/SCALE
            E(I,3) = E(I-1,3) + E(I,2)
  220    CONTINUE
         DO 260 J = 4, 6
            DO 240 I = IND0 + 1, IVAR
               E(I,J) = 0.0D0
  240       CONTINUE
  260    CONTINUE
         WK(IND0) = LOG(E(IND0,1))
         WK(IVAR+IND0) = E(IND0,1)
         DO 280 I = IND0 - 1, 1, -1
            WK(I) = WK(I+1) + LOG(E(I,1))
            WK(I+IVAR) = WK(I+IVAR+1) + E(I,1)
  280    CONTINUE
         RN = (WSUM-1.0D0) - DBLE(2*IVAR+5)/6.0D0
         DO 300 I = 1, IND0 - 1
            R = DBLE(IND0-I+1)
            E(I,4) = RN*(R*LOG(WK(I+IVAR)/R)-WK(I))
            E(I,5) = 0.5D0*(R-1.0D0)*(R+2.0D0)
            IF (E(I,4).LE.0.0D0 .OR. MATRIX.EQ.'C' .OR. MATRIX.EQ.'c')
     *          THEN
               E(I,6) = 0.0D0
            ELSE
               IFAULT = 1
               E(I,6) = G01ECF('UPPER',E(I,4),E(I,5),IFAULT)
            END IF
  300    CONTINUE
      END IF
  320 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99998 FORMAT (' ** On entry, N.lt.2 : N = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, LDV.lt.N : LDV = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, LDP.lt.NVAR : LDP = ',I16,' NVAR = ',I16)
99994 FORMAT (' ** On entry, LDE.lt.NVAR : LDE = ',I16,' NVAR = ',I16)
99993 FORMAT (' ** On entry, MATRIX is not valid : MATRIX = ',A1)
99992 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99991 FORMAT (' ** On entry, STD is not valid : STD = ',A1)
99990 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99989 FORMAT (' ** On entry, ',I16,' values in ISX.gt.0',/'    there s',
     *       'hould be NVAR = ',I16)
99988 FORMAT (' ** Number of selected variables .ge. effective number ',
     *       'of observations')
99987 FORMAT (' ** On entry, S(',I16,').le.0.0')
99986 FORMAT (' ** SVD has failed to converge')
99985 FORMAT (' ** On entry, N.le.NVAR : N = ',I16,' NVAR = ',I16)
99984 FORMAT (' ** All eigenvalues are zero')
99983 FORMAT (' ** On entry, NVAR.lt.1 : NVAR = ',I16)
99982 FORMAT (' ** On entry, NVAR.gt.M : NVAR = ',I16,' M = ',I16)
      END
