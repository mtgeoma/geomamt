      SUBROUTINE G02DEF(WEIGHT,N,IP,Q,LDQ,P,WT,X,RSS,TOL,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G02DEF adds a new variable to a regression. The R and Q'Y matrices
C     are updated.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS, TOL
      INTEGER           IFAIL, IP, LDQ, N
      CHARACTER*1       WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  P(IP+1), Q(LDQ,IP+2), WT(*), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQTW, ZETA
      INTEGER           I, IERROR, IFAULT, IP2, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP(1)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QDF, F06FRF, F06FTF, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (N.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (IP.LT.0) THEN
         WRITE (P01REC(1),FMT=99998) IP
      ELSE IF (IP.GE.N) THEN
         WRITE (P01REC(1),FMT=99993) IP, N
      ELSE IF (LDQ.LT.N) THEN
         WRITE (P01REC(1),FMT=99996) LDQ, N
      ELSE IF (TOL.LE.0.0D0) THEN
         WRITE (P01REC(1),FMT=99992) TOL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         IP2 = IP + 2
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            IF (IP.EQ.0) THEN
C
C       CHECK WEIGHTS IF INITIAL CALL TO ROUTINE
C
               DO 20 I = 1, N
                  IF (WT(I).LT.0.0D0) THEN
                     GO TO 40
C
                  ELSE
                     SQTW = SQRT(WT(I))
                     Q(I,2) = X(I)*SQTW
                     Q(I,1) = Q(I,1)*SQTW
                  END IF
   20          CONTINUE
               GO TO 100
C
   40          IERROR = 2
               WRITE (P01REC(1),FMT=99994) I
               GO TO 120
C
            ELSE
               DO 60 I = 1, N
                  IF (WT(I).LT.0.0D0) THEN
                     GO TO 80
C
                  ELSE
                     Q(I,IP2) = SQRT(WT(I))*X(I)
                  END IF
   60          CONTINUE
               GO TO 100
C
   80          IERROR = 2
               WRITE (P01REC(1),FMT=99994) I
               GO TO 120
C
            END IF
         ELSE IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
            CALL DCOPY(N,X,1,Q(1,IP2),1)
         ELSE
            IERROR = 1
            WRITE (P01REC(1),FMT=99997) WEIGHT
            GO TO 120
C
         END IF
  100    IF (IP.NE.0) THEN
            IFAULT = 1
            CALL F01QDF('T','S',N,IP,Q(1,2),LDQ,P,1,Q(1,IP2),N,WKSP,
     *                  IFAULT)
         END IF
         CALL F06FRF(N-IP-1,Q(IP+1,IP2),Q(IP2,IP2),1,TOL,ZETA)
         IF (ZETA.LT.1.0D0) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99995)
            P(IP+1) = 0.0D0
         ELSE
            P(IP+1) = ZETA
            CALL F06FTF(N-IP-1,Q(IP+1,1),Q(IP2,1),1,ZETA,Q(IP2,IP2),1)
            RSS = DDOT(N-IP-1,Q(IP2,1),1,Q(IP2,1),1)
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, IP.lt.0 : IP = ',I16)
99997 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99996 FORMAT (' ** On entry, LDQ.lt.N : LDQ = ',I16,' N = ',I16)
99995 FORMAT (' ** X variable is a linear combination of existing mode',
     *       'l terms')
99994 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99993 FORMAT (' ** On entry, IP.ge.N : IP = ',I16,' N = ',I16)
99992 FORMAT (' ** On entry, TOL.le.0.0 : TOL = ',D13.5)
      END
