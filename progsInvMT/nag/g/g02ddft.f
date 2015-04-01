      SUBROUTINE G02DDF(N,IP,Q,LDQ,RSS,IDF,B,SE,COV,SVD,IRANK,P,TOL,WK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     COMPUTES REGRESSION PARAMETERS ETC FROM R MATRIX AND Q'Y VECTOR
C     THESE ARE STORED IN Q
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS, TOL
      INTEGER           IDF, IFAIL, IP, IRANK, LDQ, N
      LOGICAL           SVD
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV(IP*(IP+1)/2), P(IP*IP+2*IP),
     *                  Q(LDQ,IP+1), SE(IP), WK(IP*IP+(IP-1)*5)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, RMS
      INTEGER           I, IERROR, IFAULT, IJ, IP2, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ, DDOT
      INTEGER           F06KLF, P01ABF
      EXTERNAL          F02WDZ, DDOT, F06KLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02WUF, F06FCF, G02AAX, G02AAY, G02AAZ, DCOPY,
     *                  DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (N.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) IP
      ELSE IF (LDQ.LT.IP) THEN
         WRITE (P01REC(1),FMT=99997) LDQ, IP
      ELSE IF (RSS.LE.0.0D0 .AND. LDQ.LT.N) THEN
         NREC = 2
         WRITE (P01REC,FMT=99994) LDQ, N
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99993) TOL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         COND = F02WDZ(IP,Q(1,2),LDQ,WK)
         IF (COND*TOL.GT.1.0D0) THEN
            SVD = .TRUE.
            DO 20 I = 2, IP + 1
               CALL DCOPY(I-1,Q(1,I),1,P(I*IP+1),1)
   20       CONTINUE
            CALL DCOPY(IP,Q,1,SE,1)
            IP2 = IP*IP
            IFAULT = 1
            CALL F02WUF(IP,P(2*IP+1),IP,1,SE,IP,.TRUE.,WK,IP,P(IP+1),
     *                  .TRUE.,WK(IP2+1),IFAULT)
            IF (IFAULT.LT.0) THEN
               IERROR = 3
               WRITE (P01REC(1),FMT=99996)
            ELSE
               IRANK = F06KLF(IP,P(IP+1),1,TOL)
               DO 40 I = 1, IRANK
                  WK(IP2+I) = 1.0D0/P(IP+I)
   40          CONTINUE
               DO 60 I = 2, IP + 1
                  CALL F06FCF(IRANK,WK(IP2+1),1,P(I*IP+1),1)
   60          CONTINUE
               CALL DGEMV('T',IRANK,IP,1.0D0,P(2*IP+1),IP,SE,1,0.0D0,B,
     *                    1)
               IJ = 1
               DO 80 I = 1, IP
                  CALL DCOPY(IRANK,P(IP*I+IP+1),1,WK(IP2+1),1)
                  CALL DGEMV('T',IRANK,I,1.0D0,P(2*IP+1),IP,WK(IP2+1),1,
     *                       0.0D0,COV(IJ),1)
                  IJ = IJ + I
   80          CONTINUE
               IF (RSS.LE.0.0D0) RSS = DDOT(IP-IRANK,SE(IRANK+1),1,
     *                                 SE(IRANK+1),1) + DDOT(N-IP,
     *                                 Q(IP+1,1),1,Q(IP+1,1),1)
            END IF
         ELSE
            SVD = .FALSE.
            IRANK = IP
            CALL DCOPY(IP,Q,1,B,1)
            CALL DTRSV('U','N','N',IP,Q(1,2),LDQ,B,1)
            IJ = 1
            DO 100 I = 1, IP
               CALL DCOPY(I,Q(1,I+1),1,COV(IJ),1)
               IJ = IJ + I
  100       CONTINUE
            CALL G02AAZ('U','N',IP,COV)
            CALL G02AAX('U',IP,COV,WK)
            CALL G02AAY('L','N',IP,WK)
            CALL G02AAX('L',IP,WK,COV)
            IF (RSS.LE.0.0D0) RSS = DDOT(N-IP,Q(IP+1,1),1,Q(IP+1,1),1)
         END IF
         IDF = N - IRANK
         IF (IDF.LE.0) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99995)
         ELSE
C
            RMS = RSS/DBLE(IDF)
            DO 120 I = 1, (IP*IP+IP)/2
               COV(I) = COV(I)*RMS
  120       CONTINUE
            DO 140 J = 1, IP
               SE(J) = SQRT(COV((J*J+J)/2))
  140       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99997 FORMAT (' ** On entry, LDQ.lt.IP : LDQ = ',I16,' IP = ',I16)
99996 FORMAT (' ** SVD solution failed to converge')
99995 FORMAT (' ** Degrees of freedom for error are less than or equal',
     *       ' to 0')
99994 FORMAT (' ** On entry, LDQ.lt.N and RSS.le.0.0 :',/'            ',
     *       '  LDQ = ',I16,' IP = ',I16)
99993 FORMAT (' ** On entry, TOL.lt.0.0: TOL = ',D13.5)
      END
