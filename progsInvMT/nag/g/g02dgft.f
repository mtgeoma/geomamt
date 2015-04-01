      SUBROUTINE G02DGF(WEIGHT,N,WT,RSS,IP,IRANK,COV,Q,LDQ,SVD,P,Y,B,SE,
     *                  RES,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES A REGRESSION FOR NEW Y VARIABLE
C     G02DGF SHOULD BE CALLED AFTER G02DAF
C     NOTE: THE FIRST COLUMN OF Q IS OVERWRITTEN
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           IFAIL, IP, IRANK, LDQ, N
      LOGICAL           SVD
      CHARACTER*1       WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV((IP*IP+IP)/2), P(*), Q(LDQ,IP+1),
     *                  RES(N), SE(IP), WK(5*(IP-1)+IP*IP), WT(*), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  RCF
      INTEGER           I, IERROR, IFAULT, IP2, J, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP(1)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QDF, DSCAL, F06FBF, DCOPY, DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (IP.LT.1) THEN
         WRITE (P01REC,FMT=99999) IP
      ELSE IF (N.LT.IP) THEN
         WRITE (P01REC,FMT=99998) N, IP
      ELSE IF (IRANK.LE.0) THEN
         WRITE (P01REC,FMT=99991) IRANK
      ELSE IF ( .NOT. SVD .AND. IRANK.NE.IP) THEN
         NREC = 2
         WRITE (P01REC,FMT=99992) IRANK, IP
      ELSE IF (SVD .AND. IRANK.GT.IP) THEN
         NREC = 2
         WRITE (P01REC,FMT=99993) IRANK, IP
      ELSE IF (LDQ.LT.N) THEN
         WRITE (P01REC,FMT=99996) LDQ, N
      ELSE IF (RSS.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99997) RSS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            DO 20 I = 1, N
               IF (WT(I).LT.0.0D0) THEN
                  GO TO 40
C
               ELSE
                  Q(I,1) = SQRT(WT(I))*Y(I)
               END IF
   20       CONTINUE
            GO TO 60
C
   40       IERROR = 2
            WRITE (P01REC(1),FMT=99995) I
            GO TO 100
C
         ELSE IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
            CALL DCOPY(N,Y,1,Q,1)
         ELSE
            IERROR = 1
            WRITE (P01REC(1),FMT=99994) WEIGHT
            GO TO 100
C
         END IF
   60    IFAULT = 1
         CALL F01QDF('T','S',N,IP,Q(1,2),LDQ,P,1,Q,N,WKSP,IFAULT)
         IF (SVD) THEN
            CALL DGEMV('T',IP,IP,1.0D0,WK,IP,Q,1,0.0D0,SE,1)
            IP2 = IP*IP
            CALL DGEMV('T',IRANK,IP,1.0D0,P(2*IP+1),IP,SE,1,0.0D0,B,1)
            CALL DGEMV('N',IP,IP-IRANK,1.0D0,WK(IRANK*IP+1),IP,
     *                 SE(IRANK+1),1,0.0D0,RES,1)
            CALL DCOPY(N-IP,Q(IP+1,1),1,RES(IP+1),1)
            IFAULT = 1
            CALL F01QDF('N','S',N,IP,Q(1,2),LDQ,P,1,RES,N,WKSP,IFAULT)
         ELSE
            CALL DCOPY(IP,Q,1,B,1)
            CALL DTRSV('U','N','N',IP,Q(1,2),LDQ,B,1)
            CALL F06FBF(IP,0.0D0,RES,1)
            CALL DCOPY(N-IP,Q(IP+1,1),1,RES(IP+1),1)
            IFAULT = 1
            CALL F01QDF('N','S',N,IP,Q(1,2),LDQ,P,1,RES,N,WKSP,IFAULT)
         END IF
         RCF = RSS
         RSS = DDOT(N,RES,1,RES,1)
         RCF = RSS/RCF
         CALL DSCAL((IP*IP+IP)/2,RCF,COV,1)
         DO 80 J = 1, IP
            SE(J) = SQRT(COV((J*J+J)/2))
   80    CONTINUE
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99998 FORMAT (' ** On entry, N.lt.IP : N = ',I16,' IP = ',I16)
99997 FORMAT (' ** On entry, RSS.le.0.0 : RSS = ',D13.5)
99996 FORMAT (' ** On entry, LDQ.lt.N : LDQ = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99994 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99993 FORMAT (' ** On entry, IRANK.gt.IP and SVD is .TRUE. :',/'      ',
     *       '        IRANK = ',I16,' IP = ',I16)
99992 FORMAT (' ** On entry, IRANK.ne.IP and SVD is .FALSE. :',/'     ',
     *       '         IRANK = ',I16,' IP = ',I16)
99991 FORMAT (' ** On entry, IRANK.le.0: IRANK = ',I16)
      END
