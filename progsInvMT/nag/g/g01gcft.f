      DOUBLE PRECISION FUNCTION G01GCF(X,DF,RLAMDA,TOL,MAXIT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1032 (JUN 1993).
C
C     Computes the lower tail probability of X for the non-central
C     CHI-Squared distribution with degrees of freedom DF.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01GCF')
      DOUBLE PRECISION                 ZERO, ONE, TWO
      PARAMETER                        (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, RLAMDA, TOL, X
      INTEGER                          IFAIL, MAXIT
C     .. Local Scalars ..
      DOUBLE PRECISION                 A2, A2J, ALGF, ALGF0, ALGM, C2,
     *                                 GM, GM0, LP, P, PJ, PJGM, PREC,
     *                                 PSUM, Q, SUM, UFLO, X2
      INTEGER                          IERR, IFAULT, IMAXIT, J, M
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S14ABF, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         S14ABF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G01BKF, S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE, EXP, LOG
C     .. Executable Statements ..
      G01GCF = ZERO
      SUM = ZERO
      IERR = 1
      IF (DF.LT.ZERO) THEN
         WRITE (REC,FMT=99999) DF
      ELSE IF (RLAMDA.LT.ZERO) THEN
         WRITE (REC,FMT=99998) RLAMDA
      ELSE IF (X.LT.ZERO) THEN
         WRITE (REC,FMT=99997) X
      ELSE IF (RLAMDA.EQ.ZERO .AND. DF.EQ.ZERO) THEN
         WRITE (REC,FMT=99991)
      ELSE IF (MAXIT.LT.1) THEN
         WRITE (REC,FMT=99993) MAXIT
      ELSE
         IERR = 0
      END IF
      IF (IERR.EQ.0) THEN
         PREC = X02AJF()*10.0D0
         IF (TOL.GT.PREC .AND. TOL.LT.1.0D0) PREC = TOL
         A2 = DF/TWO
         X2 = X/TWO
         IF (RLAMDA.EQ.ZERO) THEN
            IFAULT = 1
            CALL S14BAF(A2,X2,PREC,P,Q,IFAULT)
            IF (IFAULT.NE.0) THEN
               WRITE (REC,FMT=99994)
               IERR = 5
               GO TO 80
            END IF
            SUM = P
         ELSE IF (X.NE.ZERO) THEN
            C2 = RLAMDA/TWO
            M = C2
            UFLO = LOG(X02AMF())
C
C           Set up values for term with highest Poisson weight
C
            IF (A2.EQ.ZERO) THEN
               PSUM = ONE
               IF (M.EQ.0) M = 1
            ELSE
               IFAULT = 1
               CALL S14BAF(A2,X2,PREC,PSUM,Q,IFAULT)
               IF (IFAULT.NE.0) THEN
                  WRITE (REC,FMT=99994)
                  IERR = 5
                  GO TO 80
               END IF
            END IF
            IFAULT = 1
            LP = DBLE(M)*LOG(C2) - S14ABF(DBLE(M+1),IFAULT) - C2
            IF (LP.LT.UFLO .OR. IFAULT.EQ.2) THEN
               IERR = 2
               WRITE (REC,FMT=99996)
               GO TO 80
            END IF
            P = EXP(LP)
            PJ = P
            J = M
            A2J = A2 + M
            IFAULT = 1
            CALL S14BAF(A2J,X2,PREC,GM0,Q,IFAULT)
            GM = GM0
            IF (IFAULT.NE.0) THEN
               WRITE (REC,FMT=99994)
               IERR = 5
               GO TO 80
            END IF
            ALGF0 = S14ABF(A2J+ONE,IFAULT)
            ALGF = ALGF0
            PJGM = PJ*GM
C
C           Sum down from highest Poisson Weight
C
   20       SUM = SUM + PJGM
            J = J - 1
            IF (J.GE.0) THEN
               ALGF = ALGF - LOG(A2J)
               PJ = PJ/(C2/DBLE(J+1))
               A2J = A2J - ONE
               IFAULT = 1
               ALGM = A2J*LOG(X2) - X2 - ALGF
               IF (ALGM.GT.UFLO) THEN
                  IF (IFAULT.EQ.2) THEN
                     IERR = 4
                     WRITE (REC,FMT=99992)
                     GO TO 80
                  ELSE
                     GM = GM + EXP(ALGM)
                  END IF
               END IF
               PJGM = PJ*GM
               IF (PJ*PSUM.GT.PREC*SUM) GO TO 20
            END IF
C
C           Restart from highest Poisson weight
C
            PJ = P
            IFAULT = 1
            CALL G01BKF(C2,M,Q,PSUM,GM,IFAULT)
            GM = GM0
            A2J = A2 + M
            ALGF = ALGF0
            IMAXIT = MAXIT + M - J
            J = M
C
C           Sum up from highest Poisson weight
C
   40       J = J + 1
            IF (J.GT.IMAXIT) THEN
               IERR = 3
               WRITE (REC,FMT=99995) MAXIT
               GO TO 60
            END IF
            PJ = PJ*(C2/DBLE(J))
            PSUM = PSUM - PJ
            IF (PSUM*GM.GT.PREC*SUM) THEN
               A2J = A2J + ONE
               IFAULT = 1
               ALGM = (A2J-ONE)*LOG(X2) - X2 - ALGF
               IF (ALGM.GT.UFLO) THEN
                  IF (IFAULT.EQ.2) THEN
                     IERR = 4
                     WRITE (REC,FMT=99992)
                     GO TO 80
                  ELSE
                     GM = GM - EXP(ALGM)
                  END IF
               END IF
               ALGF = ALGF + LOG(A2J)
               PJGM = PJ*GM
               SUM = SUM + PJGM
               GO TO 40
            END IF
         END IF
      END IF
   60 G01GCF = SUM
   80 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, DF.lt.0.0 : DF = ',1P,D13.5)
99998 FORMAT (1X,'** On entry, RLAMDA.lt.0.0 : RLAMDA = ',1P,D13.5)
99997 FORMAT (1X,'** On entry, X.lt.0.0 : X = ',1P,D13.5)
99996 FORMAT (1X,'** Poisson weight too small to calculate.')
99995 FORMAT (1X,'** Convergence has not been achieved in MAXIT steps:',
     *       ' MAXIT = ',I16)
99994 FORMAT (1X,'** Central Chi squared function fails to converge.')
99993 FORMAT (1X,'** On entry, MAXIT.lt.1 : MAXIT = ',I16)
99992 FORMAT (1X,'** The values of X and RLAMDA are too large for accu',
     *       'rate evaluation.')
99991 FORMAT (1X,'** On entry, DF.eq.0.0 and RLAMDA.eq.0.0')
      END
