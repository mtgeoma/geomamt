      SUBROUTINE G10ACF(METHOD,WEIGHT,N,X,Y,WT,YHAT,C,LDC,RSS,RDF,RES,H,
     *                  CRIT,RHO,U,TOL,MAXCAL,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     H IS USED TO STORE WWT UNTIL G10ABW
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G10ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CRIT, RDF, RHO, RSS, TOL, U
      INTEGER           IFAIL, LDC, MAXCAL, N
      CHARACTER         METHOD, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), H(N), RES(N), WK(7*(N+2)), WT(*),
     *                  X(N), Y(N), YHAT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AVH, AVRDF, P, Q, UX, WTSUM
      INTEGER           I, IERROR, IND
C     .. Local Arrays ..
      DOUBLE PRECISION  D(17)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G10ABW, G10ABX, G10ABY, G10ABZ, G10ACX, G10ACZ
C     .. Executable Statements ..
C
C     Input checking
C
      IERROR = 0
      IF (N.LT.3) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (LDC.LT.N-1) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) LDC, N
      ELSE IF (METHOD.NE.'D' .AND. METHOD.NE.'d' .AND. METHOD.NE.
     *         'C' .AND. METHOD.NE.'c' .AND. METHOD.NE.'G' .AND.
     *         METHOD.NE.'g') THEN
         IERROR = 1
         WRITE (REC,FMT=99997) METHOD
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         IERROR = 1
         WRITE (REC,FMT=99996) WEIGHT
      ELSE IF ((METHOD.EQ.'D' .OR. METHOD.EQ.'d') .AND. CRIT.LE.2.0D0)
     *         THEN
         IERROR = 1
         WRITE (REC,FMT=99995) CRIT
      ELSE IF ((METHOD.EQ.'D' .OR. METHOD.EQ.'d') .AND. CRIT.GT.N) THEN
         IERROR = 1
         WRITE (REC,FMT=99994) CRIT
      ELSE IF (TOL.LT.X02AJF()) THEN
         IERROR = 1
         WRITE (REC,FMT=99993) TOL
      ELSE IF (U.LE.TOL) THEN
         IERROR = 1
         WRITE (REC,FMT=99992) U, TOL
      ELSE IF (MAXCAL.LT.3) THEN
         IERROR = 1
         WRITE (REC,FMT=99991) MAXCAL
      ELSE
C
         CALL G10ABZ(WEIGHT,X,AVH,Y,WT,AVRDF,N,H,YHAT,C,LDC,WK,WK(3*N+7)
     *               ,IERROR)
         IF (IERROR.EQ.2) THEN
            WRITE (REC,FMT=99990)
         ELSE IF (IERROR.EQ.3) THEN
            WRITE (REC,FMT=99989)
         ELSE
C
            IF (METHOD.EQ.'D' .OR. METHOD.EQ.'d') THEN
C
C
               CALL G10ACX(U,TOL,MAXCAL,CRIT,RHO,X,AVH,H,N,P,Q,YHAT,C,
     *                     LDC,WK,RES,D,IERROR)
               IF (IERROR.EQ.4) THEN
                  WRITE (REC,FMT=99988) U
                  GO TO 60
               ELSE IF (IERROR.EQ.5) THEN
                  WRITE (REC,FMT=99987) TOL
               ELSE IF (IERROR.EQ.6) THEN
                  WRITE (REC,FMT=99986)
               END IF
            ELSE
               IF (METHOD.EQ.'G' .OR. METHOD.EQ.'g') THEN
                  IND = 0
               ELSE
                  IND = -1
               END IF
C
C              Find local minimum
C
               UX = U
               CALL G10ACZ(U,TOL,MAXCAL,IND,RHO,CRIT,X,AVH,H,N,P,Q,YHAT,
     *                     C,LDC,WK,WK(3*N+7),WK(5*N+11),RES,IERROR)
               IF (IERROR.EQ.6) THEN
                  WRITE (REC,FMT=99986)
               ELSE IF (IERROR.EQ.7) THEN
                  WRITE (REC,FMT=99985) UX
               END IF
            END IF
C
C           Calculate spline coefficients
C
            CALL G10ABY(X,AVH,H,N,RHO,P,Q,RSS,RDF,YHAT,C,LDC,WK,
     *                  WK(3*N+7),WK(5*N+11),RES)
C
            WTSUM = 0.0D0
            DO 20 I = 1, N
               YHAT(I) = Y(I) - H(I)*RES(I)
               RES(I) = RES(I)/AVRDF
               WTSUM = WTSUM + WT(I)
   20       CONTINUE
            RSS = RSS/(AVRDF*AVRDF)
            IF (METHOD.NE.'D' .AND. METHOD.NE.'d') THEN
               IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
                  WTSUM = 0.0D0
                  DO 40 I = 1, N
                     WTSUM = WTSUM + WT(I)
   40             CONTINUE
               ELSE
                  WTSUM = N
               END IF
               CRIT = CRIT/(AVRDF*AVRDF)
               CRIT = CRIT/WTSUM
            END IF
C
C
C           Calculate coefficients
C
            CALL G10ABX(X,AVH,Y,H,N,Q,YHAT,C,LDC,WK(5*N+11),WK(6*N+13))
C
C           Calculate influences
C
            CALL G10ABW(X,AVH,N,WK,P,H,H)
C
         END IF
      END IF
C
C            Check for error condition
C
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.3: N = ',I16)
99998 FORMAT (1X,'** On entry, LDC.lt.N-1: LDC = ',I16,' N = ',I16)
99997 FORMAT (1X,'** On entry, METHOD is not valid: METHOD = ',A1)
99996 FORMAT (1X,'** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99995 FORMAT (1X,'** On entry, METHOD = D and CRIT.le.2.0: CRIT = ',
     *       D13.5)
99994 FORMAT (1X,'** On entry, METHOD = D and CRIT.gt.N: CRIT = ',D13.5)
99993 FORMAT (1X,'** On entry, TOL.lt.machine precision: TOL = ',D13.5)
99992 FORMAT (1X,'** On entry, U.le.TOL: U = ',D13.5,' TOL = ',D13.5)
99991 FORMAT (1X,'** On entry, MAXCAL.lt.3: MAXCAL = ',I16)
99990 FORMAT (1X,'** On entry, at least one element of WT.le.0.0')
99989 FORMAT (1X,'** On entry, X is not a strictly ordered array')
99988 FORMAT (1X,'** For the specified degrees of freedom, RHO.gt.U: U',
     *       ' = ',D13.5)
99987 FORMAT (1X,'** Accuracy of TOL cannot be achieved: TOL = ',D13.5)
99986 FORMAT (1X,'** MAXCAL iterations have been performed')
99985 FORMAT (1X,'** Optimum value of RHO lies above U: U = ',D13.5)
      END
