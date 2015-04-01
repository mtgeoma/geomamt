      SUBROUTINE G10ABF(MODE,WEIGHT,N,X,Y,WT,RHO,YHAT,C,LDC,RSS,DF,RES,
     *                  H,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Fits a cubic smoothing spline
C     Option MODE allows different stages to be performed
C     Note: WK is used to store results between calls.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G10ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DF, RHO, RSS
      INTEGER           IFAIL, LDC, N
      CHARACTER         MODE, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), H(N), RES(N), WK(9*N+14), WT(*), X(N),
     *                  Y(N), YHAT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AVDF, AVH, P, Q
      INTEGER           I, IERROR
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G10ABW, G10ABX, G10ABY, G10ABZ
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
      ELSE IF (RHO.LT.0.0D0) THEN
         IERROR = 1
         WRITE (REC,FMT=99997) RHO
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         IERROR = 1
         WRITE (REC,FMT=99996) WEIGHT
      ELSE IF (MODE.NE.'Q' .AND. MODE.NE.'q' .AND. MODE.NE.'P' .AND.
     *         MODE.NE.'p' .AND. MODE.NE.'F' .AND. MODE.NE.'f') THEN
         IERROR = 1
         WRITE (REC,FMT=99995) MODE
      END IF
C
C     Initialise arrays for spline fitting
C
      IF (IERROR.EQ.0) THEN
         IF (MODE.NE.'Q' .AND. MODE.NE.'q') THEN
            CALL G10ABZ(WEIGHT,X,AVH,Y,WT,AVDF,N,WK(7*(N+2)+1),
     *                  WK(7*(N+2)+N+1),C,LDC,WK,WK(3*N+7),IERROR)
            IF (IERROR.EQ.0) THEN
               C(1,1) = AVH
               C(1,2) = AVDF
            ELSE IF (IERROR.EQ.2) THEN
               WRITE (REC,FMT=99994)
            ELSE IF (IERROR.EQ.3) THEN
               WRITE (REC,FMT=99993)
            END IF
         ELSE IF (MODE.EQ.'Q' .OR. MODE.EQ.'q') THEN
            AVH = C(1,1)
            AVDF = C(1,2)
         END IF
         IF (IERROR.EQ.0) THEN
C
C           Fit spline
C
            CALL G10ABY(X,AVH,WK(7*(N+2)+1),N,RHO,P,Q,RSS,DF,WK(7*(N+2)
     *                  +N+1),C,LDC,WK,WK(3*N+7),WK(5*N+11),RES)
C
            DO 20 I = 1, N
               YHAT(I) = Y(I) - WK(7*(N+2)+I)*RES(I)
               RES(I) = RES(I)/AVDF
   20       CONTINUE
            RSS = RSS/(AVDF*AVDF)
C
C           Calculate influences
C
            CALL G10ABW(X,AVH,N,WK,P,H,WK(7*(N+2)+1))
C
C           Calculate coefficients
C
            IF (MODE.EQ.'F' .OR. MODE.EQ.'f') THEN
               CALL G10ABX(X,AVH,Y,WT,N,Q,YHAT,C,LDC,WK(5*N+11),
     *                     WK(6*N+13))
            END IF
         END IF
      END IF
C
C     Check for error condition
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.3: N = ',I16)
99998 FORMAT (1X,'** On entry, LDC.lt.N-1: LDC = ',I16,' N = ',I16)
99997 FORMAT (1X,'** On entry, RHO.lt.0.0: RHO = ',D13.5)
99996 FORMAT (1X,'** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99995 FORMAT (1X,'** On entry, MODE is not valid: MODE = ',A1)
99994 FORMAT (1X,'** On entry, at least one element of WT.lt.0.0')
99993 FORMAT (1X,'** On entry, X is not a strictly ordered array')
      END
