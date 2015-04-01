      SUBROUTINE G08AKF(N1,N2,TAIL,RANKS,U,P,WRK,LWRK,IWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1665 (JUN 1995).
C
C     THIS ROUTINE CALCULATES THE TAIL PROBABILITY P FOR
C     THE WILCOXON-MANN-WHINEY STATISTIC U FOR SAMPLE SIZES
C     N1 AND N2 FOR THE CASE OF TIES IN THE POOLED SAMPLE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AKF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, U
      INTEGER           IFAIL, LWRK, N1, N2
      CHARACTER*1       TAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  RANKS(N1+N2), WRK(LWRK)
      INTEGER           IWRK(2*(N1+N2+1))
C     .. Local Scalars ..
      INTEGER           I, IERROR, IF2, IV, N, NM, NREC, NSUM, NWRK
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AKZ, M01CBF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, NINT
C     .. Executable Statements ..
      NREC = 1
      IF (N1.LT.1 .OR. N2.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N1, N2
      ELSE IF (TAIL.NE.'T' .AND. TAIL.NE.'t' .AND. TAIL.NE.'U' .AND.
     *         TAIL.NE.'u' .AND. TAIL.NE.'L' .AND. TAIL.NE.'l') THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) TAIL
      ELSE IF (U.LT.0.0D0) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99997) U
      ELSE
         IERROR = 0
         N = MIN(N1,N2)
         NSUM = N1 + N2
         NWRK = N + N*(N+1)*NSUM - (N*(N+1)*(2*N+1))/3 + 1
         IF (LWRK.LT.NWRK) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99996) LWRK, NWRK
         ELSE
            DO 20 I = 1, NSUM
               IWRK(I) = NINT(2.0D0*RANKS(I))
   20       CONTINUE
            IF2 = 1
            CALL M01CBF(IWRK,1,NSUM,'A',IF2)
            NM = 2*N1*N2
            IV = NINT(2.0D0*U)
            IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
               CALL G08AKZ(N1,N2,IWRK,NSUM,IV,P,WRK,LWRK)
            ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
               IV = NM-IV
               CALL G08AKZ(N2,N1,IWRK,NSUM,IV,P,WRK,LWRK)
            ELSE IF (TAIL.EQ.'T' .OR. TAIL.EQ.'t') THEN
               IF ((2*IV).LE.NM) THEN
                  CALL G08AKZ(N1,N2,IWRK,NSUM,IV,P,WRK,LWRK)
               ELSE
                  IV = NM - IV
                  CALL G08AKZ(N2,N1,IWRK,NSUM,IV,P,WRK,LWRK)
               END IF
               P = 2.0D0*P
               P = MIN(P,1.0D0)
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (1X,'** On entry, N1.lt.1 or N2.lt.1 : N1 = ',I16,'  and ',
     *       'N2 = ',I16)
99998 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99997 FORMAT (1X,'** On entry, U.lt.0.0 : U = ',D13.5)
99996 FORMAT (1X,'** On entry, the workspace provided in WRK is not la',
     *       'rge enough.',/4X,'LWRK was ',I16,' but should be at leas',
     *       't ',I16)
      END
