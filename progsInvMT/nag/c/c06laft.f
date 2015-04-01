      SUBROUTINE C06LAF(FUN,N,T,VALINV,ERREST,RELERR,ALPHAB,TFAC,MXTERM,
     *                  NTERMS,NA,ALOW,AHIGH,NFEVAL,WORK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Estimates values of the inverse Laplace transform of a given
C     function using a Fourier series approximation.
C     Real and imaginary parts of the function, and a bound on the
C     exponential order of the inverse, are required.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06LAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AHIGH, ALOW, ALPHAB, RELERR, TFAC
      INTEGER           IFAIL, MXTERM, N, NA, NFEVAL, NTERMS
C     .. Array Arguments ..
      DOUBLE PRECISION  ERREST(N), T(N), VALINV(N), WORK(4*MXTERM+2)
C     .. Subroutine Arguments ..
      EXTERNAL          FUN
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSFIN, ATEST, FINCHK, GTZERO, PI, RELERX, TAU,
     *                  TMAX
      INTEGER           I, IA1, IA2, IERR, IERR5, LINC, NCURR, NF, NREC
      LOGICAL           FIRST
C     .. Local Arrays ..
      DOUBLE PRECISION  A(2)
      INTEGER           IEV(2), IW(2)
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06LAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, MAX
C     .. Executable Statements ..
C
C     Check arguments, and initialize ERREST (also VALINV for safety,
C     although this is not essential).
C     IERR5 is used only for IFAIL = 5 which is an error return from
C     C06LAZ.
C
      IERR = 0
      IERR5 = 0
      IF (N.LT.1 .OR. MXTERM.LT.1 .OR. RELERR.LT.0.0D0 .OR. RELERR.GE.
     *    1.0D0 .OR. TFAC.LE.0.5D0) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) N, MXTERM, RELERR, TFAC
         NREC = 3
         GO TO 120
      ELSE IF (T(1).LT.0.0D0) THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) T(1)
         NREC = 1
         GO TO 120
      END IF
      VALINV(1) = 0.0D0
      ERREST(1) = 1.0D0
      DO 20 I = 2, N
         IF (T(I-1).GE.T(I)) GO TO 140
         VALINV(I) = 0.0D0
         ERREST(I) = 1.0D0
   20 CONTINUE
C
C     PI is required by C06LAZ.  The value of GTZERO is used only to
C     calculate TAU when T(N) is small, and is not critical.
C
      PI = X01AAF(0.0D0)
      GTZERO = 0.01D0
C
C     Initialize.
C     LINC is used to prevent looping if a has to
C     be increased successively.  NCURR is the current value of N.  IA1,
C     IA2 hold the index of the current low and high values of a,
C     respectively.  IW holds the start positions of the two halves of
C     the WORK array for use with C06LAZ.
C
      NFEVAL = 0
      NTERMS = 0
      NREC = 0
      NA = 2
      LINC = 0
      NCURR = N
      IW(1) = 1
      IW(2) = 2*MXTERM + 2
      IA1 = 1
      IA2 = 2
      RELERX = MAX(RELERR,100.0D0*X02AJF())
C
C     TAU and ALOW are not changed again in the algorithm.
C     ALOW holds the smallest value of a used, not the current low
C     value.
C
      TMAX = T(NCURR)
      IF (TMAX.EQ.0.0D0) TMAX = GTZERO
      TAU = TFAC*TMAX
      ALOW = ALPHAB - LOG(0.1D0*RELERX)/(2.0D0*TAU)
      AHIGH = ALOW + 1.0D0/TAU
      A(IA1) = ALOW
      A(IA2) = AHIGH
C
C     Test that exponentials can be taken in C06LAZ.
C
      ATEST = AHIGH*T(NCURR)
      IF (ATEST.LT.LOG(X02AMF()) .OR. ATEST.GT.-LOG(X02AMF())) THEN
         IERR = 3
         WRITE (P01REC,FMT=99997) T(N)
         NREC = 1
         GO TO 100
      END IF
C
C     These calls to C06LAZ are with FIRST = .TRUE. and use calls to FUN
C     only.  Subsequent calls, with FIRST = .FALSE., will use already
C     calculated values where possible.
C
      FIRST = .TRUE.
      CALL C06LAZ(T(NCURR),FUN,TAU,A(IA1),RELERX,FIRST,MXTERM,FINCHK,
     *            WORK(IW(IA1)),PI,IEV(IA1),NF,IERR5)
      NFEVAL = NFEVAL + NF
      NTERMS = MAX(NTERMS,IEV(IA1))
      IF (IERR5.GT.0) GO TO 100
C
   40 FIRST = .TRUE.
      CALL C06LAZ(T(NCURR),FUN,TAU,A(IA2),RELERX,FIRST,MXTERM,
     *            VALINV(NCURR),WORK(IW(IA2)),PI,IEV(IA2),NF,IERR5)
      NFEVAL = NFEVAL + NF
      NTERMS = MAX(NTERMS,IEV(IA2))
      IF (IERR5.GT.0) GO TO 100
C
C     Calculate absolute error estimate.  Convert to relative error if
C     function value is large enough.
C
      ERREST(NCURR) = ABS(VALINV(NCURR)-FINCHK)
      ABSFIN = ABS(VALINV(NCURR))
      IF (RELERX.LT.ABSFIN) ERREST(NCURR) = ERREST(NCURR)/ABSFIN
C
C     If RELERX is not achieved but 10.0*RELERX is achieved then
C     increase a.  An increment of 1.0/TAU should give about a factor
C     of 10 improvement in relative error.
C
      IF (RELERX.LT.ERREST(NCURR)) THEN
         IF (RELERX.GT.0.1D0*ERREST(NCURR) .AND. LINC.NE.NCURR) THEN
            NA = NA + 1
            LINC = NCURR
            IA1 = IA2
            IA2 = 3 - IA1
            AHIGH = A(IA1) + 1.0D0/TAU
            A(IA2) = AHIGH
C
C           If exponentials cannot be calculated at this stage,
C           give the IFAIL error that would have resulted if a had
C           not been changed.
C
            ATEST = AHIGH*T(NCURR)
            IF (ATEST.LT.LOG(X02AMF()) .OR. ATEST.GT.-LOG(X02AMF()))
     *          THEN
               IERR = 6
               IF (NCURR.LT.N) GO TO 60
               IERR = 4
               GO TO 100
            END IF
            GO TO 40
         ELSE
C
C           If this happens for first T value, then wrong value of
C           ALPHAB is a strong possibility, so give a more serious
C           IFAIL error.
C
            IERR = 6
            IF (NCURR.LT.N) GO TO 60
            IERR = 4
            GO TO 100
         END IF
      END IF
C
C     Results have been obtained for at least one value of T for the
C     current a values.
C     Attempt to repeat for the remaining values of T.
C
   60 DO 80 I = NCURR - 1, 1, -1
         CALL C06LAZ(T(I),FUN,TAU,A(IA1),RELERX,FIRST,MXTERM,FINCHK,
     *               WORK(IW(IA1)),PI,IEV(IA1),NF,IERR5)
         NFEVAL = NFEVAL + NF
         NTERMS = MAX(NTERMS,IEV(IA1))
         IF (IERR5.GT.0) GO TO 100
         CALL C06LAZ(T(I),FUN,TAU,A(IA2),RELERX,FIRST,MXTERM,VALINV(I),
     *               WORK(IW(IA2)),PI,IEV(IA2),NF,IERR5)
         NFEVAL = NFEVAL + NF
         NTERMS = MAX(NTERMS,IEV(IA2))
         IF (IERR5.GT.0) GO TO 100
         ERREST(I) = ABS(VALINV(I)-FINCHK)
         ABSFIN = ABS(VALINV(I))
         IF (RELERX.LT.ABSFIN) ERREST(I) = ERREST(I)/ABSFIN
         IF (RELERX.LT.ERREST(I)) THEN
            IF (RELERX.GT.0.1D0*ERREST(I)) THEN
               NA = NA + 1
               NCURR = I
               IA1 = IA2
               IA2 = 3 - IA1
               AHIGH = A(IA1) + 1.0D0/TAU
               A(IA2) = AHIGH
               ATEST = AHIGH*T(NCURR)
               IF (ATEST.LT.LOG(X02AMF()) .OR. ATEST.GT.-LOG(X02AMF()))
     *             THEN
                  IERR = 6
                  GO TO 80
               END IF
               GO TO 40
            ELSE
               IERR = 6
               GO TO 80
            END IF
         END IF
   80 CONTINUE
C
  100 IF (IERR5.GT.0) THEN
         IERR = 5
         WRITE (P01REC,FMT=99996) MXTERM
         NREC = 1
      END IF
  120 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
  140 WRITE (P01REC,FMT=99995)
      IERR = 2
      NREC = 1
      GO TO 120
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'         N =',I16,'  MXTERM =',I16,/'    R',
     *  'ELERR =',1P,D13.5,'       TFAC =',1P,D13.5)
99998 FORMAT (' ** On entry, T(1) is negative: T(1) =',1P,D13.5)
99997 FORMAT (' ** T(N) is too large: T(N) =',1P,D13.5)
99996 FORMAT (' ** Convergence failure in epsilon algorithm: MXTERM =',
     *  I8)
99995 FORMAT (' ** On entry, the elements of T are not strictly increa',
     *  'sing')
      END
