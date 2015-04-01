      SUBROUTINE G03FAF(ROOTS,N,D,NDIM,X,LDX,EVAL,WK,IWK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Finds principal coordinates of distance
C     matrix D. Returned in X
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03FAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDX, N, NDIM
      CHARACTER         ROOTS
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N*(N-1)/2), EVAL(N), WK(N*(N+17)/2-1),
     *                  X(LDX,NDIM)
      INTEGER           IWK(5*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  R, TOTAL
      INTEGER           I, IERROR, IJ, IJI, J, M, NN, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06EDF, G03FAY, G03FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (NDIM.LT.1) THEN
         WRITE (P01REC,FMT=99999) NDIM
      ELSE IF (N.LE.NDIM) THEN
         WRITE (P01REC,FMT=99998) N, NDIM
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC,FMT=99997) LDX, N
      ELSE IF (ROOTS.NE.'A' .AND. ROOTS.NE.'a' .AND. ROOTS.NE.'L' .AND.
     *         ROOTS.NE.'l') THEN
         WRITE (P01REC,FMT=99996) ROOTS
      ELSE
         IERROR = 0
C
C          Form E = -1/2d**2 in WK
C
         IJ = 0
         IJI = 0
         DO 40 I = 1, N
            DO 20 J = 1, I - 1
               IJ = IJ + 1
               IF (D(IJ).LT.0.0D0) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99995) IJ
                  GO TO 100
               END IF
               IJI = IJI + 1
               WK(IJI) = -0.5D0*D(IJ)*D(IJ)
   20       CONTINUE
            IJI = IJI + 1
            WK(IJI) = 0.0D0
   40    CONTINUE
         NN = IJI
C
C        Standardise E to give F
C
         CALL G03FAZ(N,WK,WK(NN+1),TOTAL)
         IF (TOTAL.EQ.0.0D0) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99994)
            GO TO 100
         END IF
C
C           Find eigenvalues and eigenvectors of Q
C
         CALL G03FAY(ROOTS,WK,N,NDIM,EVAL,M,WK(NN+1),WK(NN+N+1),
     *               WK(NN+2*N),IWK,IWK(N+1),X,LDX,WK(NN+3*N-1),
     *               IWK(2*N+1),IERROR)
C
C        Scale eigenvectors
C
         IF (IERROR.EQ.0) THEN
            DO 60 I = 1, NDIM
               IF (EVAL(I).LE.0.0D0) THEN
                  IERROR = 3
                  WRITE (P01REC,FMT=99993)
                  GO TO 80
               ELSE
                  R = SQRT(EVAL(I))
                  CALL F06EDF(N,R,X(1,I),1)
               END IF
   60       CONTINUE
   80       CONTINUE
C
C           Normalize eigenvalues
C
            R = -1.0D0/(DBLE(N)*TOTAL)
            CALL F06EDF(M,R,EVAL,1)
         ELSE
            IERROR = 4
            WRITE (P01REC,FMT=99992)
         END IF
      END IF
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NDIM .lt. 1 : NDIM = ',I16)
99998 FORMAT (' ** On entry, N .le. NDIM: N = ',I16,' NDIM = ',I16)
99997 FORMAT (' ** On entry, LDX .lt. N: LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, ROOTS is not valid: ROOTS = ',A1)
99995 FORMAT (' ** On entry, the ',I16,' th element of D .lt. 0.0')
99994 FORMAT (' ** On entry, all the elements of D = 0.0')
99993 FORMAT (' ** There are less than NDIM eigenvalues .gt. 0.0')
99992 FORMAT (' ** The computation of the eigenvalues or eigenvectors ',
     *       'has failed')
      END
