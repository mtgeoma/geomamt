      SUBROUTINE F01ZDF(JOB,M,N,KL,KU,A,LDA,B,LDB,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Converts a complex banded matrix from packed storage to unpacked,
C     or vice-versa.
C
C     Parameters:
C
C     JOB - CHARACTER*1.
C       JOB = 'P' or 'p' : 'Pack' - the contents of matrix A are to be
C       packed into matrix B.
C       JOB = 'U' or 'u' : 'Unpack' - the contents of matrix B are to
C       be unpacked into matrix A.
C
C     M - INTEGER.
C       On entry, M specifies the number of rows of the matrix A.
C
C     N - INTEGER.
C       On entry, N specifies the number of columns of the matrix A.
C
C     KL - INTEGER.
C       On entry, KL specifies the number of sub-diagonals of the
C       matrix A. KL >= 0.
C
C     KU - INTEGER.
C       On entry, KU specifies the number of super-diagonals of the
C       matrix A. KU >= 0.
C
C     A(LDA, N) - complex array.
C       If JOB = 'P' or 'p', then on entry the leading M by N part
C       of array A must contain the band matrix stored in unpacked
C       form. Elements of A that lie above the KUth super-diagonal
C       and below the KLth sub-diagonal are not referenced and need
C       not be defined.
C       If JOB = 'U' or 'u', then on entry the array A may be
C       undefined. On exit, A contains the band matrix stored in
C       unpacked form. Elements of the leading M by N part of A that
C       lie above the KUth super-diagonal or below the KLth
C       sub-diagonal are assigned the value zero.
C
C     LDA - INTEGER.
C       On entry, LDA specifies the first dimension of A as declared
C       in the calling (sub)program. LDA >= M.
C
C     B(LDA, MIN(M+KU,N)) - complex array.
C       If JOB = 'U' or 'u', then on entry the leading (KL + KU + 1)
C       by min(M+KU,N) part of array B must contain the band matrix,
C       packed column by column, with the leading diagonal of the
C       matrix in row (KU + 1) of the array, the first super-diagonal
C       starting at position 2 in row KU, the first sub-diagonal
C       starting at position 1 in row (KU + 2), and so on.
C       Elements in the array B that do not correspond to elements
C       in the band matrix (such as the top left KU by KU triangle)
C       are not referenced.
C       If JOB = 'P', then on entry B may be undefined, and on exit
C       B contains the band matrix, packed as described above.
C       Elements in the array B that do not correspond to elements
C       in the band matrix (such as the top left KU by KU triangle)
C       are unreferenced.
C
C     LDB - INTEGER.
C       On entry, LDB specifies the first dimension of B as declared
C       in the calling (sub)program. LDB >= (KL + KU + 1).
C
C     IFAIL - INTEGER.
C       Error flag. On exit,
C         IFAIL = 1 if JOB <> 'P', 'p', 'U' or 'u'.
C         IFAIL = 2 if KL < 0.
C         IFAIL = 3 if KU < 0.
C         IFAIL = 4 if LDA < M.
C         IFAIL = 5 if LDB < KL + KU + 1.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01ZDF')
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      INTEGER           IFAIL, KL, KU, LDA, LDB, M, N
      CHARACTER         JOB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,N), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, K, NREC
      LOGICAL           PACK
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      IERR = 0
      PACK = JOB .EQ. 'P' .OR. JOB .EQ. 'p'
      IF ( .NOT. PACK .AND. JOB.NE.'U' .AND. JOB.NE.'u') THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) JOB
      ELSE IF (KL.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (REC,FMT=99998) KL
      ELSE IF (KU.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (REC,FMT=99997) KU
      ELSE IF (LDA.LT.M) THEN
         IERR = 4
         NREC = 1
         WRITE (REC,FMT=99996) LDA, M
      ELSE IF (LDB.LT.KL+KU+1) THEN
         IERR = 5
         NREC = 2
         WRITE (REC,FMT=99995) LDB, KL, KU
      ELSE
         IF (PACK) THEN
C         Convert from unpacked form to packed form.
            DO 40 J = 1, N
               K = KU + 1 - J
               DO 20 I = MAX(1,J-KU), MIN(M,J+KL)
                  B(K+I,J) = A(I,J)
   20          CONTINUE
   40       CONTINUE
         ELSE
C         Convert from packed form to unpacked form.
            DO 120 J = 1, N
               K = KU + 1 - J
               DO 60 I = MAX(1,J-KU), MIN(M,J+KL)
                  A(I,J) = B(K+I,J)
   60          CONTINUE
C            Zero-fill unused elements.
               DO 80 I = 1, J - KU - 1
                  A(I,J) = ZERO
   80          CONTINUE
               DO 100 I = J + KL + 1, M
                  A(I,J) = ZERO
  100          CONTINUE
  120       CONTINUE
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, JOB is not equal to ''P'' or ''U'': JOB',
     *       ' = ''',A,'''')
99998 FORMAT (1X,'** On entry, KL.lt.0: KL =',I16)
99997 FORMAT (1X,'** On entry, KU.lt.0: KU =',I16)
99996 FORMAT (1X,'** On entry, LDA.lt.M: LDA = ',I16,', M = ',I16)
99995 FORMAT (1X,'** On entry, LDB.lt.KL+KU+1: LDB =',I16,',',/4X,
     *       'KL =',I16,', KU =',I16)
      END
