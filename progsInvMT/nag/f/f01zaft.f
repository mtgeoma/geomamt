      SUBROUTINE F01ZAF(JOB,UPLO,DIAG,N,A,LDA,B,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Copies a real triangular matrix from packed vector storage to
C     a two-dimensional array, or vice-versa.
C
C     Parameters:
C
C     JOB - CHARACTER*1.
C       JOB = 'P' or 'p' : 'Pack' - the contents of matrix A are to be
C       packed into vector B.
C       JOB = 'U' or 'u' : 'Unpack' - the contents of vector B are to
C       be unpacked into matrix A.
C
C     UPLO - CHARACTER*1.
C       UPLO = 'U' or 'u' : 'Upper' - the matrix is upper triangular,
C       and stored by column in the packed vector.
C       UPLO = 'L' or 'l' : 'Upper' - the matrix is lower triangular,
C       and stored by column in the packed vector.
C
C     DIAG - CHARACTER*1.
C       DIAG = 'N' or 'n' : 'Non-unit' - the triangular matrix has
C       non-unit diagonal elements, which are packed or unpacked with
C       the rest of the matrix.
C       DIAG = 'U' or 'u' : 'Unit' - the triangular matrix is assumed
C       to have unit diagonal elements which are not stored, in the
C       format being copied from, but are inserted into the format
C       being copied to.
C       DIAG = ' ' : 'Blank' - the matrix is assumed to be strictly
C       upper or strictly lower triangular. No reference is made
C       either to the diagonal elements of matrix A or to the
C       corresponding elements of vector B, although space must be
C       allocated for them in the packed vector.
C
C     N - INTEGER.
C       The order of the matrix. N > 0.
C
C     A(LDA, N) - real array.
C       If JOB = 'P' or 'p', then on entry A must contain the matrix
C       stored in the upper triangle if UPLO = 'U' or 'u', and stored
C       in the lower triangle if UPLO = 'L' or 'l'. The opposite
C       triangle of A may be undefined.
C       If JOB = 'U' or 'u', then on exit A will contain the matrix
C       stored in the upper triangle if UPLO = 'U' or 'u', and stored
C       in the lower triangle if UPLO = 'L' or 'l'. The opposite
C       triangle is left untouched.
C
C     LDA - INTEGER.
C       The leading dimension of array A as declared in the calling
C       (sub)program.
C
C     B((N*(N+1))/2) - real array.
C       If JOB = 'U' or 'u', then on entry B must contain the triangular
C       matrix packed by column.
C       If JOB = 'P' or 'p', then on exit B will contain the triangular
C       matrix packed by column.
C
C     IFAIL - INTEGER.
C       Error flag. On exit,
C         IFAIL = 1 if JOB <> 'P', 'p', 'U' or 'u'.
C         IFAIL = 2 if UPLO <> 'U', 'u', 'L' or 'l'.
C         IFAIL = 3 if DIAG <> 'N', 'n', 'U', 'u', 'B' or 'b'.
C         IFAIL = 4 if N < 1.
C         IFAIL = 5 if LDA < N.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01ZAF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, N
      CHARACTER         DIAG, JOB, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), B((N*(N+1))/2)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, K, NREC
      LOGICAL           LOWER, PACK, UNIT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IERR = 0
      PACK = JOB .EQ. 'P' .OR. JOB .EQ. 'p'
      LOWER = UPLO .EQ. 'L' .OR. UPLO .EQ. 'l'
      UNIT = DIAG .EQ. 'U' .OR. DIAG .EQ. 'u'
      IF ( .NOT. PACK .AND. JOB.NE.'U' .AND. JOB.NE.'u') THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) JOB
      ELSE IF ( .NOT. LOWER .AND. UPLO.NE.'U' .AND. UPLO.NE.'u') THEN
         IERR = 2
         NREC = 1
         WRITE (REC,FMT=99998) UPLO
      ELSE IF ( .NOT. UNIT .AND. DIAG.NE.'N' .AND. DIAG.NE.'n' .AND.
     *         DIAG.NE.'B' .AND. DIAG.NE.'b') THEN
         IERR = 3
         NREC = 1
         WRITE (REC,FMT=99997) DIAG
      ELSE IF (N.LT.1) THEN
         IERR = 4
         NREC = 1
         WRITE (REC,FMT=99996) N
      ELSE IF (LDA.LT.N) THEN
         IERR = 5
         NREC = 1
         WRITE (REC,FMT=99995) LDA, N
      END IF
C
      IF (IERR.EQ.0) THEN
         IF (PACK) THEN
C           Pack the triangular matrix from A into B.
            IF (LOWER) THEN
C              Lower triangular matrix.
               IF (UNIT) THEN
                  K = 0
                  DO 40 J = 1, N
                     K = K + 1
                     B(K) = ONE
                     DO 20 I = J + 1, N
                        K = K + 1
                        B(K) = A(I,J)
   20                CONTINUE
   40             CONTINUE
               ELSE IF ((DIAG.EQ.'B') .OR. (DIAG.EQ.'b')) THEN
                  K = 0
                  DO 80 J = 1, N
                     K = K + 1
                     DO 60 I = J + 1, N
                        K = K + 1
                        B(K) = A(I,J)
   60                CONTINUE
   80             CONTINUE
               ELSE
                  K = 0
                  DO 120 J = 1, N
                     DO 100 I = J, N
                        K = K + 1
                        B(K) = A(I,J)
  100                CONTINUE
  120             CONTINUE
               END IF
            ELSE
C              Upper triangular matrix.
               IF (UNIT) THEN
                  K = 0
                  DO 160 J = 1, N
                     DO 140 I = 1, J - 1
                        K = K + 1
                        B(K) = A(I,J)
  140                CONTINUE
                     K = K + 1
                     B(K) = ONE
  160             CONTINUE
               ELSE IF ((DIAG.EQ.'B') .OR. (DIAG.EQ.'b')) THEN
                  K = 0
                  DO 200 J = 1, N
                     DO 180 I = 1, J - 1
                        K = K + 1
                        B(K) = A(I,J)
  180                CONTINUE
                     K = K + 1
  200             CONTINUE
               ELSE
                  K = 0
                  DO 240 J = 1, N
                     DO 220 I = 1, J
                        K = K + 1
                        B(K) = A(I,J)
  220                CONTINUE
  240             CONTINUE
               END IF
            END IF
         ELSE
C           Unpack the triangular matrix from B into A.
            IF (LOWER) THEN
C              Lower triangular matrix.
               IF (UNIT) THEN
                  K = 0
                  DO 280 J = 1, N
                     K = K + 1
                     A(J,J) = ONE
                     DO 260 I = J + 1, N
                        K = K + 1
                        A(I,J) = B(K)
  260                CONTINUE
  280             CONTINUE
               ELSE IF ((DIAG.EQ.'B') .OR. (DIAG.EQ.'b')) THEN
                  K = 0
                  DO 320 J = 1, N
                     K = K + 1
                     DO 300 I = J + 1, N
                        K = K + 1
                        A(I,J) = B(K)
  300                CONTINUE
  320             CONTINUE
               ELSE
                  K = 0
                  DO 360 J = 1, N
                     DO 340 I = J, N
                        K = K + 1
                        A(I,J) = B(K)
  340                CONTINUE
  360             CONTINUE
               END IF
            ELSE
C              Upper triangular matrix.
               IF (UNIT) THEN
                  K = 0
                  DO 400 J = 1, N
                     DO 380 I = 1, J - 1
                        K = K + 1
                        A(I,J) = B(K)
  380                CONTINUE
                     K = K + 1
                     A(J,J) = ONE
  400             CONTINUE
               ELSE IF ((DIAG.EQ.'B') .OR. (DIAG.EQ.'b')) THEN
                  K = 0
                  DO 440 J = 1, N
                     DO 420 I = 1, J - 1
                        K = K + 1
                        A(I,J) = B(K)
  420                CONTINUE
                     K = K + 1
  440             CONTINUE
               ELSE
                  K = 0
                  DO 480 J = 1, N
                     DO 460 I = 1, J
                        K = K + 1
                        A(I,J) = B(K)
  460                CONTINUE
  480             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, JOB is not equal to ''P'' or ''U'': JOB',
     *       ' = ''',A,'''')
99998 FORMAT (1X,'** On entry, UPLO is not equal to ''U'' or ''L'': UP',
     *       'LO = ''',A,'''')
99997 FORMAT (1X,'** On entry, DIAG is not equal to ''N'', ''U'' or ''',
     *       'B'': DIAG = ''',A,'''')
99996 FORMAT (1X,'** On entry, N.lt.1: N = ',I16)
99995 FORMAT (1X,'** On entry, LDA.lt.N: LDA = ',I16,', N = ',I16)
      END
