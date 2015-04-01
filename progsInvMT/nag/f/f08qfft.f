      SUBROUTINE F08QFF(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DTREXC(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,WORK,INFO)
C
C  Purpose
C  =======
C
C  DTREXC reorders the real Schur factorization of a real matrix
C  A = Q*T*Q', so that the diagonal block of T with row index IFST is
C  moved to row ILST.
C
C  The real Schur form T is reordered by an orthogonal similarity
C  transformation Z'*T*Z, and optionally the matrix Q of Schur vectors
C  is updated by postmultiplying it with Z.
C
C  T must be in Schur canonical form, that is, block upper triangular
C  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
C  has its diagonal elemnts equal and its off-diagonal elements of
C  opposite sign.
C
C  Arguments
C  =========
C
C  COMPQ   (input) CHARACTER*1
C          = 'V': update the matrix Q of Schur vectors;
C          = 'N': do not update Q.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C          On entry, the upper quasi-triangular matrix T, in Schur
C          Schur canonical form.
C          On exit, T is overwritten by the reordered matrix T, again in
C          Schur canonical form..
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          orthogonal transformation matrix Z which reorders T.
C          If COMPQ = 'N', Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= 1; and if COMPQ = 'V', LDQ >= max(1,N).
C
C  IFST    (input/output) INTEGER
C  ILST    (input/output) INTEGER
C          Specify the re-ordering of the diagonal blocks of T.
C          The block with row index IFST is moved to row ILST, by a
C          sequence of transpositions between adjacent blocks.
C          On exit, if IFST pointed on entry to the second row of a
C          2-by-2 block, it is changed to point to the first row; ILST
C          always points to the first row of the block in its final
C          position (which may differ from its input value by +1 or -1).
C          1 <= IFST <= N; 1 <= ILST <= N.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C          = 1: two adjacent blocks were too close to swap (the problem
C               is very ill-conditioned); T may have been partially
C               reordered, and ILST points to the first row of the
C               current position of the block being moved.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFST, ILST, INFO, LDQ, LDT, N
      CHARACTER         COMPQ
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*), T(LDT,*), WORK(*)
C     .. Local Scalars ..
      INTEGER           HERE, II, J, NBF, NBL, NBNEXT
      LOGICAL           WANTQ
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DSCAL, F06AAZ, F08QFZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Decode and test the input arguments.
C
      INFO = 0
      WANTQ = (COMPQ.EQ.'V' .OR. COMPQ.EQ.'v')
      IF ( .NOT. WANTQ .AND. .NOT. (COMPQ.EQ.'N' .OR. COMPQ.EQ.'n'))
     *    THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDT.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LDQ.LT.1 .OR. (WANTQ .AND. LDQ.LT.MAX(1,N))) THEN
         INFO = -6
      ELSE IF (IFST.LT.1 .OR. IFST.GT.N) THEN
         INFO = -7
      ELSE IF (ILST.LT.1 .OR. ILST.GT.N) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QFF/DTREXC',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.1) RETURN
C
C     Determine the first row of specified block
C     and find out it is 1 by 1 or 2 by 2.
C
      IF (IFST.GT.1) THEN
         IF (T(IFST,IFST-1).NE.ZERO) IFST = IFST - 1
      END IF
      NBF = 1
      IF (IFST.LT.N) THEN
         IF (T(IFST+1,IFST).NE.ZERO) NBF = 2
      END IF
C
C     Determine the first row of the final block
C     and find out it is 1 by 1 or 2 by 2.
C
      IF (ILST.GT.1) THEN
         IF (T(ILST,ILST-1).NE.ZERO) ILST = ILST - 1
      END IF
      NBL = 1
      IF (ILST.LT.N) THEN
         IF (T(ILST+1,ILST).NE.ZERO) NBL = 2
      END IF
C
      IF (IFST.EQ.ILST) RETURN
C
      IF (IFST.LT.ILST) THEN
C
C        Update ILST
C
         IF (NBF.EQ.2 .AND. NBL.EQ.1) ILST = ILST - 1
         IF (NBF.EQ.1 .AND. NBL.EQ.2) ILST = ILST + 1
C
         HERE = IFST
C
   20    CONTINUE
C
C        Swap block with next one below
C
         IF (NBF.EQ.1 .OR. NBF.EQ.2) THEN
C
C           Current block either 1 by 1 or 2 by 2
C
            NBNEXT = 1
            IF (HERE+NBF+1.LE.N) THEN
               IF (T(HERE+NBF+1,HERE+NBF).NE.ZERO) NBNEXT = 2
            END IF
            CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,NBF,NBNEXT,WORK,INFO)
            IF (INFO.NE.0) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
C
C           Test if 2 by 2 block breaks into two 1 by 1 blocks
C
            IF (NBF.EQ.2) THEN
               IF (T(HERE+1,HERE).EQ.ZERO) NBF = 3
            END IF
C
         ELSE
C
C           Current block consists of two 1 by 1 blocks each of which
C           must be swapped individually
C
            NBNEXT = 1
            IF (HERE+3.LE.N) THEN
               IF (T(HERE+3,HERE+2).NE.ZERO) NBNEXT = 2
            END IF
            CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE+1,1,NBNEXT,WORK,INFO)
            IF (INFO.NE.0) THEN
               ILST = HERE
               RETURN
            END IF
            IF (NBNEXT.EQ.1) THEN
C
C              Swap two 1 by 1 blocks, no problems possible
C
               CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,1,NBNEXT,WORK,INFO)
               HERE = HERE + 1
            ELSE
C
C              Recompute NBNEXT in case 2 by 2 split
C
               IF (T(HERE+2,HERE+1).EQ.ZERO) NBNEXT = 1
               IF (NBNEXT.EQ.2) THEN
C
C                 2 by 2 Block did not split
C
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,1,NBNEXT,WORK,
     *                        INFO)
                  IF (INFO.NE.0) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
C
C                 2 by 2 Block did split
C
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,1,1,WORK,INFO)
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE+1,1,1,WORK,INFO)
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF (HERE.LT.ILST) GO TO 20
C
      ELSE
C
         HERE = IFST
   40    CONTINUE
C
C        Swap block with next one above
C
         IF (NBF.EQ.1 .OR. NBF.EQ.2) THEN
C
C           Current block either 1 by 1 or 2 by 2
C
            NBNEXT = 1
            IF (HERE.GE.3) THEN
               IF (T(HERE-1,HERE-2).NE.ZERO) NBNEXT = 2
            END IF
            CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE-NBNEXT,NBNEXT,NBF,WORK,
     *                  INFO)
            IF (INFO.NE.0) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
C
C           Test if 2 by 2 block breaks into two 1 by 1 blocks
C
            IF (NBF.EQ.2) THEN
               IF (T(HERE+1,HERE).EQ.ZERO) NBF = 3
            END IF
C
         ELSE
C
C           Current block consists of two 1 by 1 blocks each of which
C           must be swapped individually
C
            NBNEXT = 1
            IF (HERE.GE.3) THEN
               IF (T(HERE-1,HERE-2).NE.ZERO) NBNEXT = 2
            END IF
            CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE-NBNEXT,NBNEXT,1,WORK,
     *                  INFO)
            IF (INFO.NE.0) THEN
               ILST = HERE
               RETURN
            END IF
            IF (NBNEXT.EQ.1) THEN
C
C              Swap two 1 by 1 blocks, no problems possible
C
               CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,NBNEXT,1,WORK,INFO)
               HERE = HERE - 1
            ELSE
C
C              Recompute NBNEXT in case 2 by 2 split
C
               IF (T(HERE,HERE-1).EQ.ZERO) NBNEXT = 1
               IF (NBNEXT.EQ.2) THEN
C
C                 2 by 2 Block did not split
C
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE-1,2,1,WORK,INFO)
                  IF (INFO.NE.0) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
C
C                 2 by 2 Block did split
C
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE,1,1,WORK,INFO)
                  CALL F08QFZ(WANTQ,N,T,LDT,Q,LDQ,HERE-1,1,1,WORK,INFO)
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF (HERE.GT.ILST) GO TO 40
      END IF
      ILST = HERE
C
      IF (WANTQ) THEN
C
C        Normalize Schur vectors so that element of largest absolute
C        value is positive.
C
         DO 60 J = 1, N
            II = IDAMAX(N,Q(1,J),1)
            IF (Q(II,J).LT.ZERO) THEN
               CALL DSCAL(N,-ONE,Q(1,J),1)
               IF (J.GT.1) THEN
                  CALL DSCAL(J-1,-ONE,T(1,J),1)
                  T(J,J-1) = -T(J,J-1)
               END IF
               IF (J.LT.N) THEN
                  CALL DSCAL(N-J,-ONE,T(J,J+1),LDT)
                  T(J+1,J) = -T(J+1,J)
               END IF
            END IF
   60    CONTINUE
      END IF
      RETURN
C
C     End of F08QFF (DTREXC)
C
      END
