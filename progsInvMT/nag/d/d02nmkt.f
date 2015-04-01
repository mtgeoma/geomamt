      SUBROUTINE D02NMK(N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME JGROUP
C
C-----------------------------------------------------------------------
C THIS SUBROUTINE CONSTRUCTS GROUPINGS OF THE COLUMN INDICES OF
C THE JACOBIAN MATRIX, USED IN THE NUMERICAL EVALUATION OF THE
C JACOBIAN BY FINITE DIFFERENCES.
C
C INPUT..
C N      = THE ORDER OF THE MATRIX.
C IA,JA  = SPARSE STRUCTURE DESCRIPTORS OF THE MATRIX BY ROWS.
C MAXG   = LENGTH OF AVAILABLE STORATE IN THE IGP ARRAY.
C
C OUTPUT..
C NGRP   = NUMBER OF GROUPS.
C JGP    = ARRAY OF LENGTH N CONTAINING THE COLUMN INDICES BY GROUPS.
C IGP    = POINTER ARRAY OF LENGTH NGRP + 1 TO THE LOCATIONS IN JGP
C          OF THE BEGINNING OF EACH GROUP.
C IER    = ERROR INDICATOR.  IER = 0 IF NO ERROR OCCURRED, OR 1 IF
C          MAXG WAS INSUFFICIENT.
C
C INCL AND JDONE ARE WORKING ARRAYS OF LENGTH N.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IER, MAXG, N, NGRP
C     .. Array Arguments ..
      INTEGER           IA(*), IGP(*), INCL(N), JA(*), JDONE(N), JGP(N)
C     .. Local Scalars ..
      INTEGER           I, J, K, KMAX, KMIN, NCOL, NG
C     .. Executable Statements ..
      IER = 0
      DO 20 J = 1, N
         JDONE(J) = 0
   20 CONTINUE
      NCOL = 1
      DO 120 NG = 1, MAXG
         IGP(NG) = NCOL
         DO 40 I = 1, N
            INCL(I) = 0
   40    CONTINUE
         DO 100 J = 1, N
C REJECT COLUMN J IF IT IS ALREADY IN A GROUP.--------------------------
            IF (JDONE(J).EQ.1) GO TO 100
            KMIN = IA(J)
            KMAX = IA(J+1) - 1
            DO 60 K = KMIN, KMAX
C REJECT COLUMN J IF IT OVERLAPS ANY COLUMN ALREADY IN THIS GROUP.------
               I = JA(K)
               IF (INCL(I).EQ.1) GO TO 100
   60       CONTINUE
C ACCEPT COLUMN J INTO GROUP NG.----------------------------------------
            JGP(NCOL) = J
            NCOL = NCOL + 1
            JDONE(J) = 1
            DO 80 K = KMIN, KMAX
               I = JA(K)
               INCL(I) = 1
   80       CONTINUE
  100    CONTINUE
C STOP IF THIS GROUP IS EMPTY (GROUPING IS COMPLETE).-------------------
         IF (NCOL.EQ.IGP(NG)) GO TO 140
  120 CONTINUE
C ERROR RETURN IF NOT ALL COLUMNS WERE CHOSEN (MAXG TOO SMALL).---------
      IF (NCOL.LE.N) GO TO 160
      NG = MAXG
  140 NGRP = NG - 1
      RETURN
  160 IER = 1
      RETURN
      END
