      SUBROUTINE E04NBW(MODE,N,NZ,NFREE,NQ,UNITQ,KX,V,ZY,WRK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-590 (MAR 1988).
C
C     ******************************************************************
C     E04NBW  transforms the vector  v  in various ways using the
C     matrix  Q = ( Z  Y )  defined by the input parameters.
C
C        MODE               result
C        ----               ------
C
C          1                v = Z v
C          2                v = Y v
C          3                v = Q v
C
C     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
C     on output, v  is a full n-vector.
C
C
C          4                v = Z'v
C          5                v = Y'v
C          6                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is ordered as  ( v(free)  v(fixed) ).
C
C          7                v = Y'v
C          8                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is as in modes 5 and 6 except that v(fixed) is not
C     set.
C
C     Modes  1, 4, 7 and 8  do not involve  v(fixed).
C     Original F66 version  April 1983.
C     Fortran 77 version written  9-February-1985.
C     Level 2 BLAS added 10-June-1986.
C     This version of E04NBW dated 10-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NFREE, NQ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  V(N), WRK(N), ZY(NQ,*)
      INTEGER           KX(N)
C     .. Local Scalars ..
      INTEGER           J, J1, J2, K, L, LENV, NFIXED
C     .. External Subroutines ..
      EXTERNAL          F06FBF, DCOPY, DGEMV
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      J1 = 1
      J2 = NFREE
      IF (MODE.EQ.1 .OR. MODE.EQ.4) J2 = NZ
      IF (MODE.EQ.2 .OR. MODE.EQ.5 .OR. MODE.EQ.7) J1 = NZ + 1
      LENV = J2 - J1 + 1
      IF (MODE.LE.3) THEN
C        ===============================================================
C        Mode = 1, 2  or  3.
C        ===============================================================
C
         IF (NFREE.GT.0) CALL F06FBF(NFREE,ZERO,WRK,1)
C
C        Copy  v(fixed)  into the end of  wrk.
C
         IF (MODE.GE.2 .AND. NFIXED.GT.0) CALL DCOPY(NFIXED,V(NFREE+1),
     *       1,WRK(NFREE+1),1)
C
C        Set  WRK  =  relevant part of  ZY * V.
C
         IF (LENV.GT.0) THEN
            IF (UNITQ) THEN
               CALL DCOPY(LENV,V(J1),1,WRK(J1),1)
            ELSE
               CALL DGEMV('N',NFREE,J2-J1+1,ONE,ZY(1,J1),NQ,V(J1),1,ONE,
     *                    WRK,1)
            END IF
         END IF
C
C        Expand  WRK  into  V  as a full n-vector.
C
         CALL F06FBF(N,ZERO,V,1)
         DO 20 K = 1, NFREE
            J = KX(K)
            V(J) = WRK(K)
   20    CONTINUE
C
C        Copy  WRK(fixed)  into the appropriate parts of  V.
C
         IF (MODE.GT.1) THEN
            DO 40 L = 1, NFIXED
               J = KX(NFREE+L)
               V(J) = WRK(NFREE+L)
   40       CONTINUE
         END IF
C
      ELSE
C        ===============================================================
C        Mode = 4, 5, 6, 7  or  8.
C        ===============================================================
C        Put the fixed components of  V  into the end of  WRK.
C
         IF (MODE.EQ.5 .OR. MODE.EQ.6) THEN
            DO 60 L = 1, NFIXED
               J = KX(NFREE+L)
               WRK(NFREE+L) = V(J)
   60       CONTINUE
         END IF
C
C        Put the free  components of  V  into the beginning of  WRK.
C
         IF (NFREE.GT.0) THEN
            DO 80 K = 1, NFREE
               J = KX(K)
               WRK(K) = V(J)
   80       CONTINUE
C
C           Set  V  =  relevant part of  ZY' * WRK.
C
            IF (LENV.GT.0) THEN
               IF (UNITQ) THEN
                  CALL DCOPY(LENV,WRK(J1),1,V(J1),1)
               ELSE
                  CALL DGEMV('T',NFREE,J2-J1+1,ONE,ZY(1,J1),NQ,WRK,1,
     *                       ZERO,V(J1),1)
               END IF
            END IF
         END IF
C
C        Copy the fixed components of  WRK  into the end of  V.
C
         IF (NFIXED.GT.0 .AND. (MODE.EQ.5 .OR. MODE.EQ.6))
     *       CALL DCOPY(NFIXED,WRK(NFREE+1),1,V(NFREE+1),1)
      END IF
C
      RETURN
C
C     End of  E04NBW. (CMQMUL)
C
      END
