      SUBROUTINE G13EBF(TRANSF,N,M,L,A,LDS,B,STQ,Q,LDQ,C,LDM,R,S,K,H,U,
     *                  TOL,IWK,WK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     The algorithm calculates a combined measurement and time update
C     of one iteration of the time invariant Kalman filter. This
C     requires the matrix pair (A,C) to be in lower observer Hessenberg.
C     This transformation is carried out if TRANSF = 'T'. The update is
C     given for the square root covariance filter.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13EBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, L, LDM, LDQ, LDS, M, N
      LOGICAL           STQ
      CHARACTER         TRANSF
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDS,N), B(LDS,L), C(LDM,N), H(LDM,M),
     *                  K(LDS,M), Q(LDQ,*), R(LDM,M), S(LDS,N),
     *                  U(LDS,N), WK((N+M)*(N+M+L))
      INTEGER           IWK(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  TOLER
      INTEGER           I, IERROR, IFAULT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DSCAL, F01QCF, F06QFF,
     *                  F06QHF, G13EBY, G13EBZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (N.LT.1) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC,FMT=99998) M
      ELSE IF (L.LT.1) THEN
         WRITE (P01REC,FMT=99997) L
      ELSE IF (LDS.LT.N) THEN
         WRITE (P01REC,FMT=99996) LDS, N
      ELSE IF (LDM.LT.M) THEN
         WRITE (P01REC,FMT=99995) LDM, M
      ELSE IF ( .NOT. STQ .AND. LDQ.LT.L) THEN
         WRITE (P01REC,FMT=99994) LDQ, L
      ELSE IF (STQ .AND. LDQ.LT.1) THEN
         WRITE (P01REC,FMT=99993) LDQ
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC,FMT=99992) TOL
      ELSE IF (TRANSF.NE.'T' .AND. TRANSF.NE.'t' .AND. TRANSF.NE.
     *         'H' .AND. TRANSF.NE.'h') THEN
         WRITE (P01REC,FMT=99991) TRANSF
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         TOLER = DBLE(M*M)*X02AJF()
         IF (TOLER.LT.TOL) TOLER = TOL
         IF (TRANSF.EQ.'T' .OR. TRANSF.EQ.'t') THEN
            CALL F06QHF('G',N,N,0.0D0,1.0D0,U,LDS)
            CALL G13EBY(N,M,A,LDS,C,LDM,U,LDS,.TRUE.,.FALSE.)
            DO 20 I = 1, N
               CALL DCOPY(N-I+1,S(I,I),1,WK,1)
               CALL DGEMV('N',N,N-I+1,1.0D0,U(1,I),LDS,WK,1,0.0D0,
     *                    WK(N+I),N)
   20       CONTINUE
            IFAULT = 0
            CALL F01QCF(N,N,WK(N+1),N,WK,IFAULT)
            DO 40 I = 1, N
               CALL DCOPY(I,WK(N+(I-1)*N+1),1,S(I,1),LDS)
   40       CONTINUE
            CALL DGEMM('N','N',N,L,N,1.0D0,U,LDS,B,LDS,0.0D0,WK,N)
            CALL F06QFF('G',N,L,WK,N,B,LDS)
         END IF
         CALL G13EBZ(N,L,M,S,LDS,A,LDS,B,LDS,Q,LDQ,C,LDM,R,LDM,K,LDS,H,
     *               LDM,IWK,WK,TOLER,STQ,IFAULT)
         IF (IFAULT.EQ.2) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99990)
         ELSE IF (H(1,1).LT.0.0D0) THEN
            DO 60 I = 1, M
               CALL DSCAL(M-I+1,-1.0D0,H(I,I),1)
   60       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N .lt. 1: N = ',I16)
99998 FORMAT (' ** On entry, M .lt. 1: M = ',I16)
99997 FORMAT (' ** On entry, L .lt. 1: L = ',I16)
99996 FORMAT (' ** On entry, LDS .lt. N: LDS = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, LDM .lt. M: LDM = ',I16,' M = ',I16)
99994 FORMAT (' ** On entry, LDQ .lt. L: LDQ = ',I16,' L = ',I16)
99993 FORMAT (' ** On entry, LDQ .lt. 1: LDQ = ',I16)
99992 FORMAT (' ** On entry, TOL .lt. 0.0: TOL = ',E13.5)
99991 FORMAT (' ** On entry, TRANSF is not valid: TRANSF = ',A1)
99990 FORMAT (' ** The matrix H is singular')
      END
