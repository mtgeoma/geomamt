      SUBROUTINE F08PEY(A,B,C,D,RT1R,RT1I,RT2R,RT2I,CS,SN)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLANV2(A,B,C,D,RT1R,RT1I,RT2R,RT2I,CS,SN)
C
C  Purpose
C  =======
C
C  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
C  matrix in standard form:
C
C       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
C       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
C
C  where either
C  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
C  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
C  conjugate eigenvalues.
C
C  Arguments
C  =========
C
C  A       (input/output) DOUBLE PRECISION
C  B       (input/output) DOUBLE PRECISION
C  C       (input/output) DOUBLE PRECISION
C  D       (input/output) DOUBLE PRECISION
C          On entry, the elements of the input matrix.
C          On exit, they are overwritten by the elements of the
C          standardised Schur form.
C
C  RT1R    (output) DOUBLE PRECISION
C  RT1I    (output) DOUBLE PRECISION
C  RT2R    (output) DOUBLE PRECISION
C  RT2I    (output) DOUBLE PRECISION
C          The real and imaginary parts of the eigenvalues. If the
C          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
C          eigenvalues are a complex conjugate pair, RT1I > 0.
C
C  CS      (output) DOUBLE PRECISION
C  SN      (output) DOUBLE PRECISION
C          Parameters of the rotation matrix.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, BB, CC, CS1, DD, P, SAB, SAC, SIGMA, SN1,
     *                  TAU, TEMP
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF
      EXTERNAL          F06BNF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
C
C     Initialize CS and SN
C
      CS = ONE
      SN = ZERO
C
      IF (C.EQ.ZERO) THEN
         GO TO 20
C
      ELSE IF (B.EQ.ZERO) THEN
C
C        Swap rows and columns
C
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 20
      ELSE IF (A.EQ.D .AND. SIGN(ONE,B).NE.SIGN(ONE,C)) THEN
         GO TO 20
      ELSE
C
C        Make diagonal elements equal
C
         TEMP = A - D
         P = HALF*TEMP
         SIGMA = B + C
         TAU = F06BNF(SIGMA,TEMP)
         CS1 = SQRT(HALF*(ONE+ABS(SIGMA)/TAU))
         SN1 = -(P/(TAU*CS1))*SIGN(ONE,SIGMA)
C
C        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
C                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
C
         AA = A*CS1 + B*SN1
         BB = -A*SN1 + B*CS1
         CC = C*CS1 + D*SN1
         DD = -C*SN1 + D*CS1
C
C        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
C                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
C
         A = AA*CS1 + CC*SN1
         B = BB*CS1 + DD*SN1
         C = -AA*SN1 + CC*CS1
         D = -BB*SN1 + DD*CS1
C
C        Accumulate transformation
C
         TEMP = CS*CS1 - SN*SN1
         SN = CS*SN1 + SN*CS1
         CS = TEMP
C
         TEMP = HALF*(A+D)
         A = TEMP
         D = TEMP
C
         IF (C.NE.ZERO) THEN
            IF (B.NE.ZERO) THEN
               IF (SIGN(ONE,B).EQ.SIGN(ONE,C)) THEN
C
C                 Real eigenvalues: reduce to upper triangular form
C
                  SAB = SQRT(ABS(B))
                  SAC = SQRT(ABS(C))
                  P = SIGN(SAB*SAC,C)
                  TAU = ONE/SQRT(ABS(B+C))
                  A = TEMP + P
                  D = TEMP - P
                  B = B - C
                  C = ZERO
                  CS1 = SAB*TAU
                  SN1 = SAC*TAU
                  TEMP = CS*CS1 - SN*SN1
                  SN = CS*SN1 + SN*CS1
                  CS = TEMP
               END IF
            ELSE
               B = -C
               C = ZERO
               TEMP = CS
               CS = -SN
               SN = TEMP
            END IF
         END IF
      END IF
C
   20 CONTINUE
C
C     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
C
      RT1R = A
      RT2R = D
      IF (C.EQ.ZERO) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT(ABS(B))*SQRT(ABS(C))
         RT2I = -RT1I
      END IF
      RETURN
C
C     End of F08PEY (DLANV2)
C
      END
