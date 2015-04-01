      SUBROUTINE G02DKZ(IP,ICONST,P,IMP,C,LDC,B,S,SE,COV,WK,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES LEAST SQUARES PARAMETER ESTIMATES FOR GIVEN SET OF
C     LINEAR CONSTRAINTS GIVEN THE SVD SOLUTION FROM G02DAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S
      INTEGER           ICONST, IMP, IND, IP, LDC
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), C(LDC,ICONST), COV((IP*IP+IP)/2), P(*),
     *                  SE(IP), WK(2*IP*IP+IP*ICONST+2*ICONST*ICONST+4*
     *                  ICONST)
C     .. Local Scalars ..
      DOUBLE PRECISION  VARB
      INTEGER           I, IC2, IJ, IP2, IPC1, IRANK, LDP
C     .. External Subroutines ..
      EXTERNAL          F04AEF, DCOPY, DGEMV, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IRANK = IP - ICONST
      IPC1 = (IP+1)*ICONST
      IC2 = ICONST*ICONST
      IF (IMP.EQ.0) THEN
         IP2 = 2*IP + 1
         LDP = IP
      ELSE
         IP2 = 1
         LDP = IMP
      END IF
C
C       FIND C'P0 AND C'B
C
      DO 20 I = 1, ICONST
         CALL DGEMV('T',IP,ICONST,1.0D0,C,LDC,P(IP2+IRANK+I-1),LDP,
     *              0.0D0,WK(IPC1+1+(I-1)*ICONST),1)
   20 CONTINUE
      CALL DGEMV('T',IP,ICONST,1.0D0,C,LDC,B,1,0.0D0,WK(IPC1+1+IC2),1)
C
C        FIND INV(C'P0)*B AND INV(C'P0)*C'
C
      DO 40 I = 1, ICONST
         CALL DCOPY(IP,C(1,I),1,WK(IPC1+IC2+ICONST+I),ICONST)
   40 CONTINUE
      IND = 1
      CALL F04AEF(WK(IPC1+1),ICONST,WK(IPC1+IC2+1),ICONST,ICONST,IP+1,
     *            WK,ICONST,WK(2*IPC1+IC2+1),WK(2*IPC1+IC2+ICONST+1),
     *            ICONST,WK(2*(IPC1+IC2)+ICONST+1),ICONST,IND)
      IF (IND.EQ.0) THEN
C
C        FIND NEW B AND COV(NEW B) = D*P1*P1'*D'*S
C                      WHERE D=I-P0*INV(C'P0)*C'
C
         CALL DGEMV('T',ICONST,IP,-1.0D0,P(IP2+IRANK),LDP,WK,1,1.0D0,B,
     *              1)
         DO 60 I = 1, IP
            CALL DGEMV('T',ICONST,IP,-1.0D0,P(IP2+IRANK),LDP,
     *                 WK(I*ICONST+1),1,0.0D0,WK(IPC1+(I-1)*IP+1),1)
            WK(IPC1+(I-1)*IP+I) = WK((I-1)*IP+IPC1+I) + 1.0D0
   60    CONTINUE
         DO 80 I = 1, IP
            CALL DGEMV('N',IRANK,IP,1.0D0,P(IP2),LDP,WK(IPC1+I),IP,
     *                 0.0D0,WK(IPC1+IP*IP+(I-1)*IRANK+1),1)
   80    CONTINUE
         IJ = 1
         DO 100 I = 1, IP
            CALL DCOPY(IRANK,WK(IPC1+IP*IP+(I-1)*IRANK+1),1,WK,1)
            CALL DGEMV('T',IRANK,I,1.0D0,WK(IPC1+IP*IP+1),IRANK,WK,1,
     *                 0.0D0,COV(IJ),1)
            IJ = IJ + I
  100    CONTINUE
         CALL DSCAL((IP*(IP+1)/2),S,COV,1)
         DO 120 I = 1, IP
            VARB = COV((I*I+I)/2)
            IF (VARB.GT.0.0D0) THEN
               SE(I) = SQRT(VARB)
            ELSE
               SE(I) = 0.0D0
            END IF
  120    CONTINUE
      END IF
      END
