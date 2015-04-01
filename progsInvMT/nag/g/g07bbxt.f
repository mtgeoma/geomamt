      SUBROUTINE G07BBX(XMU,XSIG,N,T,S,L1,L2,L11,L12,L22,WK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     EVALUATES FIRST AND SECOND DERIVATIVES OF THE LOG LIKELIHOOD
C     FOR CENSORED NORMAL DATA.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  L1, L11, L12, L2, L22, S, T, XMU, XSIG
C     .. Array Arguments ..
      DOUBLE PRECISION  WK(*)
      INTEGER           N(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, D, E, F, G, LL, P, PE, RR, XVAR, Z, ZZ
      INTEGER           IFAULT, J, N1, N2, N3
C     .. External Functions ..
      DOUBLE PRECISION  G01MAZ, S15ABF
      EXTERNAL          G01MAZ, S15ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      N1 = N(1)
      N2 = N1 + N(2)
      N3 = N2 + N(3) + N(3)
      IFAULT = 0
      XVAR = XSIG**2
      PE = DBLE(N(4))
      A = 0.0D0
      B = 0.0D0
      C = 0.0D0
      D = 0.0D0
      E = 0.0D0
      F = 0.0D0
      G = 0.0D0
      DO 20 J = 1, N1
         RR = (WK(J)-XMU)/XSIG
         P = 1.0D0 - S15ABF(RR,IFAULT)
         IF (P.GT.0.0D0) THEN
            Z = G01MAZ(RR)/(XSIG*P)
            A = A + Z
            E = E + Z*Z
            F = F + Z*WK(J)*Z
            Z = WK(J)*Z
            B = B + Z
            G = G + Z**2
            Z = WK(J)*Z
            C = C + Z
            Z = WK(J)*Z
            D = D + Z
         END IF
   20 CONTINUE
      DO 40 J = N1 + 1, N2
         LL = (WK(J)-XMU)/XSIG
         P = S15ABF(LL,IFAULT)
         IF (P.GT.0.0D0) THEN
            Z = G01MAZ(LL)/(XSIG*P)
            A = A - Z
            E = E + Z*Z
            F = F + Z*WK(J)*Z
            Z = WK(J)*Z
            B = B - Z
            G = G + Z**2
            Z = WK(J)*Z
            C = C - Z
            Z = WK(J)*Z
            D = D - Z
         END IF
   40 CONTINUE
      DO 60 J = N2 + 1, N3, 2
         LL = (WK(J+1)-XMU)/XSIG
         RR = (WK(J)-XMU)/XSIG
         P = S15ABF(LL,IFAULT) - S15ABF(RR,IFAULT)
         IF (P.GT.0.0D0) THEN
            Z = G01MAZ(RR)/XSIG
            ZZ = G01MAZ(LL)/XSIG
            A = A + (Z-ZZ)/P
            E = E + ((Z-ZZ)/P)*((Z-ZZ)/P)
            F = F + ((Z-ZZ)/P)*((WK(J)*Z-WK(J+1)*ZZ)/P)
            Z = WK(J)*Z
            ZZ = WK(J+1)*ZZ
            B = B + (Z-ZZ)/P
            G = G + ((Z-ZZ)/P)**2
            Z = WK(J)*Z
            ZZ = WK(J+1)*ZZ
            C = C + (Z-ZZ)/P
            Z = WK(J)*Z
            ZZ = WK(J+1)*ZZ
            D = D + (Z-ZZ)/P
         END IF
   60 CONTINUE
C
      L1 = A - (PE*XMU-T)/XVAR
      L2 = (B-XMU*A-PE)/XSIG + (S+XMU*(PE*XMU-2.0D0*T))/(XVAR*XSIG)
      L11 = (B-XMU*A-PE)/XVAR - E
      L12 = (XMU*E-F-A)/XSIG + (C-2.0D0*T+XMU*(A*XMU-2.0D0*B+2.0D0*PE))
     *      /(XVAR*XSIG)
      Z = PE - 2.0D0*B - G + XMU*(2.0D0*F-XMU*E+2.0D0*A)
      ZZ = D - 3.0D0*S + XMU*(6.0D0*T-3.0D0*C+XMU*
     *     (3.0D0*B-XMU*A-3.0D0*PE))
      Z = Z + ZZ/XVAR
      L22 = Z/XVAR
      RETURN
      END
