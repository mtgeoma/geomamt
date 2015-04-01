      SUBROUTINE G13AJY(MR,NP,ND,NQ,NPS,NDS,NQS,NS,NPD,NDD,NQD,MPQS,
     *                  NPAR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     .. Scalar Arguments ..
      INTEGER           ND, NDD, NDS, NP, NPAR, NPD, NPS, NQ, NQD, NQS,
     *                  NS
C     .. Array Arguments ..
      INTEGER           MPQS(4), MR(7)
C     .. Executable Statements ..
      NP = MR(1)
      ND = MR(2)
      NQ = MR(3)
      NPS = MR(4)
      NDS = MR(5)
      NQS = MR(6)
      NS = MR(7)
      NPD = NP + NPS*NS
      NDD = ND + NDS*NS
      NQD = NQ + NQS*NS
      MPQS(1) = NP
      MPQS(2) = NQ
      MPQS(3) = NPS
      MPQS(4) = NQS
      NPAR = NP + NQ + NPS + NQS
      RETURN
      END
