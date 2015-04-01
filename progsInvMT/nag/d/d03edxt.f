      SUBROUTINE D03EDX(U,V,LEV,NGU,NGV,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     U on GRID(LEV)=prolongation of V on GRID(LEV-1), LEV>1
C
C     .. Scalar Arguments ..
      INTEGER           LEV, NGU, NGV
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NGU), V(NGV)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           J, K, NPC, NPF, NPFO, NX, NX2
C     .. Executable Statements ..
      NX2 = 2*NGRIDX(LEV)
      NPC = NGP(LEV) - NGRIDX(LEV-1)
      NPFO = NGP(LEV) - NGRIDX(LEV)*NGRIDY(LEV) - NX2
      NPF = NPFO
      DO 40 K = 1, NGRIDY(LEV-1)
         NPC = NPC + NGRIDX(LEV-1)
         NPF = NPF + NX2
         DO 20 J = 1, NGRIDX(LEV-1)
            U(NPF+2*J-1) = V(NPC+J)
   20    CONTINUE
   40 CONTINUE
      NPF = NPFO
      DO 80 K = 1, NGRIDY(LEV), 2
         NPF = NPF + NX2
         DO 60 J = 2, NGRIDX(LEV), 2
            U(NPF+J) = (U(NPF+J-1)+U(NPF+J+1))/2.0D0
   60    CONTINUE
   80 CONTINUE
      NX = NGRIDX(LEV)
      NPF = NPFO + NX
      DO 120 K = 2, NGRIDY(LEV), 2
         NPF = NPF + NX2
         DO 100 J = 1, NX, 2
            U(NPF+J) = (U(NPF+J-NX)+U(NPF+J+NX))/2.0D0
  100    CONTINUE
  120 CONTINUE
      NPF = NPFO + NX
      DO 160 K = 2, NGRIDY(LEV), 2
         NPF = NPF + NX2
         DO 140 J = 2, NX, 2
            U(NPF+J) = (U(NPF+J+NX-1)+U(NPF+J-NX+1))/2.0D0
  140    CONTINUE
  160 CONTINUE
      RETURN
      END
