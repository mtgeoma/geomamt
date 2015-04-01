      SUBROUTINE D03EDW(U,V,LEV,NGU,NGV,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     U on GRID(LEV)=restriction of V on GRID(LEV+1)
C
C     .. Scalar Arguments ..
      INTEGER           LEV, NGU, NGV
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NGU), V(NGV)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           J, JEND, K, KEND, NF, NPC, NPCC, NPF, NPFF, NXC,
     *                  NXF
C     .. Executable Statements ..
      NXC = NGRIDX(LEV)
      NXF = NGRIDX(LEV+1)
      NPCC = NGP(LEV) - NXC*NGRIDY(LEV)
      NPFF = NGP(LEV+1) - NXF*NGRIDY(LEV+1)
      NPC = NPCC
      NPF = NPFF - 1
      JEND = NGRIDY(LEV) - 1
      KEND = NXC - 1
      DO 40 J = 2, JEND
         NPC = NPC + NXC
         NPF = NPF + 2*NXF
         DO 20 K = 2, KEND
            NF = NPF + 2*K
            U(NPC+K) = V(NF) + (V(NF+1)+V(NF-1)+V(NF+NXF)+V(NF-NXF)
     *                 +V(NF+NXF-1)+V(NF-NXF+1))/2.0D0
   20    CONTINUE
   40 CONTINUE
      NPC = NPCC
      NPF = NPFF + 1
      DO 60 K = 2, KEND
         NPF = NPF + 2
         U(NPC+K) = V(NPF) + (V(NPF+1)+V(NPF-1)+V(NPF+NXF)+V(NPF+NXF-1))
     *              /2.0D0
   60 CONTINUE
      NPC = NPCC + NXC*(NGRIDY(LEV)-1)
      NPF = NPFF + (NGRIDY(LEV+1)-1)*NXF + 1
      DO 80 K = 2, KEND
         NPF = NPF + 2
         U(NPC+K) = V(NPF) + (V(NPF+1)+V(NPF-1)+V(NPF-NXF)+V(NPF-NXF+1))
     *              /2.0D0
   80 CONTINUE
      NPC = NPCC + 1
      NPF = NPFF + 1
      DO 100 J = 2, JEND
         NPF = NPF + 2*NXF
         NPC = NPC + NXC
         U(NPC) = V(NPF) + (V(NPF+1)+V(NPF+NXF)+V(NPF-NXF)+V(NPF-NXF+1))
     *            /2.0D0
  100 CONTINUE
      NPC = NPCC + NXC
      NPF = NPFF + NXF
      DO 120 J = 2, JEND
         NPC = NPC + NXC
         NPF = NPF + 2*NXF
         U(NPC) = V(NPF) + (V(NPF-1)+V(NPF+NXF)+V(NPF-NXF)+V(NPF+NXF-1))
     *            /2.0D0
  120 CONTINUE
      NPC = NPCC + 1
      NPF = NPFF + 1
      U(NPC) = V(NPF) + (V(NPF+NXF)+V(NPF+1))/2.0D0
      NPC = NPCC + NXC
      NPF = NPFF + NXF
      U(NPC) = V(NPF) + (V(NPF-1)+V(NPF+NXF)+V(NPF+NXF-1))/2.0D0
      NPC = NGP(LEV) - NXC + 1
      NPF = NGP(LEV+1) - NXF + 1
      U(NPC) = V(NPF) + (V(NPF+1)+V(NPF-NXF)+V(NPF-NXF+1))/2.0D0
      NPC = NGP(LEV)
      NPF = NGP(LEV+1)
      U(NPC) = V(NPF) + (V(NPF-1)+V(NPF-NXF))/2.0D0
      RETURN
      END
