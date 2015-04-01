      SUBROUTINE D02LAZ(A,C,BPHAT,B,BP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     THIS ROUTINE RETURNS THE COEFFICIENTS FOR THE
C     DORMAND/PRINCE RKN6(4)6FD PAIR USED IN THE CODE D02LAF
C     TO COMPUTE Y AND YP.  THE COEFFICIENTS FOR THE DENSE
C     OUTPUT FORMULA ARE RETURNED BY THE ROUTINE D02LZZ TO
C     THE INTERPOLATION SUBROUTINE.
C
C     .. Array Arguments ..
      DOUBLE PRECISION  A(17,17), B(16), BP(16), BPHAT(17), C(17)
C     .. Executable Statements ..
      C(1) = 0.0D0
      C(2) = 1.29295903136704415288990053209D-1
      C(3) = 2.58591806273408830577980106418D-1
      C(4) = 6.70297082615480058310908782471D-1
      C(5) = 9.0D-1
      C(6) = 1.0D0
C
      A(2,1) = 8.35871528396802532822102346743D-3
      A(3,1) = 1.11449537119573671042946979566D-2
      A(3,2) = 2.22899074239147342085893959132D-2
      A(4,1) = 1.45474742801091785895935232316D-1
      A(4,2) = -2.29860640522647473120261456297D-1
      A(4,3) = 3.09034987202967536528726080729D-1
      A(5,1) = -2.07668262950789954335146205252D-1
      A(5,2) = 6.86366784292514312273571851042D-1
      A(5,3) = -1.99549277872349252201326358888D-1
      A(5,4) = 1.25850756530624894262900713098D-1
      A(6,1) = 7.81101614434947768281101614435D-2
      A(6,2) = 0.0D0
      A(6,3) = 2.882917411897667776841772093D-1
      A(6,4) = 1.22425537174570410182422422005D-1
      A(6,5) = 1.1172560192168035305290207251D-2
C
      BPHAT(1) = 7.81101614434947768281101614435D-2
      BPHAT(2) = 0.0D0
      BPHAT(3) = 3.88843478705982602715952208523D-1
      BPHAT(4) = 3.71320757928842267403035557523D-1
      BPHAT(5) = 1.1172560192168035305290207251D-1
      BPHAT(6) = 5.0D-2
C
      B(1) = 1.05885926037041827822001005886D0
      B(2) = -2.40675137192445205319176777632D0
      B(3) = 1.84789211155403377497175771746D0
      B(4) = 0.0D0
      B(5) = 0.0D0
      B(6) = 0.0D0
C
      BP(1) = 5.46058879392212725546058879392D-2
      BP(2) = 0.0D0
      BP(3) = 4.6126678590362684429468353488D-1
      BP(4) = 1.95880859479312656291875875208D-1
      BP(5) = 3.88246466677839226858834701972D-1
      BP(6) = -1.0D-1
C
      RETURN
      END