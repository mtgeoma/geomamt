c
c  The following routines are from the book Numerical Recipes
c
      FUNCTION BETACF(A,B,X)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.
        IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
11    CONTINUE
      PAUSE 'A or B too big, or ITMAX too small'
1     BETACF=AZ
      RETURN
      END

      FUNCTION BETAI(A,B,X)
      IF(X.LT.0..OR.X.GT.1.)PAUSE 'bad argument X in BETAI'
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *      +A*ALOG(X)+B*ALOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
      END

      SUBROUTINE CRANK(N,W,S)
      DIMENSION W(N)
      S=0.
      J=1
1     IF(J.LT.N)THEN
        IF(W(J+1).NE.W(J))THEN
          W(J)=J
          J=J+1
        ELSE
          DO 11 JT=J+1,N
            IF(W(JT).NE.W(J))GO TO 2
11        CONTINUE
          JT=N+1
2         RANK=0.5*(J+JT-1)
          DO 12 JI=J,JT-1
            W(JI)=RANK
12        CONTINUE
          T=JT-J
          S=S+T**3-T
          J=JT
        ENDIF
      GO TO 1
      ENDIF
      IF(J.EQ.N)W(N)=N
      RETURN
      END

      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END

      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

      SUBROUTINE SORT2(N,RA,RB)
      DIMENSION RA(N),RB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END

      SUBROUTINE SPEAR(DATA1,DATA2,N,WKSP1,WKSP2,D,ZD,PROBD,RS,PROBRS)
      DIMENSION DATA1(N),DATA2(N),WKSP1(N),WKSP2(N)
      DO 11 J=1,N
        WKSP1(J)=DATA1(J)
        WKSP2(J)=DATA2(J)
11    CONTINUE
      CALL SORT2(N,WKSP1,WKSP2)
      CALL CRANK(N,WKSP1,SF)
      CALL SORT2(N,WKSP2,WKSP1)
      CALL CRANK(N,WKSP2,SG)
      D=0.
      DO 12 J=1,N
        D=D+(WKSP1(J)-WKSP2(J))**2
12    CONTINUE
      EN=N
      EN3N=EN**3-EN
      AVED=EN3N/6.-(SF+SG)/12.
      FAC=(1.-SF/EN3N)*(1.-SG/EN3N)
      VARD=((EN-1.)*EN**2*(EN+1.)**2/36.)*FAC
      ZD=(D-AVED)/SQRT(VARD)
      PROBD=ERFCC(ABS(ZD)/1.4142136)
      RS=(1.-(6./EN3N)*(D+0.5*(SF+SG)))/FAC
      T=RS*SQRT((EN-2.)/((1.+RS)*(1.-RS)))
      DF=EN-2.
      PROBRS=BETAI(0.5*DF,0.5,DF/(DF+T**2))
      RETURN
      END


      FUNCTION RTBIS(FUNC,X1,X2,XACC,P,NU)
      PARAMETER (JMAX=40)
      FMID=FUNC(X2,P,NU)
      F=FUNC(X1,P,NU)
      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID,P,NU)
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      PAUSE 'too many bisections'
      END

      FUNCTION GAMMQ(A,X)
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END

      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

