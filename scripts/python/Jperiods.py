#!/usr/bin/env python
# coding=utf-8
import sys
import string
import math
import jformat

if (__name__=="__main__"):
    if (len(sys.argv)==1):
        print "usage: %s list.txt [dlogT]" % sys.argv[0]
        sys.exit(1)

    list=sys.argv[1]
    list=open(list,'r').read()
    list=list.splitlines()

    dlogT=0.
    if (len(sys.argv)==3):
        dlogT=float(sys.argv[2])

    Tstr = {}
    Tfloat = {}
    for line in list:
        Jfile=line.split(None,1)[0]
        Z=jformat.readJformat(Jfile)
        components=["ZXX","ZXY","ZYX","ZYY"]
        for comp in components:
            for j in range(0,len(Z[comp])):
                T="%.4e" % Z[comp][j][0]
                if (not T in Tstr):
                    Tstr[T]=0
                    Tfloat[Z[comp][j][0]]=T
                Tstr[T]+=1

    T0=0
    Tlist = {}
    for T in sorted(Tfloat):
        if(T0==0):
            # start a new Tlist
            Tlist = {}
            Tlist[Tstr[Tfloat[T]]] = [Tfloat[T]]
            T0=T
        else:
            if (math.log(T/T0)<=dlogT):
                if (Tstr[Tfloat[T]] in Tlist):
                    Tlist[Tstr[Tfloat[T]]].append(Tfloat[T])
                else:
                    Tlist[Tstr[Tfloat[T]]] = [Tfloat[T]]
            else:
                # print previous Tlist
                if (len(Tlist)>0):
                    P=""
                    C=""
                    for n in reversed(sorted(Tlist)):
                        for t in Tlist[n]:
                            P="%s %s" % (P,t)
                            C="%s %d" % (C,n)

                    print "%s #%s" % (P[1:],C)
                # start a new Tlist
                Tlist = {}
                Tlist[Tstr[Tfloat[T]]] = [Tfloat[T]]
                T0=T

    # print previous Tlist
    if (len(Tlist)>0):
        P=""
        C=""
        for n in reversed(sorted(Tlist)):
            for t in Tlist[n]:
                P="%s %s" % (P,t)
                C="%s %d" % (C,n)

        print "%s #%s" % (P[1:],C)
