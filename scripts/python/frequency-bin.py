#!/usr/bin/env python
# coding=utf-8

# calculate frequency bin (n)

import sys
import string
import math

if (__name__=="__main__"):
    n=0
    N=0
    sp=0.0
    sf=0.0
    T=0.0
    f=0.0
    if(len(sys.argv)==1):
        print "usage: %s parameter=value" % sys.argv[0]
        print "parameter could be:"
        print "sp=sampling period or sf=sampling frequency"
        print "T=period or f=frequency"
        print "n=frequency bin"
        print "N=window size"
        print "define 3 independents parameters and the 4rd will be calculated"
        print "couldn't define T and f together (are dependents)"
        print "couldn't define sp and sf together (are dependents)"
        print "Example:"
        print "to find the freq. bin equivalent to 60Hz when using a window"
        print "of 8192 points and data colected at a sampling freq. of 4096Hz:"
        print "%s f=60 sf=4096 N=8192" % sys.argv[0]
        sys.exit()
    # read input parameters
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="sp"):
                sp=float(parameter[1])
            elif(parameter[0]=="sf"):
                sf=float(parameter[1])
            elif(parameter[0]=="T"):
                T=float(parameter[1])
            elif(parameter[0]=="f"):
                f=float(parameter[1])
            elif(parameter[0]=="n"):
                n=int(parameter[1])
            elif(parameter[0]=="N"):
                N=int(parameter[1])
            else:
                print "error: \""+parameter[0]+"\" unknow parameter"
                sys.exit()
        else:
            print "error: all parameters must have an equal symbol [=]"
            sys.exit()
    # only f or T must be defined
    if(f>0.0 and T > 0.0):
        print "error: only f or T must be defined"
        sys.exit()
    if(f>0.0):
        T=1.0/f
    elif(T>0.0):
        f=1.0/T
    # only sf or sp must be defined
    if(sf>0.0 and sp > 0.0):
        print "error: only sf or sp must be defined"
        sys.exit()
    elif(sf>0.0):
        sp=1.0/sf
    else:
        sf=1.0/sp

    # f = n*sf/N; T=1/f; sp=1/sf
    # defined 3 independent parameters, the 4rd will be calculated
    if(n>0 and N > 0 and sf > 0.0):
        f=float(n)*sf/float(N)
        print "f=%f T=%f" % (f,1.0/f)
    elif(n>0 and N > 0 and f > 0.0):
        sf=float(N)*f/float(n)
        print "sf=%f sp=%f" % (sf,1.0/sf)
    elif(n>0 and sf > 0.0 and T > 0.0):
        N=int(float(n)*sf*T+0.5)
        print "N=%d" % N
    elif(N>0 and sp > 0.0 and f > 0.0):
        n=int(float(N)*f*sp+0.5)
        print "n=%d" % n
    else:
        print "unexpected combination"
