#!/usr/bin/env python
# coding=utf-8
import sys
import string
import numpy

if (__name__=="__main__"):
    if (len(sys.argv)==1):
        print "usage: %s parameter=value" % sys.argv[0]
        print "parameter could be:"
        print "model: must be a 2D mackie model (mandatory)"
        print "scale: multiply block sizes (default is 1)"
        sys.exit()

    model=""
    scale=1.0
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="model"):
                model=parameter[1]
            elif(parameter[0]=="scale"):
                scale=float(parameter[1])
            else:
                print "error: \""+parameter[0]+"\" unknow parameter"
                sys.exit()
        else:
            print "error: all parameters must have an equal symbol [=]"
            sys.exit()

    if(len(model)==0):
        print "must define model parameter"
        sys.exit()

    model=open(model,'r').read()
    model=model.splitlines()

    ln=0 # line number
    col=model[ln].split()
    Ny=int(col[0])
    Nz=int(col[1])
    rho="scalar"
    if(len(col)==3 and col[2]=="LOGE"):
        rho="loge"

    YB=[]
    while (len(YB)!=Ny):
        ln+=1
        col=model[ln].split()
        for c in col:
            YB.append(float(c)*scale)

    ZB=[]
    while (len(ZB)!=Nz):
        ln+=1
        col=model[ln].split()
        for c in col:
            ZB.append(float(c)*scale)

    ln+=1
    col=model[ln].split()
    mode=int(col[0])
    if (not(mode == 0 or mode ==1)):
        print "mode must be 0 or 1"
        sys.exit(1)

    M=[]
    while (len(M)!=Ny*Nz):
        ln+=1
        col=model[ln].split()
        for c in col:
            if (mode==0):
                if (rho=="loge"):
                    M.append(numpy.exp(float(c)))
                else:
                    M.append(float(c))
            else:
                M.append(int(c))

    if (mode==1):
        ln+=1
        col=model[ln].split()
        for i in range(0,len(M)):
            M[i]=float(col[M[i]])

    M=numpy.array(M)

    print "# rho min= %f; rho max= %f" % (numpy.amin(M),numpy.amax(M))
    M.shape = (Nz,Ny)
    z=0
    for i in range(0,Nz):
        y=0
        for j in range(0,Ny):
            print "> -Z%f" % M[i,j]
            print "%f %f" % (y,z)
            print "%f %f" % (y+YB[j],z)
            print "%f %f" % (y+YB[j],z+ZB[i])
            print "%f %f" % (y,z+ZB[i])
            y+=YB[j]
        z+=ZB[i]
