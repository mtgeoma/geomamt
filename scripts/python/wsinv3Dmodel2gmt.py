#!/usr/bin/env python
# coding=utf-8
import sys
import string
import numpy
import wsinv3Dmodel

if (__name__=="__main__"):
    if (len(sys.argv)==1):
        print "usage: %s parameter=value" % sys.argv[0]
        print "parameter could be:"
        print "model: must be a 3D wsinv3d model (mandatory)"
        print "scale: multiply block sizes (default is 1)"
        print "slice: select slice in format [x|y|z]distance (default is z0)"
        print "rho  : select rho scale between linear or log10 (default linear)"
        sys.exit(1)

    model=""
    scale=1.0
    Slice="z0"
    rho="linear"
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="model"):
                model=parameter[1]
            elif(parameter[0]=="scale"):
                scale=float(parameter[1])
            elif(parameter[0]=="slice"):
                Slice=parameter[1]
            elif(parameter[0]=="rho"):
                rho=parameter[1]
            else:
                print "error: \""+parameter[0]+"\" unknow parameter"
                sys.exit(1)
        else:
            print "error: all parameters must have an equal symbol [=]"
            sys.exit(1)

    if(len(model)==0):
        print "must define model parameter"
        sys.exit(1)

    if( not (rho == "linear" or rho == "log10") ):
        print "parameter rho must be linear or log10"
        sys.exit(1)

    model=wsinv3Dmodel.wsinv3Dmodel(model)

    model.orig*=scale
    model.x*=scale
    model.y*=scale
    model.z*=scale

    Nx=len(model.x)
    Ny=len(model.y)
    Nz=len(model.z)

    if (Slice[0]=='z'):
        pos=float(Slice[1:])
        z=model.orig[2]
        k=0
        while (z<=pos and k<Nz):
            z+=model.z[k]
            k+=1

        if (pos<z or k==Nz):
            k=k-1

        SM=model.m[k,:,:]
        if (rho=="log10"):
            SM=numpy.log10(SM)
        print "# rho(%s) min= %f; rho max= %f" % (rho,numpy.amin(SM),numpy.amax(SM))
        x=model.orig[0]
        for i in range(0,Nx):
            y=model.orig[1]
            for j in range(0,Ny):
                print "> -Z%f" % SM[j,i]
                print "%f %f" % (y,x)
                print "%f %f" % (y+model.y[j],x)
                print "%f %f" % (y+model.y[j],x+model.x[i])
                print "%f %f" % (y,x+model.x[i])
                y+=model.y[j]
            x+=model.x[i]
    elif (Slice[0]=='x'):
        pos=float(Slice[1:])
        x=model.orig[0]
        i=0
        while (x<=pos and i<Nx):
            x+=model.x[i]
            i+=1

        if (pos<x or i==Nx):
            i=i-1

        SM=model.m[:,:,i]
        if (rho=="log10"):
            SM=numpy.log10(SM)
        print "# rho(%s) min= %f; rho max= %f" % (rho,numpy.amin(SM),numpy.amax(SM))
        z=model.orig[2]
        for k in range(0,Nz):
            y=model.orig[1]
            for j in range(0,Ny):
                print "> -Z%f" % SM[k,j]
                print "%f %f" % (y,z)
                print "%f %f" % (y+model.y[j],z)
                print "%f %f" % (y+model.y[j],z+model.z[k])
                print "%f %f" % (y,z+model.z[k])
                y+=model.y[j]
            z+=model.z[k]
    elif (Slice[0]=='y'):
        pos=float(Slice[1:])
        y=model.orig[1]
        j=0
        while (y<=pos and j<Ny):
            y+=model.y[j]
            j+=1

        if (pos<y or j==Ny):
            j=j-1

        SM=model.m[:,j,:]
        if (rho=="log10"):
            SM=numpy.log10(SM)
        print "# rho(%s) min= %f; rho max= %f" % (rho,numpy.amin(SM),numpy.amax(SM))
        z=model.orig[2]
        for k in range(0,Nz):
            x=model.orig[0]
            for i in range(0,Nx):
                print "> -Z%f" % SM[k,i]
                print "%f %f" % (x,z)
                print "%f %f" % (x+model.x[i],z)
                print "%f %f" % (x+model.x[i],z+model.z[k])
                print "%f %f" % (x,z+model.z[k])
                x+=model.x[i]
            z+=model.z[k]
    else:
        print "option slice must start with 'x' or 'y' or 'z'."
        sys.exit(1)
