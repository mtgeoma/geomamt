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
        print "slice: select slice in format [x|y|z]distance (default is all model)"
        sys.exit(1)

    model=""
    scale=1.0
    Slice="all"
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="model"):
                model=parameter[1]
            elif(parameter[0]=="scale"):
                scale=float(parameter[1])
            elif(parameter[0]=="slice"):
                Slice=parameter[1]
            else:
                print "error: \""+parameter[0]+"\" unknow parameter"
                sys.exit(1)
        else:
            print "error: all parameters must have an equal symbol [=]"
            sys.exit(1)

    if(len(model)==0):
        print "must define model parameter"
        sys.exit(1)

    model=wsinv3Dmodel.wsinv3Dmodel(model)
    if(len(model.m)==0):
        print "fail to get model"
        sys.exit(1)

    model.orig*=scale
    model.x*=scale
    model.y*=scale
    model.z*=scale

    print "rotation: %f" % model.rotation
    print "Yorig/Xorig: %f/%f" % (model.orig[1],model.orig[0])
    print "Ymin/Ymax/Xmin/Xmax: %f/%f/%f/%f" % (
     0+model.orig[1],numpy.sum(model.y)+model.orig[1],
     0+model.orig[0],numpy.sum(model.x)+model.orig[0])

    if (Slice == "all"):
        print "# rho_min= %f rho_max= %f" % (min(model.m.flat),max(model.m.flat))
    else:
        if (Slice[0]=='z'):
            pos=float(Slice[1:])
            z=0.0
            k=0
            while (z<=pos and k<len(model.z)):
                z+=model.z[k]
                k+=1

            if (pos<z or k==len(model.z)):
                k=k-1

            SM=model.m[k,:,:]
            print "# rho_min= %f rho_max= %f" % (min(SM.flat),max(SM.flat))
        elif (Slice[0]=='x'):
            pos=float(Slice[1:])
            x=0.0
            i=0
            while (x<=pos and i<len(model.x)):
                x+=model.x[i]
                i+=1

            if (pos<x or i==len(model.x)):
                i=i-1

            SM=model.m[:,i,:]
            print "# rho_min= %f rho_max= %f" % (min(SM.flat),max(SM.flat))
        elif (Slice[0]=='y'):
            pos=float(Slice[1:])
            y=0.0
            j=0
            while (y<=pos and j<len(model.y)):
                y+=model.y[j]
                j+=1

            if (pos<y or j==len(model.y)):
                j=j-1

            SM=model.m[:,:,j]
            print "# rho_min= %f rho_max= %f" % (min(SM.flat),max(SM.flat))
        else:
            print "option slice must start with 'x' or 'y' or 'z'."
            sys.exit(1)
