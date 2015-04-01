#!/usr/bin/env python
# coding=utf-8
import sys
import string
import wsinv3Dmodel

if (__name__=="__main__"):
    if (len(sys.argv)==1):
        print "usage: %s parameter=value" % sys.argv[0]
        print "parameter could be:"
        print "model: must be a 3D wsinv3d model (mandatory)"
        sys.exit(1)

    model=""
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="model"):
                model=parameter[1]
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

    print "+-----------------------------------------------------------------------------+"
    print "| This file defines model covariance for a recursive autoregression scheme.   |"
    print "| The model space may be divided into distinct areas using integer masks.     |"
    print "| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |"
    print "| air, ocean and the rest of the model is turned off automatically. You can   |"
    print "| also define exceptions to override smoothing between any two model areas.   |"
    print "| To turn off smoothing set it to zero. This header is 16 lines long.         |"
    print "| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |"
    print "| 2. Smoothing in the X direction (NzEarth real values)                       |"
    print "| 3. Smoothing in the Y direction (NzEarth real values)                       |"
    print "| 4. Vertical smoothing (1 real value)                                        |"
    print "| 5. Number of times the smoothing should be applied (1 integer >= 0)         |"
    print "| 6. Number of exceptions (1 integer >= 0)                                    |"
    print "| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |"
    print "| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|"
    print "+-----------------------------------------------------------------------------+"

    print "\n%d %d %d" % (len(model.x),len(model.y),len(model.z))
    line=""
    for k in range(0,len(model.z)):
        line+=" 0.3"
    print "\n%s" % line
    print "\n%s" % line
    print "\n0.3"
    print "\n1\n\n0\n"

    print "%d %d" % (1,len(model.z))
    line=""
    for j in range(0,len(model.y)):
        line+=" 1"
    for i in range(0,len(model.x)):
        print line
