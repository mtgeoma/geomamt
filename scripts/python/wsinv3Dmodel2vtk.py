#!/usr/bin/env python
# coding=utf-8
import sys
import string
import numpy as np
import wsinv3Dmodel
from evtk.hl import gridToVTK

if (__name__=="__main__"):
    if (len(sys.argv)==1):
        print "usage: %s parameter=value" % sys.argv[0]
        print "parameter could be:"
        print "model: must be a 3D wsinv3d model (mandatory)"
        print "scale: multiply block sizes (default is 1)"
        sys.exit(1)

    model_file=""
    scale=1.0
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parameter=sys.argv[i].split("=")
            if(parameter[0]=="model"):
                model_file=parameter[1]
            elif(parameter[0]=="scale"):
                scale=float(parameter[1])
            else:
                print "error: \""+parameter[0]+"\" unknow parameter"
                sys.exit(1)
        else:
            print "error: all parameters must have an equal symbol [=]"
            sys.exit(1)

    if(len(model_file)==0):
        print "must define model parameter"
        sys.exit(1)

    model=wsinv3Dmodel.wsinv3Dmodel(model_file)
    if(len(model.m)==0):
        print "fail to get model"
        sys.exit(1)

    Nx=len(model.x)
    Ny=len(model.y)
    Nz=len(model.z)

    Xgrid=(np.add.accumulate(np.append([0.],model.x))+model.orig[0])*scale
    Ygrid=(np.add.accumulate(np.append([0.],model.y))+model.orig[1])*scale
    Zgrid=(np.add.accumulate(np.append([0.],model.z))+model.orig[2])*scale

    # in wsinv3Dmodel.py: self.m.shape = (Nz,Ny,Nx)
    gridToVTK(model_file, Zgrid, Ygrid, Xgrid,
              cellData = {"resistivity" : model.m.copy()})

    # print "# vtk DataFile Version 3.0"
    # print "%s" % model.comment_line
    # print "ASCII"
    # print "DATASET RECTILINEAR_GRID"
    # print "DIMENSIONS %d %d %d" % (Nx,Ny,Nz)
    # print "X_COORDINATES %d float" % Nx
    # print "%s" % print_grid(model.x,model.orig[0])
    # print "Y_COORDINATES %d float" % Ny
    # print "%s" % print_grid(model.y,model.orig[1])
    # print "Z_COORDINATES %d float" % Nz
    # print "%s" % print_grid(model.z,model.orig[2])
    # print "CELL_DATA 1"
    # print "POINT_DATA %d" % (Nx*Ny*Nz)
    # print "SCALARS resistivity float"
    # print "LOOKUP_TABLE default"
