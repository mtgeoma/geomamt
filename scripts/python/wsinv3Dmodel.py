#!/usr/bin/env python
# coding=utf-8

import sys
import traceback
import string
import numpy

class wsinv3Dmodel:
    "read and write model wsinv3DMT"
    def __read__(self,N,ln,model):
        a=[]
        while (len(a)!=N):
            ln+=1
            if(ln>=len(model)):
                break
            col=model[ln].split()
            for c in col:
                a.append(c)
        return ln,a
        
    def __init__(self,model):
        # variables to keep model data
        self.comment_line=""
        self.rotation=0
        self.orig=numpy.array([])
        self.x=numpy.array([])
        self.y=numpy.array([])
        self.z=numpy.array([])
        self.m=numpy.array([])

        # read wsinv3DMT from file model
        try:
            model=open(model,'r').read()
        except:
            (ErrorType,ErrorValue,ErrorTB)=sys.exc_info()
            traceback.print_exc(ErrorTB)
            sys.exit(1)
        model=model.splitlines()

        ln=0 # line number
        self.comment_line=model[ln]
        ln+=1
        col=model[ln].split()
        Nx=int(col[0])
        Ny=int(col[1])
        Nz=int(col[2])
        Nr=int(col[3])
        rho="scalar"
        if(len(col)>=5 and col[4]=="LOGE"):
            rho="loge"

        ln,self.x=self.__read__(Nx,ln,model)
        if(len(self.x)!=Nx):
            print "couldn't find all X in model file"
            sys.exit(1)
        self.x=numpy.array([float(s) for s in self.x])

        ln,self.y=self.__read__(Ny,ln,model)
        if(len(self.y)!=Ny):
            print "couldn't find all Y in model file"
            sys.exit(1)
        self.y=numpy.array([float(s) for s in self.y])

        ln,self.z=self.__read__(Nz,ln,model)
        if(len(self.z)!=Nz):
            print "couldn't find all Z in model file"
            sys.exit(1)
        self.z=numpy.array([float(s) for s in self.z])

        if ( Nr == 0 ):
            ln,self.m=self.__read__(Nx*Ny*Nz,ln,model)
            if(len(self.m)!=Nx*Ny*Nz):
                print "couldn't find all matrix elements in model file"
                sys.exit(1)
            self.m=numpy.array([float(s) for s in self.m])
            if (rho=="loge"):
                self.m=numpy.exp(self.m)
            self.m.shape = (Nz,Ny,Nx)
            # invert axes X order
            self.m=self.m[:,:,::-1]
        else:
            print "for now, only works with Nr==0"
            sys.exit(1)

        ln,self.orig=self.__read__(3,ln,model)
        if(len(self.orig)==0):
            self.orig=("%f %f %f" % (-numpy.sum(self.x)/2.,
                                      -numpy.sum(self.y)/2.,0.)).split()
        elif(len(self.orig)!=3):
            print "couldn't find all elements of model origin"
            sys.exit(1)
        self.orig=numpy.array([float(s) for s in self.orig])

        a=[]
        ln,a=self.__read__(2,ln,model)
        if(len(a)==0):
            self.rotation=float(0)
        elif(len(a)>1):
            print "find unexpected elements at the end of model"
            sys.exit(1)
        else:
            self.rotation=float(a[0])

    def write(self):
        print self.comment_line
        print " %4d %4d %4d %4d LOGE" % (len(self.x),len(self.y),len(self.z),0)
        self.m=numpy.log(self.m)
        # invert axes X order
        self.m=self.m[:,:,::-1]

        out=""
        for d in self.x:
            out+=" %11.3f" % d
        print out
        out=""
        for d in self.y:
            out+=" %11.3f" % d
        print out
        out=""
        for d in self.z:
            out+=" %11.3f" % d
        print out

        for k in range(0,len(self.z)):
            print ""
            for j in range(0,len(self.y)):
                out=""
                for i in range(0,len(self.x)):
                    out+=" %12.5E" % self.m[k,j,i]
                print out
        out=""
        for d in self.orig:
            out+=" %15.3f" % d
        print ""
        print out
        print "%9.3f" % self.rotation
