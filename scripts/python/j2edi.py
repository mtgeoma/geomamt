#!/usr/bin/env python
# coding=utf-8
import sys
import string
from datetime import date
import jformat

def deg2hms(deg):
    deg=deg.split(".")
    deg[1]="0."+deg[1]
    frac=float(deg[1])
    min=int(frac*60)
    sec=frac*60*60-min*60
    saida="%s:%02d:%02.0f" % (deg[0],min,sec)
    return saida

if (__name__=="__main__"):
    arq_jformat=sys.argv[1]
    Z=jformat.readJformat(arq_jformat)
    componentes=("ZXX","ZXY","ZYX","ZYY")
    # verifica se as componentes tem os mesmos períodos
    for comp in componentes:
        for j in range(0,len(Z[comp])):
            if(Z[comp][j][0]!=Z["ZXX"][j][0]):
                print "componentes ZXX e %s com %d-esimo período diferentes:" % (comp,j+1)
                print "%s!=%s" % (Z["ZXX"][j][0],Z[comp][j][0])
                sys.exit()
    # saída do HEAD
    print ">HEAD"
    print "  DATAID=\"%s\"" % Z["HEAD"]["STATION"]
    print "  ACQBY=geoma"
    print "  HLEBY=geoma"
    print "  ACQDATE=%s" % date.today().strftime("%m/%d/%y")
    print "  FILEDATE=%s" % date.today().strftime("%m/%d/%y")
    print "  LAT=%s" % deg2hms(Z["HEAD"]["LATITUDE"])
    print "  LONG=%s" % deg2hms(Z["HEAD"]["LONGITUDE"])
    print "  ELEV=%s" % Z["HEAD"]["ELEVATION"]
    print "  UNITS=M"

    print "\n>INFO"

    print "\n>=DEFINEMEAS"
    print "  UNITS=M"
    print "  REFTYPE=CART"
    print "  REFLOC=%s" % Z["HEAD"]["STATION"]
    print "  REFLAT=%s" % deg2hms(Z["HEAD"]["LATITUDE"])
    print "  REFLONG=%s" % deg2hms(Z["HEAD"]["LONGITUDE"])
    print "  REFELEV=%s" % Z["HEAD"]["ELEVATION"]
    print ">! only CHTYPE and AZM are relevant and correct!"
    print ">HMEAS ID=1 CHTYPE=HX X=0.0 Y=0.0 AZM=%.1f" % (0.0)
    print ">HMEAS ID=2 CHTYPE=HY X=0.0 Y=0.0 AZM=%.1f" % (90.0)
    print ">EMEAS ID=3 CHTYPE=EX X=0.0 Y=0.0 X2=0.0 Y2=0.0"
    print ">EMEAS ID=4 CHTYPE=EY X=0.0 Y=0.0 X2=0.0 Y2=0.0"

    print "\n>=MTSECT"
    print "  NFREQ=%d" % len(Z[componentes[0]])
    print "\n>FREQ  //%d" % len(Z[componentes[0]])
    saida=""
    for j in range(0,len(Z[componentes[0]])):
        T=Z[componentes[0]][j][0]
        saida+=" %16.9E" % (1.0/T)
        if((j+1)%5==0):
            print saida
            saida=""
    if(len(Z[componentes[0]])%5!=0):
        print saida

    print "\n>ZROT  //%d" % len(Z[componentes[0]])
    saida=""
    for j in range(0,len(Z[componentes[0]])):
        saida+=" %16.9E" % (float(Z["HEAD"]["AZIMUTH"]))
        if((j+1)%5==0):
            print saida
            saida=""
    if(len(Z[componentes[0]])%5!=0):
        print saida

    for comp in componentes:
        print "\n>%sR  //%d" % (comp,len(Z[comp]))
        saida=""
        for j in range(0,len(Z[comp])):
            R=Z[comp][j][1]
            saida+=" %16.9E" % (R)
            if((j+1)%5==0):
                print saida
                saida=""
        if(len(Z[comp])%5!=0):
            print saida

        print "\n>%sI  //%d" % (comp,len(Z[comp]))
        saida=""
        for j in range(0,len(Z[comp])):
            I=Z[comp][j][2]
            saida+=" %16.9E" % (I)
            if((j+1)%5==0):
                print saida
                saida=""
        if(len(Z[comp])%5!=0):
            print saida

        print "\n>%s.VAR  //%d" % (comp,len(Z[comp]))
        saida=""
        for j in range(0,len(Z[comp])):
            STD=Z[comp][j][3]
            saida+=" %16.9E" % (STD**2)
            if((j+1)%5==0):
                print saida
                saida=""
        if(len(Z[comp])%5!=0):
            print saida
    print ">END"
