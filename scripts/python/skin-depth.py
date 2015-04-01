#!/usr/bin/env python
# coding=utf-8

# para calcular o skin-depth (sd)
# ./skin-depth.py rho=100.0 T=3.94784176044E-4
# para calcular o rho
# ./skin-depth.py  sd=100.0 T=3.94784176044E-4
# para calcular o período
# ./skin-depth.py  sd=100.0 rho=100.0

import sys
import string
import math

if (__name__=="__main__"):
    # sd = C*sqrt(rho*T); onde C=1000.0*sqrt(10.0)/(2*PI); no SI
    C=1.0E3*math.sqrt(10.0)/(2.0*math.pi)
    rho=0.0
    sd=0.0
    T=0.0
    f=0.0
    # le parametros de entrada
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parametro=sys.argv[i].split("=")
            if(parametro[0]=="rho"):
                rho=float(parametro[1])
            elif(parametro[0]=="T"):
                T=float(parametro[1])
                f=1./T
            elif(parametro[0]=="f"):
                f=float(parametro[1])
                T=1./f
            elif(parametro[0]=="sd"):
                sd=float(parametro[1])
            else:
                print "erro: \""+parametro[0]+"\" parâmetro desconhecido"
        else:
            print "erro: todos os parâmetros devem ter um sinal de igual"
    # encontrando dois parâmetros maiores que zero, assume que se quer calcular o terceiro parâmetro
    if(rho>0.0 and T > 0.0):
        print C*math.sqrt(rho*T)
    elif(rho>0.0 and sd > 0.0):
        print (sd/C)**2/rho
    else:
        print (sd/C)**2/T
