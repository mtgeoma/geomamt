#!/usr/bin/env python
# coding=utf-8

import string

def readJformat(arq_jformat):
    componentes={}
    arq=open(arq_jformat,'r').read()
    arq=arq.splitlines()
    i=0
    # pula comentários
    while(len(arq[i])>0 and arq[i][0]=="#"):
        i+=1
    # lê cabeçalho
    dados={}
    while(len(arq[i])>0 and arq[i][0]==">"):
        col=arq[i].split()
        if(col[0]==">STATION"):
            dados["STATION"]=col[1][1:]
        elif(col[0]==">AZIMUTH"):
            dados["AZIMUTH"]=col[2]
        elif(col[0]==">LATITUDE"):
            dados["LATITUDE"]=col[2]
        elif(col[0]==">LONGITUDE"):
            dados["LONGITUDE"]=col[2]
        elif(col[0]==">ELEVATION"):
            dados["ELEVATION"]=col[2]
        else:
            print "%s: rótulo desconhecido no cabeçalho" % col[0]
            sys.exit()
        i+=1
    componentes["HEAD"]=dados
    # lê blocos de dados de cada componente
    while(i<(len(arq)-1)):
        i+=1
        col=arq[i].split()
        rotulo=col[0]
        i+=1
        N=int(arq[i])
        dados=[]
        for n in range(0,N):
            i+=1
            col=arq[i].split()
            dados.append([float(s) for s in col])
            if(dados[n][0]<0.):
                dados[n][0]=-1./dados[n][0]

        componentes[rotulo]=dados
    return componentes
