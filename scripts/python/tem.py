#!/usr/bin/env python
# uso: ./tem.py captureTEM0700.txt [block=5 frequency=8 window=152.6u delay=26u current=15]

import sys

def le_parametros(argv):
    # le parametros de entrada
    # arquivo dos dados e
    # define o que se deve buscar
    chaves={"block":1,"frequency":1,"current":1,"window":1,"delay":1}
    parametros={}
    parametros["arquivo"]=""
    for i in range(1,len(argv)):
        if(argv[i].find("=")!=-1):
            parametro=argv[i].split("=")
            if(chaves.get(parametro[0],None)==None):
                print "erro: chave desconhecida"
            else:
                if(parametro[0]=="block"):
                    parametros[parametro[0]]=int(parametro[1])
                else:
                    parametros[parametro[0]]=parametro[1]
        elif (len(parametros["arquivo"])==0):
            parametros["arquivo"]=argv[i]
        else:
            print "erro: um e apenas um parametro (nome do arquivo) deve ser sem sinal de igual"
    return parametros

def tofloat(bloco_dados):
    for i in range(0,len(bloco_dados)):
        for j in range (0,3):
            if(bloco_dados[i][j][-1]=="M"):
                bloco_dados[i][j]=float(bloco_dados[i][j][0:-1])*1e6
            elif(bloco_dados[i][j][-1]=="K"):
                bloco_dados[i][j]=float(bloco_dados[i][j][0:-1])*1e3
            elif(bloco_dados[i][j][-1]=="m"):
                bloco_dados[i][j]=float(bloco_dados[i][j][0:-1])*1e-3
            elif(bloco_dados[i][j][-1]=="u"):
                bloco_dados[i][j]=float(bloco_dados[i][j][0:-1])*1e-6
            else:
                bloco_dados[i][j]=float(bloco_dados[i][j])
        
parametros=le_parametros(sys.argv)
arquivo=open(parametros["arquivo"],'r').read()
arquivo=arquivo.splitlines()
for line in range(0,len(arquivo)):
    if(arquivo[line].find("Wn")!=-1 and arquivo[line].find("Mag")!=-1 and arquivo[line].find("Rho")!=-1):
        #encontrou um bloco de dados
        # coleta dados do cabecalho
        cabecalho={}
        cabecalho["block"]=int(arquivo[line-5].split()[0])
        cabecalho["frequency"]=(arquivo[line-2].split()[0])
        cabecalho["current"]=(arquivo[line-2].split()[6])
        cabecalho["window"]=(arquivo[line-2].split()[7])
        cabecalho["delay"]=(arquivo[line-2].split()[8])
        # coleta dados do bloco
        bloco_dados=[]
        line+=1
        dados=arquivo[line].split()
        while(len(dados)==3):
            bloco_dados.append(dados)
            line+=1
            dados=arquivo[line].split()
        tofloat(bloco_dados)
        chaves=[]
        for i in range (1,len(parametros)):
            # busca por chaves ativas
            if(parametros.get("block",None)!=None):
                chaves.append("block")
            elif(parametros.get("frequency",None)!=None):
                chaves.append("frequency")
            elif(parametros.get("current",None)!=None):
                chaves.append("current")
            elif(parametros.get("window",None)!=None):
                chaves.append("window")
            elif(parametros.get("delay",None)!=None):
                chaves.append("delay")
        chaves_validas=0
        for i in range(0,len(chaves)):
            if(parametros[chaves[i]]==cabecalho[chaves[i]]):
                chaves_validas+=1
        if(len(chaves)==0 or len(chaves)==chaves_validas):
            print "> ",cabecalho
            for i in range(0,len(bloco_dados)):
                print bloco_dados[i][0],abs(bloco_dados[i][1])
