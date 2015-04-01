#!/usr/bin/env python
# coding=utf-8
# le os arquivos x e y (assume que aparecem nesta ordem) dados em arq com o tamanho das grades em metros
#
# ATENÇÃO! Inverte o sentido e origem do eixo Y (Y'=Ymax-Y) mas não muda a numeração das células e a ordem em que são apresentadas!
# usos:
# 1) para obter a coordenada do centro da celula ao longo da coluna 3 do modelo:
# grade.py x=3 arq=x.dat,y.dat
# 2) para obter a coordenada do centro da celula ao longo da linha 5 do modelo:
# grade.py y=5 arq=x.dat,y.dat
# 3) para obter a coordenada do centro da celula da coluna 3 e linha 5 modelo:
# grade.py x=3 y=5 arq=x.dat,y.dat
# 4) para obter a coordenada do centro de todas as celulas do modelo:
# grade.py arq=x.dat,y.dat
# 5) para obter grade do modelo, ja separada com ">" no começo de cada looping interno
# grade.py saida=grade arq=x.dat,y.dat
# para as opções 4 e 5, o looping interno é do y, para ser o contrario, use "LI=x"
# para obter a grade rotacionada 30 graus em torno do seu centro, use "azimute=30"

import sys
import string
import math

def rotaciona(X,Y,Xcentro,Ycentro,azimute):
    if(azimute==0.0):
        return (X,Y)
    else:
        X-=Xcentro
        Y-=Ycentro
        x=X*math.cos(azimute)+Y*math.sin(azimute)
        y=X*(-1.0)*math.sin(azimute)+Y*math.cos(azimute)
        x+=Xcentro
        y+=Ycentro
        return (x,y)

if (__name__=="__main__"):
    x=-1
    y=-1
    LI="y"
    saida=""
    azimute=0.0
    # le parametros de entrada
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parametro=sys.argv[i].split("=")
            if(parametro[0]=="x"):
                x=int(parametro[1])
            elif(parametro[0]=="y"):
                y=int(parametro[1])
            elif(parametro[0]=="arq"):
                arq=parametro[1].split(",")
            elif(parametro[0]=="saida"):
                saida=parametro[1]
            elif(parametro[0]=="LI"):
                LI=parametro[1]
            elif(parametro[0]=="azimute"):
                azimute=float(parametro[1])*(math.pi/180.0)
            else:
                print "erro: \""+parametro[0]+"\" parâmetro desconhecido"
        else:
            print "erro: todos os parâmetros devem ter um sinal de igual"
    # lê as larguras das células
    largura_x=open(arq[0],'r').read().splitlines()
    largura_y=open(arq[1],'r').read().splitlines()
    X=[]
    Y=[]
    X.append(0.0)
    for largura in largura_x:
        X.append(float(largura.split(None,1)[0])+X[-1])
    Xcentro=X[-1]/2.0
    Y.append(0.0)
    for largura in largura_y:
        Y.append(float(largura.split(None,1)[0])+Y[-1])
    Ycentro=Y[-1]/2.0
    # inverte o sentido e origem do eixo Y (Y'=Ymax-Y)
    for i in range(0,len(Y)):
        Y[i]=Y[-1]-Y[i]

    if(saida=="grade"): # caso 5
        if(LI=="y"):
            for j in range (0,len(X)):
                print ">"
                for i in range (0,len(Y)):
                    p=rotaciona(X[j],Y[i],Xcentro,Ycentro,azimute)
                    print "%f %f %d %d" % (p[0],p[1],j,i)
            for j in range (0,len(Y)):
                print ">"
                for i in range (0,len(X)):
                    p=rotaciona(X[i],Y[j],Xcentro,Ycentro,azimute)
                    print "%f %f %d %d" % (p[0],p[1],i,j)
        else:
            for j in range (0,len(Y)):
                print ">"
                for i in range (0,len(X)):
                    p=rotaciona(X[i],Y[j],Xcentro,Ycentro,azimute)
                    print "%f %f %d %d" % (p[0],p[1],i,j)
            for j in range (0,len(X)):
                print ">"
                for i in range (0,len(Y)):
                    p=rotaciona(X[j],Y[i],Xcentro,Ycentro,azimute)
                    print "%f %f %d %d" % (p[0],p[1],j,i)
    else:
        if(x==-1 and y==-1): # caso 4
            if(LI=="y"):
                for j in range (1,len(X)):
                    for i in range (1,len(Y)):
                        p=rotaciona((X[j]+X[j-1])/2.0,(Y[i]+Y[i-1])/2.0,Xcentro,Ycentro,azimute)
                        print "%f %f %d %d" % (p[0],p[1],j,i)
            else:
                for j in range (1,len(Y)):
                    for i in range (1,len(X)):
                        p=rotaciona((X[i]+X[i-1])/2.0,(Y[j]+Y[j-1])/2.0,Xcentro,Ycentro,azimute)
                        print "%f %f %d %d" % (p[0],p[1],i,j)
        elif(y==-1): # caso 1
            for i in range (1,len(Y)):
                p=rotaciona((X[x]+X[x-1])/2.0,(Y[i]+Y[i-1])/2.0,Xcentro,Ycentro,azimute)
                print "%f %f %d %d" % (p[0],p[1],x,i)
        elif(x==-1): # caso 2
            for i in range (1,len(X)):
                p=rotaciona((X[i]+X[i-1])/2.0,(Y[y]+Y[y-1])/2.0,Xcentro,Ycentro,azimute)
                print "%f %f %d %d" % (p[0],p[1],i,y)
        else: # caso 3
                p=rotaciona((X[x]+X[x-1])/2.0,(Y[y]+Y[y-1])/2.0,Xcentro,Ycentro,azimute)
                print "%f %f %d %d" % (p[0],p[1],x,y)
