#!/usr/bin/env python
# coding=utf-8
# ponto=x,y marca o par ordenado. Se uma ordenada for marcada com o valor zero, assume-se que esta ordenada é indefinida
# com x e y definidos, saída será no formato Jones neste ponto
# le_d3imp.py p=1,1
# com x, y e T definidos, saída será o vetor de indução neste ponto
# le_d3imp.py p=1,1 NT=1
# com x e T definidos, saída será vetores de indução nesta linha e período
# le_d3imp.py p=1,0 NT=1
# com y e T definidos, saída será vetores de indução nesta linha e período
# le_d3imp.py p=0,1 NT=1
#
# pode-se entrar com um arquivo na variável p, as duas primeiras colunas desse arquivo serão tratados como pontos
# onde
# x e/ou y marca o ponto/linha que se quer recuperar
# NT numero do periodo que se quer analizar
# use e=2 para multiplicar por 2 o modulo do vetor de indução
# ATENÇÃO! Inverte o sentido e origem do eixo Y (Y'=Ymax-Y) mas não muda a numeração das células e a ordem em que são apresentadas!

import sys
import string
import numpy
import math

def le_bloco(arquivo_imp,line,nX,nY,Re,Im,f):
    # A primeira linha contém o rótulo do bloco, deve-se pular
    line+=1
    f.append(float(arquivo_imp[line]))
    for cy in range(0,nY):
        line+=1
        cx=0
        dados=arquivo_imp[line].split()
        NF=len(dados)
        while(NF!=0):
            if(NF%2!=0):
                print "ERRO! Linha %d com número de campos impar!" % (line)
            else:
                for c in range(0,NF):
                    if(c%2==1):
                        Im[cx,cy]=float(dados[c])
                        cx+=1
                    else:
                        Re[cx,cy]=float(dados[c])
                line+=1
                dados=arquivo_imp[line].split()
                NF=len(dados)
        if(cx!=nX):
            print "ERRO! cx!=nX"
    line+=1
    # retorna a linha com o rótulo da próxima componente
    return line

def saida_jones(saida,rotulo,T,Re,Im,i,j):
    saida.write(rotulo)
    line="%d\n" % len(T)
    saida.write(line)
    for t in range(0,len(T)):
        real=Re[t][i,j]
        imag=-1*Im[t][i,j]
        erro=0.05*math.sqrt(real**2+imag**2)
        line= "%14.4e%14.4e%14.4e%14.4e  1\n" % (T[t],real,imag,erro)
        saida.write(line)

if (__name__=="__main__"):
    NT=0
    e=0
    arq_imp=""
    arq_pontos=""
    # le parametros de entrada
    for i in range(1,len(sys.argv)):
        if(sys.argv[i].find("=")!=-1):
            parametro=sys.argv[i].split("=")
            if(parametro[0]=="p"):
                arq_pontos=parametro[1]
            elif(parametro[0]=="NT"):
                NT=int(parametro[1])
            elif(parametro[0]=="e"):
                e=float(parametro[1])
            else:
                print "erro: \""+parametro[0]+"\" parâmetro desconhecido"
        elif(len(arq_imp)==0):
            arq_imp=sys.argv[i]
        else:
            print "erro: um e apenas um parâmetro (arquivo *.imp) deve ser sem sinal de igual"
    # fator de multiplicacao do vetor
    if(e==0):
        e=1
    # cria uma lista de pontos
    pontos=[]
    if(arq_pontos.find(",")!=-1):
        ponto=arq_pontos.split(",",2)
        pontos.append([int(ponto[0]),int(ponto[1])])
    else:
        arquivo_pontos=open(arq_pontos,'r').read()
        arquivo_pontos=arquivo_pontos.splitlines()
        for line in arquivo_pontos:
            ponto=line.split(None,2)
            pontos.append([int(ponto[0]),int(ponto[1])])

    # le os dados do arquivo *.imp
    arquivo_imp=open(arq_imp,'r').read()
    arquivo_imp=arquivo_imp.splitlines()
    line=0
    # guarda as posicoes de X
    while(arquivo_imp[line].find("X-positions")==-1):
        line+=1
    line+=1
    X=[]
    while(arquivo_imp[line].find("Y-positions")==-1):
        posicoes=arquivo_imp[line].split()
        for i in range (0,len(posicoes)):
            X.append(float(posicoes[i]))
        line+=1
    nX=len(X)
    line+=1
    # guarda as posicoes de Y
    Y=[]
    while(arquivo_imp[line].find("ZXX")==-1):
        posicoes=arquivo_imp[line].split()
        for i in range (0,len(posicoes)):
            Y.append(float(posicoes[i]))
        line+=1
    nY=len(Y)
    T=[]
    RXX=[]
    IXX=[]
    RXY=[]
    IXY=[]
    RYX=[]
    IYX=[]
    RYY=[]
    IYY=[]
    RZX=[]
    IZX=[]
    RZY=[]
    IZY=[]
    
    while (line<len(arquivo_imp)):
        freq=[]
        RXX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IXX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RXX[-1],IXX[-1],freq)
        if(arquivo_imp[line].find("ZXY")==-1):
            print "ERRO! não encontrou o rótulo ZXY esperado na linha %d\n%s" % (line,arquivo_imp[line])
        RXY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IXY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RXY[-1],IXY[-1],freq)
        if(arquivo_imp[line].find("ZYX")==-1):
            print "ERRO! não encontrou o rótulo ZYX esperado na linha %d\n%s" % (line,arquivo_imp[line])
        RYX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IYX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RYX[-1],IYX[-1],freq)
        if(arquivo_imp[line].find("ZYY")==-1):
            print "ERRO! não encontrou o rótulo ZYY esperado na linha %d\n%s" % (line,arquivo_imp[line])
        RYY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IYY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RYY[-1],IYY[-1],freq)
        if(arquivo_imp[line].find("ZZX")==-1):
            print "ERRO! não encontrou o rótulo ZZX esperado na linha %d\n%s" % (line,arquivo_imp[line])
        RZX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IZX.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RZX[-1],IZX[-1],freq)
        if(arquivo_imp[line].find("ZZY")==-1):
            print "ERRO! não encontrou o rótulo ZZY esperado na linha %d\n%s" % (line,arquivo_imp[line])
        RZY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        IZY.append(numpy.reshape(numpy.array([float(s) for s in range(nX*nY)]),(nX,nY)))
        line=le_bloco(arquivo_imp,line,nX,nY,RZY[-1],IZY[-1],freq)
        # todos os blocos devem ter a mesma freqüência
        for f in range (1,len(freq)):
            if(freq[0]!=freq[f]):
                print freq
                print "ERRO! blocos com freqüências diferentes.\nFim do bloco na linha %d" % (line)
        T.append(1.0/freq[0])
    #print line
    # Todos os dados já foram coletados e armazenados
    # preparar saída dos dados pedidos

    # calcula o valor da metade da ultima camada em Y
    # os Y[] são a posição central dos blocos,
    # de modo que Y[j]-Y[j-1] é a soma da metade de cada bloco contíguo
    # exceto o primeiro, que é apenas a metade do primeiro bloco
    y0=Y[0]
    for j in range(1,nY):
        y0=Y[j]-(Y[j-1]+y0)
    # acrescenta metade da ultima camada em Y para obter referencia correta
    Y.append(Y[-1]+y0)

    if(NT!=0): # vetores de inducao
        NT-=1
        for p in range(0,len(pontos)):
            if(pontos[p][0]==0): # vetores de inducao ao longo da linha y
                j=pontos[p][1]-1
                for i in range(0,nX):
                    azim=(180.0/math.pi)*math.atan2(RZY[NT][i,j],RZX[NT][i,j])+180.0+90.0
                    #azim=(180.0/3.1415927)*atan2(RTZY[i,j,NT],RTZX[i,j,NT])+180+90
                    if(azim>=360):
                        azim-=360
                    if(azim<0):
                        azim+=360
                    mod=math.sqrt(RZX[NT][i,j]**2+RZY[NT][i,j]**2)*e
                    print "%f %f %f %f" % (X[i],Y[-1]-Y[j],azim,mod)
            elif(pontos[p][1]==0): # vetores de inducao ao longo da linha x
                i=pontos[p][0]-1
                for j in range(0,nY):
                    azim=(180.0/math.pi)*math.atan2(RZY[NT][i,j],RZX[NT][i,j])+180.0+90.0
                    if(azim>=360):
                        azim-=360
                    if(azim<0):
                        azim+=360
                    mod=math.sqrt(RZX[NT][i,j]**2+RZY[NT][i,j]**2)*e
                    print "%f %f %f %f" % (X[i],Y[-1]-Y[j],azim,mod)
            else: # vetores de inducao no ponto xy
                i=pontos[p][0]-1
                j=pontos[p][1]-1
                azim=(180.0/math.pi)*math.atan2(RZY[NT][i,j],RZX[NT][i,j])+180.0+90.0
                if(azim>=360):
                    azim-=360
                if(azim<0):
                    azim+=360
                mod=math.sqrt(RZX[NT][i,j]**2+RZY[NT][i,j]**2)*e
                print "%f %f %f %f" % (X[i],Y[nY]-Y[j],azim,mod)
    else: # formato Jones em x y
        for p in range(0,len(pontos)):
            if(pontos[p][0]!=0 and pontos[p][1]!=0):
                arq_saida="x%03dy%03d.dat" % (pontos[p][0],pontos[p][1])
                saida=open(arq_saida,'w')
                i=pontos[p][0]-1
                j=pontos[p][1]-1
        
                line= "# dados sinteticos do modelo %s\n" % (arq_imp)
                saida.write(line)
                line= ">STATION   :x%03dy%03d\n" %  (pontos[p][0],pontos[p][1])
                saida.write(line)
                line= ">AZIMUTH   =          90\n"
                saida.write(line)
                line= ">LATITUDE  =    %f\n" % (Y[-1]-Y[j])
                saida.write(line)
                line= ">LONGITUDE =    %f\n" % (X[i])
                saida.write(line)
                line= ">ELEVATION =           0\n\n"
                saida.write(line)

                saida_jones(saida,"ZXX SI units (ohms)\n",T,RXX,IXX,i,j)
                saida_jones(saida,"ZXY SI units (ohms)\n",T,RXY,IXY,i,j)
                saida_jones(saida,"ZYX SI units (ohms)\n",T,RYX,IYX,i,j)
                saida_jones(saida,"ZYY SI units (ohms)\n",T,RYY,IYY,i,j)
                saida_jones(saida,"TZX\n",T,RZX,IZX,i,j)
                saida_jones(saida,"TZY\n",T,RZY,IZY,i,j)
                saida.close()
