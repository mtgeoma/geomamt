#!/usr/bin/env python
# coding=utf-8
# uso: gmt2kml.py arquivo-gmt.dat > arquivo.kml
# arquivo-gmt.dat: longitude latitude identificação
# sendo as coordenadas em ponto decimal

import sys

def escreve_placemark(dados):
    print "\t\t<Placemark>"
    print "\t\t\t<name>%s</name>" % dados[2]
    print "\t\t\t<Point>"
    print "\t\t\t\t<coordinates>%s,%s,0</coordinates>" % (dados[0],dados[1])
    print "\t\t\t</Point>"
    print "\t\t</Placemark>"

if (__name__=="__main__"):
    #escreve cabeçario
    print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
    print "<kml xmlns=\"http://earth.google.com/kml/2.2\">"
    print "\t<Document>"

    arquivo=open(sys.argv[1],'r').read()
    arquivo=arquivo.splitlines()
    for line in range(0,len(arquivo)):
        dados=arquivo[line].split()
        if(dados[0]!=">"):
            escreve_placemark(dados)
    # fim
    print "\t</Document>"
    print "</kml>"
