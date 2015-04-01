#!/usr/bin/env python
# coding=utf-8
import sys
import string

if (__name__=="__main__"):
    if (len(sys.argv)!=2):
        print "usage: %s check_sum" % sys.argv[0]
        print "file must be the output of a sum program (md5sum, sha256sum)"
        print "where first column is sum number and second column is file name"
        print "program will check if there is duplicated sum number in check_sum"
        sys.exit()
    sum={}
    arq=sys.argv[1]
    arq=open(arq,'r').read()
    arq=arq.splitlines()
    for line in arq:
        col=line.split()
        if(sum.has_key(col[0])):
            sum[col[0]].append(col[1])
        else:
            sum[col[0]]=[col[1]]

    for s in sum.keys():
        if(len(sum[s])>1):
            out=""
            for f in sum[s]:
                out+=" %s" % f
            print out
