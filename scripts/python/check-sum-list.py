#!/usr/bin/env python
# coding=utf-8
import sys
import string

if (__name__=="__main__"):
    if (len(sys.argv)!=3):
        print "usage: %s ref_sum check_sum" % sys.argv[0]
        print "both file must be the output of a sum program (md5sum, sha256sum)"
        print "where first column is sum number and second column is file name"
        print "program will check if there is sum number in check_sum that"
        print "is not in ref_sum"
        sys.exit()
    sum={}
    ref=sys.argv[1]
    ref=open(ref,'r').read()
    ref=ref.splitlines()
    for line in ref:
        col=line.split()
        if(sum.has_key(col[0])):
            sum[col[0]].append(col[1])
        else:
            sum[col[0]]=[col[1]]

    check=sys.argv[2]
    check=open(check,'r').read()
    check=check.splitlines()
    for line in check:
        col=line.split()
        if(sum.has_key(col[0])==0):
            print "%s only in %s" % (col[1],sys.argv[2])
