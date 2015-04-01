#!/usr/bin/env python
# coding=utf-8
# usage: selection.py [start-end]
# return start end
import sys

if (__name__=="__main__"):
    selection=sys.argv[1]
    if (selection[0]!="[" and selection[-1]!="]" ):
        print "expected selection in format [start-end]"
        sys.exit(1)

    # remove []
    selection=selection[1:-1]
    selection=selection.split("-")
    if (len(selection)!=2 or not(selection[0].isdigit())
        or not(selection[1].isdigit()) or int(selection[0])>int(selection[1])):
        print "expected selection in format [start-end]"
        print "start and end as integers and start<=end"
        sys.exit(1)
    else:
        print "%s %s" % (selection[0],selection[1])
