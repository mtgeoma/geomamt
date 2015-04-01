#!/usr/bin/env python
# coding=utf-8
# uso: calc-data.py data-formato-iso-ou-european intervalo-de-dias

import sys
from datetime import date
from datetime import timedelta

if (__name__=="__main__"):
    data=sys.argv[1]
    delta=int(sys.argv[2])
    if (data.find("-")!=-1):
        data=data.split("-")
        nova_data=date(int(data[0]),int(data[1]),int(data[2]))+timedelta(days=delta)
        print "%s" % (nova_data.isoformat())
    else:
        data=data.split("/")
        nova_data=date(int(data[2]),int(data[1]),int(data[0]))+timedelta(days=delta)
        print "%s" % (nova_data.strftime("%d/%m/%Y"))
