#!/usr/bin/env python
# coding=utf-8
# uso: calc-data.py data-formato-iso-ou-european intervalo-de-dias
import sys
from datetime import date
from datetime import timedelta

def calc_pascoa(ano):
    # ver na wikipedia: computus
    a = ano%19
    b = ano/100
    c = ano%100
    d = b/4
    e = b%4
    f = (b+8)/25
    g = (b - f + 1)/3
    h = (19*a+b-d-g+15)%30
    i = c/4
    k = c%4
    L = (32+2*e+2*i-h-k)%7
    m = (a+11*h+22*L)/451

    mes = (h+L-7*m+114)/31
    dia = ((h+L-7*m+114)%31)+1
    return date(ano,mes,dia)

def feriado(data,pascoa):
    # Feriados Nacionais com data fixa sem exceções
    # 01/01 Confraternização Universal
    # 21/04 Tiradentes
    # 01/05 Dia do Trabalho
    # 07/09 Independência do Brasil
    # 12/10 Padroeira do Brasil
    # 02/11 Finados
    # 15/11 Proclamação da República
    # 25/12 Natal
    if ( (data.day==1  and data.month==1 ) or
         (data.day==1  and data.month==5 ) or
         (data.day==21 and data.month==4 ) or
         (data.day==7  and data.month==9 ) or
         (data.day==12 and data.month==10) or
         (data.day==2  and data.month==11) or
         (data.day==15 and data.month==11) or
         (data.day==25 and data.month==12) ):
        return True
    # Feriados Nacionais com data móvel
    # Segunda de Carnaval (pascoa - 48)
    # Terça   de Carnaval (pascoa - 47)
    # Paixao              (pascoa -  2)
    # Corpus Christi      (pascoa + 60)
    if ( (data == pascoa+timedelta(days=-48)) or
         (data == pascoa+timedelta(days=-47)) or
         (data == pascoa+timedelta(days= -2)) or
         (data == pascoa+timedelta(days=+60)) ):
        return True

def dias_uteis(datai,dataf):
    datas=[]
    if (datai.find("-")!=-1 and dataf.find("-")!=-1):
        datai=datai.split("-")
        datai=date(int(datai[0]),int(datai[1]),int(datai[2]))
        dataf=dataf.split("-")
        dataf=date(int(dataf[0]),int(dataf[1]),int(dataf[2]))
        data=datai
        ano=data.year
        pascoa=calc_pascoa(ano)
        while (data<=dataf):
            if( (data.weekday()<5) and (not feriado(data,pascoa)) ):
                datas.append(data.isoformat())
            data=data+timedelta(days=1)
            if (data.year!=ano):
                ano=data.year
                pascoa=calc_pascoa(ano)
    elif (datai.find("/")!=-1 and dataf.find("/")!=-1):
        datai=datai.split("/")
        datai=date(int(datai[2]),int(datai[1]),int(datai[0]))
        dataf=dataf.split("/")
        dataf=date(int(dataf[2]),int(dataf[1]),int(dataf[0]))
        data=datai
        ano=data.year
        pascoa=calc_pascoa(ano)
        while (data<=dataf):
            if( (data.weekday()<5) and (not feriado(data,pascoa)) ):
                datas.append(data.strftime("%d/%m/%Y"))
            data=data+timedelta(days=1)
            if (data.year!=ano):
                ano=data.year
                pascoa=calc_pascoa(ano)
    return datas

if (__name__=="__main__"):
    datai=sys.argv[1]
    dataf=sys.argv[2]
    datas=dias_uteis(datai,dataf)
    for data in datas:
        print data
