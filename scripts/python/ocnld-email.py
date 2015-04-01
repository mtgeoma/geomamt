#!/usr/bin/env python
# coding=utf-8
import sys
import string
import traceback
import argparse
import smtplib
from email.MIMEText import MIMEText
import os
import re
import getpass

def isEmail(email):
    return re.match(r"^[A-Za-z0-9\.\+_-]+@[A-Za-z0-9\._-]+\.[a-zA-Z]*$", email)

if (__name__ == "__main__"):
    MaxSta = 100
    Subject="Ocean Loading Tides"

    parser = argparse.ArgumentParser(
        description="make an e-mail requesting Ocean Loading Tides\n"
        "(http://www.oso.chalmers.se/~loading/hfo.html)")

    parser.add_argument("FromEmail", metavar="E-MAIL",
                        help="send from this E-MAIL")

    parser.add_argument("StaPos", metavar="sta_pos",help="use stations"
                        " in sta_pos up to a maximum of %d" % MaxSta)

    parser.add_argument("--to-email", dest="ToEmail", metavar="E-MAIL",
                        default="loading@holt.oso.chalmers.se",
                        help="send to this E-MAIL (default=%(default)s)")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--print", dest="send", action="store_false", 
                       default=False, help="print message to stdout"
                       " and subject to stderr (default).")

    group.add_argument("--send", dest="send", action="store_true",
                       help="send e-mail")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--TLS-connection", dest="TLS_connection",
                       action="store_true", default=True,
                       help="send e-mail using a TLS connection (default).")

    group.add_argument("--no-TLS-connection", dest="TLS_connection",
                       action="store_false",
                       help="send e-mail without a TLS connection.")

    parser.add_argument("--host-email", dest="host", metavar="HOST",
                      default="smtp.gmail.com",
                      help="send from this HOST (default=%(default)s)")

    parser.add_argument("--port", dest="port", type=int, default=587,
                      help="send using this PORT (default=%(default)s)")

    args = parser.parse_args()

    if not isEmail(args.FromEmail):
        print "expected an e-mail at first positional argument\n"
        parser.print_help()
        sys.exit(1)

    try:
        StaPos=open(args.StaPos,'r').read()
        StaPos=StaPos.splitlines()

        n = 0
        X = {}
        Y = {}
        Z = {}
        for line in StaPos:
            if( n < MaxSta ):
                # just in case someone follow file_formats/sta_info example:
                # JPLM 1992 07 01 00:00:00.00 1000001.00   -2493304.0630  -4655215.5490   3565497.3390 -3.20000000e-02 1.90000000e-02 6.00000000e-03 Mon Nov  9 15:07:31 PST 1992 itrf91 1992.5
                line=line.replace(':',' ')
                col=line.split()
                sta = col[0]
                # will use only the first sta found
                if( not X.has_key(sta) ):
                    X[sta]=col[8]
                    Y[sta]=col[9]
                    Z[sta]=col[10]
                    n+=1

        text  = "Header = ---- Ocean loading values follow next: ----\n"
        text += "Model = FES2004\n"
        text += "LoadingType = displacement\n"
        text += "CMC = 1\n"
        text += "OutputFormat = BLQ\n"
        text += "Stations = \n"

        for sta in sorted(X):
            text += "%-24s %16s%16s%16s\n" % (sta,X[sta],Y[sta],Z[sta])

        text += "\nMyEmail = %s" % args.FromEmail

        if (args.send):
            msg = MIMEText(text)

            msg['From'] = args.FromEmail
            msg['To'] = args.ToEmail
            msg['Subject'] = Subject


            mailServer = smtplib.SMTP(args.host, args.port)

            if(args.TLS_connection):
                mailServer.ehlo()
                mailServer.starttls()
                mailServer.ehlo()

            Password = getpass.getpass(prompt="password to %s:" % args.FromEmail)
            mailServer.login(args.FromEmail, Password)
            mailServer.sendmail(args.FromEmail, args.ToEmail, msg.as_string())
            # Should be mailServer.quit(), but that crashes...
            mailServer.close()
        else:
            
            sys.stderr.write("send the following message to %s\n"
                             "with Subject: %s\n" % (args.ToEmail,Subject))
            print text

    except:
        (ErrorType,ErrorValue,ErrorTB)=sys.exc_info()
        traceback.print_exc(ErrorTB)
