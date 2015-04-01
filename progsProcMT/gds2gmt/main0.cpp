#include <iostream>
//#include <fstream>
//#include <strstream>
#include <cmath>
#include <iomanip>
#include "gds.h"

int main(int argc, char* argv[])
{
    string instr;
    double rotate=0.0;
    bool decl=false;
    enum format_data {egbert,edi,jones};
    enum format_data format;
    enum output_data {Rzx,Izx,Rzy,Izy,Azx,Pzx,Azy,Pzy,Riv,Iiv};
    enum output_data output;
    
    if (argc<4 || argc>5) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -I<input data> -F<format data> -D<output data> [-R<angle>]"<<endl
            <<"<input data>: data file (ascii)"<<endl
            <<"<format data> could be: egbert, edi or jones"<<endl
            <<"<output data> could be:"<<endl
            <<"         Rzx: real part of Tzx"<<endl
            <<"         Izx: imag part of Tzx"<<endl
            <<"         Rzy: real part of Tzy"<<endl
            <<"         Izy: imag part of Tzy"<<endl
            <<"         Azx: amplitude of Tzx"<<endl
            <<"         Pzx: phase(deg)of Tzx"<<endl
            <<"         Azy: amplitude of Tzy"<<endl
            <<"         Pzy: phase(deg)of Tzy"<<endl
            <<"         in the format: period value standard_deviation"<<endl<<endl
            <<"         Riv: real induction vector (inverted)"<<endl
            <<"         Iiv: imag induction vector"<<endl
            <<"         in the format: period azimuth(0 to 360 from north to east)"<<endl
            <<"                        amplitude and errors"<<endl<<endl
            <<"<angle>: declination angle to rotate all vectors [-360;360]."<<endl
            <<"         If instead to choose a number, you use de word decl (-Rdecl),"<<endl
            <<"         the program will get the declination from the file."<<endl
            <<"         WARNING! use this angle as a declination, ie, is the angle between"<<endl
            <<"                  the new coordinate system (geographic north) and"<<endl
            <<"                  the old coordinate system (magnetic north)"<<endl<<endl;
        exit (1);
    }
    else {  // get the parameters
        for(int i=1; i< argc; i++) {
            if(argv[i][0]=='-' && argv[i][1]=='I')
                instr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='F') {
                string formatstr=argv[i]+2;
                if(formatstr=="egbert") format=egbert;
                else if(formatstr=="edi") format=edi;
                else if(formatstr=="jones") format=jones;
                else {
                    cerr<<"unknow format data:"<<formatstr<<endl
                        <<"run the program without the options to see the possible format data"<<endl;
                    exit(1);
                }
            }
            else if(argv[i][0]=='-' && argv[i][1]=='D') {
                string datastr=argv[i]+2;
                if(datastr=="Rzx") output=Rzx;
                else if(datastr=="Izx") output=Izx;
                else if(datastr=="Rzy") output=Rzy;
                else if(datastr=="Izy") output=Izy;
                else if(datastr=="Azx") output=Azx;
                else if(datastr=="Pzx") output=Pzx;
                else if(datastr=="Azy") output=Azy;
                else if(datastr=="Pzy") output=Pzy;
                else if(datastr=="Riv") output=Riv;
                else if(datastr=="Iiv") output=Iiv;
                else {
                    cerr<<"unknow output data:"<<datastr<<endl
                        <<"run the program without the options to see the possible output data"<<endl;
                    exit(1);
                }
            }
            else if(argv[i][0]=='-' && argv[i][1]=='R') {
                string anglestr=argv[i]+2;
                if(anglestr=="decl") {
                    decl=true;
                }
                else {
                    rotate=atof(anglestr.c_str());
                    if(rotate>360 || rotate<-360) {
                        cerr<<"angle \""<<anglestr<<"\" in option -R must be between -360 and 360 degrees"<<endl;
                        exit(1);
                    }
                    if(rotate==0) {
                        cerr<<"angle \""<<anglestr<<"\" in option -R or is zero (therefore not necessary) or is not a number"<<endl;
                        exit(1);
                    }
                }
            }
            else {
                cerr<<argv[i]<<" is an unknow parameter"<<endl;
                exit(1);
            }
        }
    }
    
    vector<tipper> Tz;
    switch(format) {
    case egbert:
            Tz=read_egbert(instr); 
            break;
   case edi:
            Tz=read_edi(instr);
            break;
   case jones:
            Tz=read_jones(instr);
            break;
    }
    
    // write to the standard output
    for(int i=0;i<Tz.size();i++) {
        cout.setf(ios::scientific);
        cout.precision(4);
        double A, P, Aerr, Perr;
        switch(output) {
        case Rzx:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<Tz[i].x.real()<<setw(14)<<sqrt(Tz[i].xvar)<<endl;
                break;
        case Izx:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<Tz[i].x.imag()<<setw(14)<<sqrt(Tz[i].xvar)<<endl;
                break;
        case Rzy:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<Tz[i].y.real()<<setw(14)<<sqrt(Tz[i].yvar)<<endl;
                break;
        case Izy:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<Tz[i].y.imag()<<setw(14)<<sqrt(Tz[i].yvar)<<endl;
                break;
        case Azx:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<abs(Tz[i].x)<<setw(14)<<sqrt(Tz[i].xvar)<<endl;
                break;
        case Pzx:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<(180.0/M_PI)*arg(Tz[i].x)<<setw(14)
                    <<(180.0/M_PI)*sqrt(Tz[i].xvar)/abs(Tz[i].x)<<endl;
                break;
        case Azy:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<abs(Tz[i].y)<<setw(14)<<sqrt(Tz[i].yvar)<<endl;
                break;
        case Pzy:
                cout<<setw(14)<<Tz[i].T<<setw(14)<<(180.0/M_PI)*arg(Tz[i].y)<<setw(14)
                    <<(180.0/M_PI)*sqrt(Tz[i].yvar)/abs(Tz[i].y)<<endl;
                break;
        case Riv:
                A=sqrt(pow(Tz[i].x.real(),2)+pow(Tz[i].y.real(),2));
                Aerr=sqrt(pow(Tz[i].x.real(),2)*Tz[i].xvar+pow(Tz[i].y.real(),2)*Tz[i].yvar)/A;
                // note that adding 180 degrees, it invert the direction and put between 0 and 360 degrees
                P=(180.0/M_PI)*atan2(Tz[i].y.real(),Tz[i].x.real())+180.0;
                // rotate
                if(rotate!=0) {
                    if(P+rotate<0)
                        P+=rotate+360.;
                    else if(P+rotate>360)
                        P+=rotate-360.;
                    else
                        P+=rotate;
                }
                Perr=(180.0/M_PI)*sqrt(pow(Tz[i].y.real(),2)*Tz[i].xvar+pow(Tz[i].x.real(),2)*Tz[i].yvar);
                Perr/=(pow(Tz[i].x.real(),2)+pow(Tz[i].y.real(),2));
                cout<<setw(14)<<Tz[i].T<<setw(14)<<P<<setw(14)<<A<<setw(14)<<Perr<<setw(14)<<Aerr<<endl;
                break;
        case Iiv:
                A=sqrt(pow(Tz[i].x.imag(),2)+pow(Tz[i].y.imag(),2));
                Aerr=sqrt(pow(Tz[i].x.imag(),2)*Tz[i].xvar+pow(Tz[i].y.imag(),2)*Tz[i].yvar)/A;
                P=(180.0/M_PI)*atan2(Tz[i].y.imag(),Tz[i].x.imag());
                if(P<0.) // force to be between 0 and 360
                    P+=360.0;
                // rotate
                if(rotate!=0) {
                    if(P+rotate<0)
                        P+=rotate+360.;
                    else if(P+rotate>360)
                        P+=rotate-360.;
                    else
                        P+=rotate;
                }
                Perr=(180.0/M_PI)*sqrt(pow(Tz[i].y.imag(),2)*Tz[i].xvar+pow(Tz[i].x.imag(),2)*Tz[i].yvar);
                Perr/=(pow(Tz[i].x.imag(),2)+pow(Tz[i].y.imag(),2));
                cout<<setw(14)<<Tz[i].T<<setw(14)<<P<<setw(14)<<A<<setw(14)<<Perr<<setw(14)<<Aerr<<endl;
                break;
        }
    }
}
