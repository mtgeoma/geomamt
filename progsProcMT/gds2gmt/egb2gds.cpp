#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <complex>
#include <cmath>
#include <iomanip>

int main(int argc, char* argv[])
{
	
    string instr;
    double rotate=0.0;
    bool decl=false;
    enum output_data {Rzx,Izx,Rzy,Izy,Azx,Pzx,Azy,Pzy,Riv,Iiv};
    enum output_data output;
    
    if (argc<3 || argc>4) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -I<input data> -D<output data> [-R<angle>]"<<endl
            <<"<input data>: egbert file (ascii)"<<endl
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

    // open input file
    ifstream ifile(instr.c_str());
    if(ifile.fail()) {
        cout<<"can't open the "<<instr<<" file"<<endl;
        exit(1);
    }

    cout.setf(ios::scientific);
    cout.precision(4);

    string line;
    size_t ln;
    // read header
    // look for declination
    if(decl) {
        while (getline(ifile,line))
 	    if((ln=line.find("declination "))!=string::npos) break;
        string declination(line,line.length()-8,8);
        rotate=atof(declination.c_str());
        if(rotate>360 || rotate<-360) {
            cerr<<"declination \""<<declination<<"\" in the Egbert file must be between -360 and 360 degrees"<<endl;
            exit(1);
        }
    }
	
    // look for number of channels
    while (getline(ifile,line))
	if((ln=line.find("number of channels "))!=string::npos) break;
	
    // check if no. of channels==3
    string channels(line,ln+19,3);
    if(atoi(channels.c_str())!=3) {
	cerr<<"actualy this program is restricted to 3 channels"<<endl;
	exit(1);
    }
	
    int nT;
    double T;
    string nfrequencies(line,line.find("number of frequencies ")+22,4);
    if((nT=atoi(nfrequencies.c_str()))<=0) {
	cerr<<"reading number of frequencies"<<endl;
	exit(1);
    }
	
    // start to read TF, S and N
	size_t n=0;
	while (getline(ifile,line)) {
		if((ln=line.find("period :"))==string::npos) continue;
		
		string strperiod(line,ln+8,12);
		if((T=atof(strperiod.c_str()))<=0) {
			cerr<<n+1<<" element with T<=0"<<endl;
			exit(1);
		}
		
		while(getline(ifile,line))
			if(line.find("Transfer Functions")!=string::npos) break;
		
		
		double Tzxre, Tzxim, Tzyre, Tzyim;
		getline(ifile,line);
		istrstream istrline(line.c_str());
		istrline>>Tzxre>>Tzxim>>Tzyre>>Tzyim;
		
		getline(ifile,line);
		if(line.find("Signal Power Matrix")==string::npos) {
			cerr<<"reading data of period "<<T<<endl;
			exit(1);
		}
		
		// Power Matrix
		complex<double> S[3];
		for (size_t i=0;i<2;i++) {
			getline(ifile,line);
			istrstream istrline(line.c_str());
	    	for (size_t j=0;j<=i;j++) {
				double real, imag;
				istrline>>real>>imag;
				S[i+j]=complex<double>(real,imag);
			}
		}
		
		getline(ifile,line);
		if(line.find("Residual Covar")==string::npos) {
			cerr<<"reading data of period "<<T<<endl;
			exit(1);
		}
		
		double Tzxerr, Tzyerr;
		getline(ifile,line);
		
		{  // problemas de escopo
		istrstream istrline(line.c_str());
		double real, imag;
		istrline>>real>>imag;
		Tzxerr=real*S[0].real();
		Tzyerr=real*S[2].real();
		}

                double A, P, Aerr, Perr;
		switch(output) {
		case Rzx:
		        cout<<setw(14)<<T<<setw(14)<<Tzxre<<setw(14)<<Tzxerr<<endl;
		        break;
		case Izx:
		        cout<<setw(14)<<T<<setw(14)<<Tzxim<<setw(14)<<Tzxerr<<endl;
		        break;
		case Rzy:
		        cout<<setw(14)<<T<<setw(14)<<Tzyre<<setw(14)<<Tzyerr<<endl;
		        break;
		case Izy:
		        cout<<setw(14)<<T<<setw(14)<<Tzyim<<setw(14)<<Tzyerr<<endl;
		        break;
		case Azx:
		        A=sqrt(Tzxre*Tzxre+Tzxim*Tzxim);
		        Aerr=Tzxerr;
		        cout<<setw(14)<<T<<setw(14)<<A<<setw(14)<<Aerr<<endl;
		        break;
		case Pzx:
		        P=(180.0/M_PI)*atan2(Tzxim,Tzxre);
		        Perr=(180.0/M_PI)*Tzxerr/sqrt(Tzxre*Tzxre+Tzxim*Tzxim);
		        cout<<setw(14)<<T<<setw(14)<<P<<setw(14)<<Perr<<endl;
		        break;
		case Azy:
		        A=sqrt(Tzyre*Tzyre+Tzyim*Tzyim);
		        Aerr=Tzyerr;
		        cout<<setw(14)<<T<<setw(14)<<A<<setw(14)<<Aerr<<endl;
		        break;
		case Pzy:
		        P=(180.0/M_PI)*atan2(Tzyim,Tzyre);
		        Perr=(180.0/M_PI)*Tzyerr/sqrt(Tzyre*Tzyre+Tzyim*Tzyim);
		        cout<<setw(14)<<T<<setw(14)<<P<<setw(14)<<Perr<<endl;
		        break;
		case Riv:
		        A=sqrt(Tzxre*Tzxre+Tzyre*Tzyre);
		        Aerr=sqrt(Tzxre*Tzxre*Tzxerr*Tzxerr+Tzyre*Tzyre*Tzyerr*Tzyerr)/A;
		        // note that adding 180 degrees, it invert the direction and put between 0 and 360 degrees
		        P=(180.0/M_PI)*atan2(Tzyre,Tzxre)+180.0;
		        // rotate
		        if(rotate!=0) {
		            if(P+rotate<0)
		                P+=rotate+360.;
		            else if(P+rotate>360)
		                P+=rotate-360.;
		            else
		                P+=rotate;
		        }
		        Perr=(180.0/M_PI)*sqrt(Tzyre*Tzyre*Tzxerr*Tzxerr+Tzxre*Tzxre*Tzyerr*Tzyerr)/(Tzxre*Tzxre+Tzyre*Tzyre);
		        cout<<setw(14)<<T<<setw(14)<<P<<setw(14)<<A<<setw(14)<<Perr<<setw(14)<<Aerr<<endl;
		        break;
		case Iiv:
		        A=sqrt(Tzxim*Tzxim+Tzyim*Tzyim);
		        Aerr=sqrt(Tzxim*Tzxim*Tzxerr*Tzxerr+Tzyim*Tzyim*Tzyerr*Tzyerr)/A;
		        P=(180.0/M_PI)*atan2(Tzyim,Tzxim);
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
		        Perr=(180.0/M_PI)*sqrt(Tzyim*Tzyim*Tzxerr*Tzxerr+Tzxim*Tzxim*Tzyerr*Tzyerr)/(Tzxim*Tzxim+Tzyim*Tzyim);
		        cout<<setw(14)<<T<<setw(14)<<P<<setw(14)<<A<<setw(14)<<Perr<<setw(14)<<Aerr<<endl;
		        break;
		}

	}
	if(!ifile.eof())
		cerr<<"Warning! The file \""<<instr<<"\" was not read until the end"<<endl;
	ifile.close();
}
