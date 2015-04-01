/***************************************************************************
                          commandline.cpp  -  description
                             -------------------
    begin                : Tue Jan 16 2001
    copyright            : (C) 2001 by Marcelo Banik de Pádua
    email                : banik@dge.inpe.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string.h>
#include "commandline.h"

CommandLine::CommandLine(int argc, char *argv[]) {
	if(argc<4) {
		cerr<<"error in parameters"<<endl;
		exit(1);
	}

	// defaults
	error_prop=TransferFunction::stodt;
	width_fmt=14;
	precision_fmt=8;
	teta=0.;

	string str;
	for(int i=1; i<argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			case 'A': // angle to rotate
				teta=atof(argv[i]+2);
				break;
			case 'I': // induction vector
			    if(argv[i][2]=='V')
			    	out_data=ind_vectors;
				else {
					cerr<<'\''<<str<<'\''
					    <<" is a unknown option"<<endl;
					exit(1);
				}
				break;
			case 'F': // file input format
				if(strstr(argv[i]+2,"jones")!=0)
					format=jones;
				else if(strstr(argv[i]+2,"egbert")!=0)
					format=egbert;
				else if(strstr(argv[i]+2,"emprog")!=0)
					format=emprog;
				else if(strstr(argv[i]+2,"pt1")!=0)
					format=pt1;
				else if(strstr(argv[i]+2,"edi")!=0)
					format=edi;
				else {
					cerr<<'\''<<str<<'\''
					    <<" is a unknown input format in option -F"<<endl;
					exit(1);
				}
				break;
			case 'D':
				if(argv[i][2]=='z')
					out_data=impedance;
				else if(argv[i][2]=='r')
					out_data=rhophi;
				else if(argv[i][2]=='m')
					out_data=mackie;
				else if(argv[i][2]=='T') {
				    out_data=dist_types;
				    break;
				}
				else {
					cerr<<'\''<<argv[i]+2<<'\''
					    <<" is a unknown input format in option -D"<<endl;
					exit(1);
				}

				if(argv[i][3]=='x')
					r=TransferFunction::x;
				else if(argv[i][3]=='y')
					r=TransferFunction::y;
				else if(argv[i][3]=='z')
					r=TransferFunction::z;
				else {
					cerr<<'\''<<argv[i]+2<<'\''
					    <<" is a unknown input format in option -D"<<endl;
					exit(1);
				}

				if(argv[i][4]=='x')
					c=TransferFunction::x;
				else if(argv[i][4]=='y')
					c=TransferFunction::y;
				else {
					cerr<<'\''<<argv[i]+2<<'\''
					    <<" is a unknown input format in option -D"<<endl;
					exit(1);
				}
				
				break;
			case 'S':
				out_data=strikes;
				break;
			case 'P':
				out_data=phasesplit;
				break;
			case 'E': //stodt, gamble, chave, jones
				if(strstr(argv[i]+2,"jones")!=0)
					error_prop=TransferFunction::jones;
				else if(strstr(argv[i]+2,"chave")!=0)
					error_prop=TransferFunction::chave;
				else if(strstr(argv[i]+2,"gamble")!=0)
					error_prop=TransferFunction::gamble;
				else if(strstr(argv[i]+2,"stodt")!=0)
					error_prop=TransferFunction::stodt;
				else {
					cerr<<'\''<<argv[i]+2<<'\''
					    <<" is a unknown input format in option -E"<<endl;
					exit(1);
				}
			}
		}
		else if(file_name.size()==0) {
			file_name=argv[i];
		}
		else {
			cerr<<"error in parameters"<<endl;
			exit(1);
		}

	}
	
}

CommandLine::~CommandLine(){
}

string CommandLine::file() {
	return file_name;
}

type_data CommandLine::file_format() {
	return format;
}

CommandLine::output_data CommandLine::extract() {
	return out_data;
}

TransferFunction::prop_error CommandLine::prop_error(void) {
	return error_prop;
}

TransferFunction::index CommandLine::row() {
	return r;
}

TransferFunction::index CommandLine::col() {
	return c;
}

int CommandLine::width() {
	return width_fmt;
}

int CommandLine::precision() {
	return precision_fmt;
}

double CommandLine::rot() {
	return teta;
}
