/***************************************************************************
                          commandline.h  -  description
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

#ifndef COMMANDLINE_H
#define COMMANDLINE_H


/**
  *@author Marcelo Banik de Pádua
  */

#include <string>
#include "inputdata.h"
#include "transferfunction.h"

using namespace DataType;

class CommandLine {
public:
	enum output_data {impedance,rhophi,skewness,dist_types,ind_vectors,mackie,phasesplit};
private:
	string file_name;
	type_data format;
	output_data out_data;
	TransferFunction::prop_error error_prop;
	TransferFunction::index r;
	TransferFunction::index c;
	int width_fmt;
	int precision_fmt;
	double teta;
public:
	CommandLine(int argc, char *argv[]);
	~CommandLine();
	string file();
	type_data file_format();
	output_data extract();
	TransferFunction::prop_error prop_error(void);
	TransferFunction::index row();
	TransferFunction::index col();
	int width();
	int precision();
	double rot();
};

#endif
