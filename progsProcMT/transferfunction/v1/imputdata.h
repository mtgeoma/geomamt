/***************************************************************************
                          imputdata.h  -  description
                             -------------------
    begin                : Sat Jan 13 2001
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

#ifndef IMPUTDATA_H
#define IMPUTDATA_H


/**
  *@author Marcelo Banik de Pádua
  */

#include "transferfunction.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <string>

// actualy don't treat missing data
namespace DataType {

enum type_data {jones, egbert, emprog, pt1, edi};

class ImputData {
	type_data itype;
	string ifile_name;
	ifstream ifile;
	vector<TransferFunction> VTF;

	void jones_format();
	void egbert_format();
	void emprog_format();
	void pt1_format();
        vector<double> getpt1data(const char* target, int& pfirst, int& plast);
	void edi_format();
        vector<double> getedidata(const char* target, int& ndata);
public:
	ImputData(string fname, type_data type);
	~ImputData();
	void getdata(vector<TransferFunction>& VTFo);
};
} //end of namespace DataType
#endif
