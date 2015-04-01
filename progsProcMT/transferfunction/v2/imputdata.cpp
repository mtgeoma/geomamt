/***************************************************************************
                          imputdata.cpp  -  description
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

#include "imputdata.h"

using namespace DataType;

ImputData::ImputData(string fname, type_data type){
	itype=type;
	ifile_name=fname;
	ifile.open(ifile_name.c_str());
	if(!ifile.is_open()) {
		cerr<<"Can't open the file \""<<ifile_name<<'\"'<<endl;
		exit(1);
	}
}

ImputData::~ImputData(){
	ifile.close();
}

void ImputData::getdata(vector<TransferFunction>& VTFo) {
	switch (itype) {
	case jones:
		jones_format();
		break;
	case emprog:
		emprog_format();
		break;
	case egbert:
		egbert_format();
		break;
	case pt1:
		pt1_format();
                break;
	case edi:
		edi_format();
	}

	VTFo=VTF;
}
