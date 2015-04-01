#include <iostream>
#include <fstream>
#include <string>
#include "bytes2.h"

long byte2V(unsigned char *byte)
{
    unsigned char signal;

    if(byte[2]>=0x80)
        signal=0xFF;
    else
        signal=0x00;

    return ((long) signal<<24 | byte[2]<<16 | byte[1]<<8 | byte[0]);
}

long conv24to32 (long datum) {
    #ifdef WORDS_BIGENDIAN
        cout<<"bigendian"<<endl;
            datum=bswap_32(datum);
    #endif
    datum=datum<<8;
    datum/=256;
    return datum;
}

int main(int argc, char* argv[]) {
    ifstream ifile(argv[1], ios::in | ios::binary);
    if(ifile.fail()) {
        cerr<<"can't open file "<<argv[1]<<endl;
        exit(1);
    }
    unsigned char byte[3];
    unsigned char taglen;
    ifile.seekg(13);
    ifile.read(reinterpret_cast <char*> (&taglen),sizeof(taglen));
    if(taglen==0) taglen=16;

    ifile.seekg(taglen+3);
    ifile.read(reinterpret_cast <char*> (byte),sizeof(byte));

    long datum=byte2V(byte);
    cout<<hex<<datum<<':'<<dec<<datum<<endl;
    short lw= (datum & 0x0000ffff);
    short hw=((datum & 0xffff0000) >> 16);
    datum=(long)(hw<<16 | lw);
    cout<<hex<<datum<<":lw:"<<lw<<":hw:"<<hw<<endl;

    ifile.seekg(taglen+3);
    ifile.read(reinterpret_cast <char*> (&datum),3);
    cout<<hex<<datum<<':';
    datum=datum<<8;
    cout<<hex<<datum/256<<':'<<dec<<datum/256<<endl;

    ifile.seekg(taglen+3);
    ifile.read(reinterpret_cast <char*> (&datum),3);
    cout<<hex<<datum<<':';
    datum=conv24to32(datum);
    cout<<hex<<datum<<':'<<dec<<datum<<endl;
}
