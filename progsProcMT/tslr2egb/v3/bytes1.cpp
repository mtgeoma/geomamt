#include "bytes1.h"

bool correct_types(void) {
    if(sizeof(int2)==2 && sizeof(int4)==4 && sizeof(int8)==8 && sizeof(float8)==8)
        return true;
    else {
        cout<<"some type haven't the correct size:"<<endl;
        cout<<"size of int2       "<<setw(4)<<sizeof(int2)<<endl;
        cout<<"size of int4       "<<setw(4)<<sizeof(int4)<<endl;
        cout<<"size of float8     "<<setw(4)<<sizeof(float8)<<endl<<endl;

        cout<<"you will need to edit the bytes.h file."<<endl;
        cout<<"the sizes of types in your machine are:"<<endl;
        cout<<"size of short      "<<setw(4)<<sizeof(short)<<endl;
        cout<<"size of int        "<<setw(4)<<sizeof(int)<<endl;
        cout<<"size of long       "<<setw(4)<<sizeof(long)<<endl;
        cout<<"size of long long  "<<setw(4)<<sizeof(long long)<<endl;
        cout<<"size of float      "<<setw(4)<<sizeof(float)<<endl;
        cout<<"size of double     "<<setw(4)<<sizeof(double)<<endl;
        cout<<"size of long double"<<setw(4)<<sizeof(long double)<<endl;
        return false;
    }
}

//#ifndef _BYTESWAP_H    // if didn't have the gcc routine
//short bswap_16(short datum) {
//    return (((datum & 0xff00) >>  8) | ((datum & 0x00ff) <<  8));
//}
//
//long bswap_32(long datum) {
//    return (((datum & 0xff000000) >> 24) | ((datum & 0x00ff0000) >>  8) | \
//            ((datum & 0x0000ff00) <<  8) | ((datum & 0x000000ff) << 24));
//}
//#endif
