#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

typedef struct
{                                      /* AMX Time/Date Structure             */
      unsigned char amtdsec;           /* seconds  (0-59)                     */
      unsigned char amtdmin;           /* minutes  (0-59)                     */
      unsigned char amtdhr;            /* hours    (0-23)                     */
      unsigned char amtddy;            /* day      (1-31)                     */
      unsigned char amtdmn;            /* month    (1-12)                     */
      unsigned char amtdyr;            /* year     (0-99)                     */
      unsigned char amtdow;            /* day of week (Mon=1 to Sun=7)        */
      unsigned char amtdcen;           /* 0 if time/date is incorrect         */
                                       /* century if time/date is correct     */
} amxtds;

struct tblinfo {
    int snum;    // serial number of mtu unit
    string site; // site ID of the measurement
    amxtds stm1; // window1 start time
    amxtds etm1; // window1 end time
    int srw1;    // sample rate for window1
    double exln; // Ex dipole length, m.
    double eyln; // Ey dipole length, m.
};


// put the file in position men, check if the label is correct
// and then put the file's pointer ready to read his field
void positioning(long men, char* label);

ifstream ifile;
int main(int argc, char* argv[])
{
    ifile.open(argv[1], ios::in | ios::binary);
    if(ifile.fail()) {
    	cerr<<"Can't open file "<<argv[1]<<endl;
    	exit(1);
    }

    tblinfo header;
    char strbuffer[13];

    // ready to get SNUM
    positioning(4L*0x19, "SNUM");
    ifile.read(&header.snum,sizeof(int));

    // ready to get SITE
    positioning(7L*0x19, "SITE");
    ifile.read(strbuffer,13);
    header.site=strbuffer;

    // ready to get SITE
    positioning(13L*0x19, "STM1");
    ifile.read(&header.stm1,sizeof(header.stm1));

    // ready to get exln
    positioning(39L*0x19, "EXLN");
    ifile.read(&header.exln,sizeof(header.exln));

    // ready to get eyln
    positioning(40L*0x19, "EYLN");
    ifile.read(&header.eyln,sizeof(header.eyln));

    cout<<"SNUM="<<header.snum<<endl;
    cout<<"SITE="<<header.site<<"(len="<<header.site.length()<<")"<<endl;
    cout<<"STM1="<<setw(2)<<(int)header.stm1.amtddy
        <<"/"<<setw(2)<<(int)header.stm1.amtdmn
        <<"/"<<setw(2)<<(int)header.stm1.amtdyr
        <<"T"<<setw(2)<<(int)header.stm1.amtdhr
        <<":"<<setw(2)<<(int)header.stm1.amtdmin
        <<":"<<setw(2)<<(int)header.stm1.amtdsec
        <<endl;
    cout<<"EXLN="<<header.exln<<endl;
    cout<<"EYLN="<<header.eyln<<endl;

    ifile.close();
    return 0;
}

void positioning(long men, char* label)
{
    ifile.seekg(men);

    char labeltmp[5];
    char skip[7];

    ifile.read(labeltmp,5);
    if(strcmp(labeltmp,label)) {
    	cerr<<"it was spected \""<<label<<"\", not "<<labeltmp<<endl;
    	exit (1);
    }

    ifile.read(skip,7);
}
