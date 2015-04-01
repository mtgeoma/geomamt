#include "count.h"
#include <stdio.h>

timecount startcount(string instr, table &tbl) {
    timecount clock;

    // open MTU-LR file
    ifstream in(instr.c_str(), ios::in | ios::binary);
    if(in.fail()) {
        cout<<"can't open the MTU-LR file "<<instr<<endl;
        exit(1);
    }

    // starttime
    amxtds  STIM;
    in.read(&STIM,sizeof(amxtds));

    if(STIM.cen) {
        clock.swin.jul=cal_to_jul(STIM.yr+100*STIM.cen,STIM.mn,STIM.dy);
        clock.swin.min=STIM.min+STIM.hr*60;
    }
    else {  // 0 if time/date is incorrect
        cout<<"First time/date is incorrect in ts file "<<instr<<endl
            <<"give the correct start time of this file:"<<endl;
        string resp;
        double julian;
        Dates_format recognized;

        for(;;) { // forever
            getline(cin,resp);
            if(parse_date(resp.c_str(),STIM.cen,50,FMT_nohint,&julian,&recognized)==EXIT_SUCCESS) {
                int year, month, day, hour, min;
                double sec;
                jul_to_cal_and_time(julian,0.5,&year,&month,&day,&hour,&min,&sec);
                clock.swin.jul=cal_to_jul(year,month,day);
                clock.swin.min=min+hour*60;
                break;
            }
            else {
                cout<<"can't recognize the date format"<<endl
                    <<"please, enter again with the correct start time:"<<endl;
            }
        }
    }

    unsigned short sernum;
    in.read(&sernum,2);
    if(bswap)
        sernum=bswap_16(sernum);

    // check with the MTU serial number in tbl file
    if(tbl.snum>0) { // if it was found a serial number in tbl file
        if(tbl.snum!=sernum) {
            cerr<<"FATAL ERROR: diferent MTU serial number:"<<endl
                <<"in ts MTU-LR file : "<<sernum<<endl
                <<"in tbl file       : "<<tbl.snum<<endl;
            exit(1);
        }
    }
    else
        tbl.snum=sernum;

    unsigned short scan;
    in.read(&scan,2);
    if(bswap)
        scan=bswap_16(scan);

    // check with the srpm in tbl file
    if(tbl.srpm>0) { // if it was found a srpm in tbl file
        if(tbl.srpm!=scan) {
            cerr<<"WARNING: diferent sample rate per minute:"<<endl
                <<"in ts MTU-LR file : "<<scan<<endl
                <<"in tbl file       : "<<tbl.srpm<<endl
                <<"will use the value in MTU-LR file"<<endl;
            tbl.srpm=scan;
        }
    }
    else
        tbl.srpm=scan;

    clock.srpm=tbl.srpm;

    unsigned char chan;
    in.read(&chan,1);
    if(chan!=5) {
        cerr<<chan<<" channels! Work only with 5 channels"<<endl;
        exit(1);
    }
    tbl.nchan=chan;
    clock.nchan=chan;

    unsigned char taglen;
    in.read(&taglen,1);
    if(!taglen)
        taglen=16u;
    clock.tagsize=taglen;

    // size of a record
    unsigned long recsize=taglen+chan*scan*3u;

    // size of file
    struct stat results;
    if(stat(instr.c_str(),&results)) {
        cerr<<"can't check the size of file "<<instr<<endl;
        exit(1);
    }

    // endtime
    amxtds  ETIM;
    in.seekg(results.st_size-recsize);
    in.read(&ETIM,sizeof(amxtds));
    if(ETIM.cen) {
        clock.ewin.jul=cal_to_jul(ETIM.yr+100*ETIM.cen,ETIM.mn,ETIM.dy);
        clock.ewin.min=ETIM.min+1+ETIM.hr*60;
    }
    else {  // get the end time from the size of file
        long nrec=results.st_size/recsize;
        clock.ewin.jul=clock.swin.jul+nrec/1440L;  // 1 day = 1440 minutes
        clock.ewin.min=clock.swin.jul+nrec%1440L;
    }

    // ask for start time
    string resp;
    Dates_format recognized;
    int year, month, day, hour, min;
    double julian;

    cout<<"Give required start time:"<<endl;
    for(;;) { // forever
        jul_to_cal(clock.swin.jul,&year,&month,&day);
        hour=clock.swin.min/60;
        min=clock.swin.min%60;

        printf("[%04u-%02u-%02uT%02u:%02u:00]\n",year,month,day,hour,min);
        getline(cin,resp);
        if(resp.size()==0) { // choose default option
            clock.sdata.jul=clock.swin.jul;
            clock.sdata.min=clock.swin.min;
            break;
        }
        else if(parse_date(resp.c_str(),STIM.cen,50,FMT_nohint,&julian,&recognized)==EXIT_SUCCESS) {
            double sec;
            jul_to_cal_and_time(julian,0.5,&year,&month,&day,&hour,&min,&sec);
            clock.sdata.jul=cal_to_jul(year,month,day);
            clock.sdata.min=min+hour*60;

            // check if sdata is equal or greater than swin
            if(clock.sdata.jul>clock.swin.jul)
                break;
            else if (clock.sdata.jul==clock.swin.jul && clock.sdata.min>=clock.swin.min)
                break;
            else
                cout<<"required start time can't be before the start windows time"<<endl;
        }
        else
            cout<<"can't recognize the date format"<<endl;

        cout<<"please, enter again with the required start time:"<<endl;
    }

    // ask for end time
    cout<<"Give required end time:"<<endl;
    for(;;) { // forever
        jul_to_cal(clock.ewin.jul,&year,&month,&day);
        hour=clock.ewin.min/60;
        min=clock.ewin.min%60;

        printf("[%04u-%02u-%02uT%02u:%02u:00]\n",year,month,day,hour,min);
        getline(cin,resp);
        if(resp.size()==0) { // choose default option
            clock.edata.jul=clock.ewin.jul;
            clock.edata.min=clock.ewin.min;
            break;
        }
        else if(parse_date(resp.c_str(),STIM.cen,50,FMT_nohint,&julian,&recognized)==EXIT_SUCCESS) {
            double sec;
            jul_to_cal_and_time(julian,0.5,&year,&month,&day,&hour,&min,&sec);
            clock.edata.jul=cal_to_jul(year,month,day);
            clock.edata.min=min+hour*60;

            // check if edata is greater than sdata
            if(clock.edata.jul<clock.sdata.jul ||
               (clock.edata.jul==clock.sdata.jul && clock.edata.min<=clock.sdata.min))
                cout<<"required end time can't be before start time"<<endl;
            // check if edata is equal or less than ewin
            else if(clock.edata.jul<clock.ewin.jul)
                break;
            else if (clock.edata.jul==clock.ewin.jul && clock.edata.min<=clock.ewin.min)
                break;
            else
                cout<<"required end time can't be after the end windows time"<<endl;
        }
        else
            cout<<"can't recognize the date format"<<endl;

        cout<<"please, enter again with the required end time:"<<endl;
    }

    // setup counters: for(long i=scount; i<ecount; i++)
    // number of first record to read
    clock.scount=(clock.sdata.jul-clock.swin.jul)*1440L;
    clock.scount+=clock.sdata.min-clock.swin.min;

    // number of last record to read
    clock.ecount=(clock.edata.jul-clock.swin.jul)*1440L;
    clock.ecount+=clock.edata.min-clock.swin.min;
/*
    cout<<endl<<"print out some tests"<<endl;
    printf("win start time : %04u-%02u-%02uT%02u:%02u:%02u\n",
            STIM.yr+STIM.cen*100,STIM.mn,STIM.dy,STIM.hr,STIM.min,STIM.sec);

    jul_to_cal(clock.sdata.jul,&year,&month,&day);
    hour=clock.sdata.min/60;
    min=clock.sdata.min%60;
    printf("data start time: %04d-%02d-%02dT%02d:%02d:00\n",
            year,month,day,hour,min);

    printf("win end time   : %04u-%02u-%02uT%02u:%02u:%02u\n",
            ETIM.yr+ETIM.cen*100,ETIM.mn,ETIM.dy,ETIM.hr,ETIM.min+1,ETIM.sec);

    jul_to_cal(clock.edata.jul,&year,&month,&day);
    hour=clock.edata.min/60;
    min=clock.edata.min%60;
    printf("data end time  : %04d-%02d-%02dT%02d:%02d:00\n\n",
            year,month,day,hour,min);

    cout<<"first record : "<<clock.scount<<endl
        <<"last record  : "<<clock.ecount<<endl;

    cout<<"sizeof record: "<<recsize<<endl
        <<"sizeof file  : "<<results.st_size<<endl
        <<"num record   : "<<results.st_size/recsize<<endl
        <<clock.ewin.jul<<" : "<<clock.swin.jul<<" : "<<clock.ewin.min<<" : "<<clock.swin.min<<endl
        <<"num minutes  : "<<(clock.ewin.jul-clock.swin.jul)*1440+(clock.ewin.min-clock.swin.min)<<endl;
*/
    return clock;
}
