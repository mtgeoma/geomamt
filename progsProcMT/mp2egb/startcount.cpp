#include "count.h"
#include <stdio.h>

timecount startcount(table &tbl, bool drop) {
    timecount clock;

    // get information from tbl
    clock.srpm=tbl.srpm;
    clock.nchan=tbl.nchan;
    clock.sizehdr=tbl.sizehdr;
    
    // calculate century (will be useful to parse_date function)
    int century=tbl.stime.year/100;
    if(!tbl.stime.year%100)
        century+=1;

    // starttime
    clock.swin.jul=cal_to_jul(tbl.stime.year, tbl.stime.month, tbl.stime.day);
    clock.swin.min=tbl.stime.min+tbl.stime.hour*60;
    clock.swin.sec=0;

    //size of a sample (without auxiliary channels)
    unsigned long sizesample=tbl.nchan*2u;

    // size of 1 minute of data (without auxiliary channels)
    unsigned long sizemin=tbl.srpm*sizesample;

    // first 10 minute mark
    clock.imark=(10-tbl.stime.min%10)*tbl.srpm;

    // DATA BLOCKS SIZE:
    //     b1: size of the first data block before the first 10 minute mark
    unsigned long b1=(clock.imark*tbl.nchan+6)*2u;
    //     bi: size of a full 10 minutes data block
    unsigned long bi=(10*tbl.srpm*tbl.nchan+6)*2u;

    //  calculate duration of data acquisition in minutes(Tmin) and seconds (Tsec)
    unsigned long Tmin,
                  Tsec;
    unsigned datasize=tbl.sizefile-tbl.sizehdr;

    if(datasize<=b1) {    // if data file finish before first 10 minutes mark
        datasize-=6*2u;  // remove size of auxiliary channels
        Tmin=datasize/sizemin;
        datasize%=sizemin;
        if(datasize%sizesample!=0)
            cout<<"WARNING1: some record without "<<tbl.nchan<<" channels recorded"<<endl;
        Tsec=(datasize/sizesample)*long(60./tbl.srpm);
    }
    else {
        Tmin=10-tbl.stime.min%10;  // minutes of first block
        datasize-=b1;              // remove size of first block
        Tmin+=(datasize/bi)*10;    // minutes of full 10 minutes data block
        datasize%=bi;

        if(datasize!=0) { // if rest some data
            datasize-=6*2u;  // remove size of auxiliary channels
            Tmin+=datasize/sizemin;
            datasize%=sizemin;
            if(datasize%sizesample!=0)
                cout<<"WARNING2: some record without "<<tbl.nchan<<" channels recorded"<<endl;
            Tsec=(datasize/sizesample)*long(60./tbl.srpm);
        }
    }


    // endtime
    Tmin+=clock.swin.min;  // add minutes of windows start time

    clock.ewin.jul=clock.swin.jul+Tmin/1440L;
    clock.ewin.min=Tmin%1440L;
    clock.ewin.sec=Tsec;

    // Set drop points
    long idrop;
    if(drop) {
        if(tbl.srpm==12*60)
            idrop=17;
        else if(tbl.srpm==4*60)
            idrop=19;
        else if(tbl.srpm==60)
            idrop=22;
        else if(tbl.srpm==60/5)
            idrop=20;
        else if(tbl.srpm==60/30) {
            idrop=19;
            clock.swin.sec=25;  // but don't enter in count calculus
        }
        else {
            cout<<tbl.srpm<<" sample rate per minute is an unknow sample rate"<<endl
                <<"it could cause problems"<<endl;
            idrop=0;
        }
        
        // correct end time
        // note that to drop the points it will be added to the startcount
        // but as it will be removed from the end time to
        // it must be added again in the endcount
        
        long secs=clock.ewin.min*60+clock.ewin.sec;
        long dropsecs=long(idrop*(60./tbl.srpm));
        if(secs>=dropsecs) {
            secs-=dropsecs;
            clock.ewin.min=secs/60;
            clock.ewin.sec=secs%60;
        }
        else {
            clock.ewin.jul--;
            secs+=86400;
            secs-=dropsecs;
            clock.ewin.min=secs/60;
            clock.ewin.sec=secs%60;
        }
        
    }
    else
        idrop=0;

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

        printf("[%04d-%02d-%02dT%02d:%02d:%02d]\n ",year,month,day,hour,min,clock.swin.sec);
        getline(cin,resp);
        if(resp.size()==0) { // choose default option
            clock.sdata.jul=clock.swin.jul;
            clock.sdata.min=clock.swin.min;
            clock.sdata.sec=clock.swin.sec;
            break;
        }
        else if(parse_date(resp.c_str(),century,50,FMT_nohint,&julian,&recognized)==EXIT_SUCCESS) {
            double dsec;
            jul_to_cal_and_time(julian,0.5,&year,&month,&day,&hour,&min,&dsec);
            clock.sdata.jul=cal_to_jul(year,month,day);
            clock.sdata.min=min+hour*60;
            clock.sdata.sec=long(dsec);

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

        printf("[%04d-%02d-%02dT%02d:%02d:%02d]\n ",year,month,day,hour,min,clock.ewin.sec);
        getline(cin,resp);
        if(resp.size()==0) { // choose default option
            clock.edata.jul=clock.ewin.jul;
            clock.edata.min=clock.ewin.min;
            clock.edata.sec=clock.ewin.sec;
            break;
        }
        else if(parse_date(resp.c_str(),century,50,FMT_nohint,&julian,&recognized)==EXIT_SUCCESS) {
            double dsec;
            jul_to_cal_and_time(julian,0.5,&year,&month,&day,&hour,&min,&dsec);
            clock.edata.jul=cal_to_jul(year,month,day);
            clock.edata.min=min+hour*60;
            clock.edata.sec=long(dsec);

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
    // number of first sample record to read
    clock.scount=(clock.sdata.jul-clock.swin.jul)*1440L;
    clock.scount+=clock.sdata.min-clock.swin.min;
    clock.scount*=tbl.srpm;
    clock.scount+=idrop;

    // number of last sample record to read
    clock.ecount=(clock.edata.jul-clock.swin.jul)*1440L;
    clock.ecount+=clock.edata.min-clock.swin.min;
    clock.ecount*=tbl.srpm;
    clock.ecount+=long(clock.edata.sec*(tbl.srpm/60.));
    clock.ecount+=idrop;

    return clock;
}
