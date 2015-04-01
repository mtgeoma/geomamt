#include "dcmp2gmt.h"

vector<datum> tofirst (vector<datum> data) {
    for(int i=0; i<data.size();i++) {
        if(data[i].val>90. && data[i].val<=270.) {
            data[i].val-=180.;
            data[i].min-=180.;
            data[i].max-=180.;
        }
        else if(data[i].val>-270. && data[i].val<=-90.) {
            data[i].val+=180.;
            data[i].min+=180.;
            data[i].max+=180.;
        }
        else if(data[i].val>=-360. && data[i].val<=-270.) {
            data[i].val+=360.;
            data[i].min+=360.;
            data[i].max+=360.;
        }
        else if(data[i].val>270. && data[i].val<=360.) {
            data[i].val-=360.;
            data[i].min-=360.;
            data[i].max-=360.;
        }
    }

    return data;
}
