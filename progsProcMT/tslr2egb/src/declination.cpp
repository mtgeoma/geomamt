#include <cmath>
#include "declination.h"
#include "meu_geomag60.h"
float declination(float latitude, float longitude, float elev, int year, int month, int day) {
    float rad2deg=180./M_PI;

    float sdate = julday(month,day,year);

    const char* modelfile="/usr/local/share/IGRF10.cof";
    float X,Y,Z;
    calc_geomagXYZ(modelfile,latitude,longitude,elev,sdate, &X, &Y, &Z);
    float H = sqrt(X*X+Y*Y);
    float D =  rad2deg*2.*atan2(Y,H+X); // tan(D/2)=tan(D)/(1+sec(D)); tan(D)=Y/X; sec(D)=H/X
    return D;
}
