#ifndef E30954AB95D40A6A0B20276B5357BEC6
#define E30954AB95D40A6A0B20276B5357BEC6

float julday(int i_month, int i_day, int i_year);
int calc_geomagXYZ(const char* mdfile,float latitude, float longitude,float alt, float sdate,
                   float *X, float *Y, float *Z);

#endif // MY_GEOMAG60_H
