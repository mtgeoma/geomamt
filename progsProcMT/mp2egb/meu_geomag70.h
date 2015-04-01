#ifndef MY_GEOMAG70_H
#define MY_GEOMAG70_H

double julday(int month, int day, int year);
double get_declination(const char* mdfile,double latitude, double longitude,double alt, double sdate);

#endif // MY_GEOMAG70_H
