#include <iostream>
#include "declination.h"

extern "C"
{
#include "meu_geomag70.h"
}

double declination(double latitude, double longitude, double elev, int year,
		   int month, int day, const std::string& model_file)
{
  if((latitude<-90)||(latitude>90))
  {
    std::cout << "latitude de " << latitude
	      << " está fora de escala\ndeclinação será nula\n";
    return 0.0;
  }
  if((longitude<-180)||(longitude>180))
  {
    std::cout << "longitude de " << longitude
	      << " está fora de escala\ndeclinação será nula\n";
    return 0.0;
  }

  double sdate = julday(month,day,year);

  return get_declination(model_file.c_str(),latitude,longitude,elev,sdate);
}
