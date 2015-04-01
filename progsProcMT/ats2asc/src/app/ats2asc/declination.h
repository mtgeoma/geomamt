#ifndef AC330A3E6429CB54A593FCA330F97738
#define AC330A3E6429CB54A593FCA330F97738

#include <string>
// latitude e longitude em graus
// elev em quilômetros
// year em yyyy
double declination(double latitude, double longitude, double elev, int year,
		   int month, int day, const std::string& model_file);

#endif // DECLINATION_H
