#ifndef DECLINATION_H
#define DECLINATION_H

#include <string>
// latitude e longitude em graus
// elev em quilômetros
// year em yyyy
double declination(double latitude, double longitude, double elev, int year, int month, int day, const std::string& model_file);

#endif // DECLINATION_H
