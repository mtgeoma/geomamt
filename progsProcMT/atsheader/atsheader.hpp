#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>

void read_header (std::fstream &ats, std::string parameter);
void read_header (std::fstream &ats);
void write_header (std::fstream &ats, std::string parameter, std::string value);
