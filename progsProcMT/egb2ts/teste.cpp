#include "boost/date_time/posix_time/posix_time.hpp"
#include <iostream>

int main()
{
  using namespace boost::gregorian;
  using namespace boost::posix_time;
  using namespace std;
  double sr;
  int year, month, day, hour, min, sec, micro;
  year=2009;
  month=7;
  day=24;
  hour=12;
  min=2;
  sr=2.4414062500e-05;
  sec=int(sr);
  micro=int((sr-sec)*1.0e6);
  date data_inicial(year,month,day);
  ptime smp_time(data_inicial,hours(hour)+minutes(min)+seconds(sec)+microseconds(micro));
  cout << to_iso_extended_string(smp_time) << "\n";
  smp_time+=seconds(sec)+microseconds(micro);
  cout << to_iso_extended_string(smp_time) << "\n";
}
