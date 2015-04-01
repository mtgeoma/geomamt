#include "udouble.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

//  This  is  the  default  conversion  from  type  double
UDouble::UDouble(const  long double  val,  const  long double var)
  : value(val), variance(var)
{
  if  (var  <  0.0)
  {
    cerr  <<  "Error:  negative  variance:  "  <<  var <<  endl;
    exit(EXIT_FAILURE);
  }
}

UDouble::UDouble(const  UDouble&  ud)
  :  value(ud.value),  variance(ud.variance)  {}

UDouble::~UDouble(void)  {}

UDouble::UDouble  operator+(UDouble::UDouble const& d)
{
  return  UDouble(d.value, d.variance);
}

UDouble::UDouble  operator-(UDouble::UDouble const& d)
{
  return  UDouble(-d.value,  d.variance);
}

long double  UDouble::mean(void)  const
{
  return value;
}

long double  UDouble::deviation(void)  const
{
  return sqrt(variance);
}

UDouble::UDouble operator+(UDouble const& x, const  UDouble& y)
{
  return UDouble(x.value + y.value,
		 x.variance + y.variance);
}

UDouble::UDouble operator-(UDouble const& x, const  UDouble& y)
{
  return UDouble(x.value - y.value,
		 x.variance + y.variance);
}

UDouble::UDouble operator*(UDouble const& x, const  UDouble& y)
{
  return UDouble(x.value * y.value,
		 pow(y.value, 2) * x.variance
		 + pow(x.value, 2) * y.variance);
}

UDouble::UDouble operator/(UDouble const& x, const  UDouble& y)
{
  return UDouble(x.value / y.value,
		 pow(1.0 / y.value, 2) * x.variance
		 + pow(x.value / pow(y.value, 2), 2) * y.variance);
}

UDouble::UDouble operator/(UDouble const& x, const long double y)
{
  return UDouble(x.value / y,
		 pow(1.0 / y, 2) * x.variance);
}

UDouble::UDouble  pow(const  UDouble& x, const  long double y)
{
  return UDouble(pow(x.value, y),
		 pow(y * pow(x.value, y - 1.0), 2) * x.variance);
}

UDouble::UDouble  sqrt(const  UDouble& x)
{
  return UDouble( sqrt(x.value),
		  x.variance / (4 * abs( x.value ) ) );
}

UDouble::UDouble  abs(UDouble::UDouble const& d)
{
  return  UDouble(abs(d.value),  d.variance);
}

UDouble::UDouble&
UDouble::operator+=(
		    UDouble const& ud)
{
  variance += ud.variance;
  value  +=  ud.value;
  return  *this;
}

std::ostream&
operator<< (std::ostream& o, const UDouble& ud)
{
  return o << ud.value << "(" << sqrt(ud.variance) << ")";
} // operator<<

//End of File
