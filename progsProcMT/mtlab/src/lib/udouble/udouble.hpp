#ifndef UDOUBLE_HPP
#define UDOUBLE_HPP

#include <ostream>
#include <istream>
#include <iostream>

//  model uncertain number using only mean and sigma (pure Gaussian)
class UDouble
{
private:
  long double value;
  long double variance;
public:
  UDouble(const long double val = 0.0, const long double var = 0.0);
  UDouble(const UDouble&);
  ~UDouble(void);
  friend UDouble operator+(UDouble const& d);
  friend UDouble operator-(UDouble const& d);
  friend UDouble operator+(UDouble const& x, UDouble const& y);
  friend UDouble operator-(UDouble const& x, UDouble const& y);
  friend UDouble operator*(UDouble const& x, UDouble const& y);
  friend UDouble operator/(UDouble const& x, UDouble const& y);
  friend UDouble operator/(UDouble const& x, const long double y);
  friend std::ostream& operator<< (std::ostream& o, const UDouble& ud);

  UDouble& operator+=(UDouble const& ud);
  UDouble& operator-=(UDouble const&);
  UDouble& operator*=(UDouble const&);
  UDouble& operator/=(UDouble const&);
  //  math library functions
  friend UDouble ceil(UDouble);
  friend UDouble floor(UDouble);
  friend UDouble ldexp(UDouble, int);
  friend UDouble modf(UDouble, long double *);
  friend UDouble frexp(UDouble, int *);
  friend UDouble fmod(UDouble const&, UDouble const&);
  friend UDouble sqrt(UDouble const& x);
  friend UDouble abs(UDouble const& x);
  friend UDouble sin(UDouble);
  friend UDouble cos(UDouble);
  friend UDouble tan(UDouble);
  friend UDouble asin(UDouble);
  friend UDouble acos(UDouble);
  friend UDouble atan(UDouble);
  friend UDouble atan2(UDouble const&, UDouble const&);
  friend UDouble exp(UDouble);
  friend UDouble log(UDouble);
  friend UDouble log10(UDouble);
  friend UDouble sinh(UDouble);
  friend UDouble cosh(UDouble);
  friend UDouble tanh(UDouble);
  friend UDouble pow(UDouble const& x, UDouble const& y);
  friend UDouble pow(UDouble const& x, const long double y);

  //  read-only access to data members
  long double mean(void) const;
  long double deviation(void) const;
  friend UDouble PropagateUncertaintiesBySlope(long double (*)(long double),
					       UDouble const&);
  friend UDouble PropagateUncertaintiesBySlope(long double (*)(long double,
							       long double),
					       UDouble const&,
					       UDouble const&);
};


/* End of File */

#endif
