#include "transfer_functions.hpp"
#include <cmath>

void
transfer_function::impedance::rotate(real_type angle)
{
  if ( angle == 0.0 || angle / 360.0 == 1.0) return;

  angle = M_PIl * angle / 180.0L; //graus -> radianos
  complex_type Zxx = complex_type(Z(0, 0).real().value, Z(0, 0).imag().value);
  complex_type Zxy = complex_type(Z(0, 1).real().value, Z(0, 1).imag().value);
  complex_type Zyx = complex_type(Z(1, 0).real().value, Z(1, 0).imag().value);
  complex_type Zyy = complex_type(Z(1, 1).real().value, Z(1, 1).imag().value);

  real_type Zxx_v = Z(0, 0).real().variance;
  real_type Zxy_v = Z(0, 1).real().variance;
  real_type Zyx_v = Z(1, 0).real().variance;
  real_type Zyy_v = Z(1, 1).real().variance;

  complex_type Zxx_ =
    ( (Zxx + Zyy) + (Zxy + Zyx)*sinl(2.0 * angle)
      +
      (Zxx - Zyy)*cosl(2.0 * angle) ) / 2.0L;

  complex_type Zxy_ =
    ( (Zxy - Zyx) - (Zxx - Zyy)*sinl(2.0 * angle)
      +
      (Zxy + Zyx) * cosl(2 * angle) ) / 2.0L;

  complex_type Zyx_ =
    ( -(Zxy - Zyx) - (Zxx - Zyy)*sinl(2.0 * angle)
      +
      (Zxy + Zyx) * cosl(2.0 * angle) ) / 2.0L;

  complex_type Zyy_ =
    ( (Zxx + Zyy) - (Zxy + Zyx)*sinl(2.0 * angle)
      -
      (Zxx - Zyy) * cosl(2.0 * angle) ) / 2.0L;

  real_type Zxx_v_ =
    pow(1.0L + cosl(2.0 * angle) , 2) * Zxx_v
    +
    pow(sinl(2.0 * angle), 2) * Zxy_v
    +
    pow(sinl(2.0 * angle), 2) * Zyx_v
    +
    pow(1.0L - cosl(2.0 * angle), 2) * Zyy_v;

  real_type Zxy_v_ =
    pow(-sinl(2.0 * angle), 2) * Zxx_v
    +
    pow(1.0L + cosl(2.0 * angle), 2) * Zxy_v
    +
    pow(-1.0L + cosl(2.0 * angle), 2) * Zyx_v
    +
    pow(sinl(2.0 * angle), 2) * Zyy_v;

  real_type Zyx_v_ =
    pow(-sinl(2.0 * angle), 2) * Zxx_v
    +
    pow(-1.0L + cosl(2.0 * angle), 2) * Zxy_v
    +
    pow(+1.0L + cosl(2.0 * angle), 2) * Zyx_v
    +
    pow(sinl(2.0 * angle), 2) * Zyy_v;

  real_type Zyy_v_ =
    pow(1.0L - cosl(2.0 * angle), 2) * Zxx_v
    +
    pow(-sinl(2.0 * angle), 2) * Zxy_v
    +
    pow(-sinl(2.0 * angle), 2) * Zyx_v
    +
    pow(1.0L + cosl(2.0 * angle), 2) * Zyy_v;

  Z(0, 0) = complex_measurement_type( real_measurement_type(Zxx_.real(),
							    Zxx_v_),
				      real_measurement_type(Zxx_.imag(),
							    Zxx_v_) );

  Z(0, 1) = complex_measurement_type( real_measurement_type(Zxy_.real(),
							    Zxy_v_),
				      real_measurement_type(Zxy_.imag(),
							    Zxy_v_) );

  Z(1, 0) = complex_measurement_type( real_measurement_type(Zyx_.real(),
							    Zyx_v_),
				      real_measurement_type(Zyx_.imag(),
							    Zyx_v_) );

  Z(1, 1) = complex_measurement_type( real_measurement_type(Zyy_.real(),
							    Zyy_v_),
				      real_measurement_type(Zyy_.imag(),
							    Zyy_v_) );
}
