#include "wal_invariants.hpp"
#include "invariants_parameters.hpp"

namespace ip = erro_manual::base_xi_eta::invariant_parameter;
namespace tf = transfer_function;

namespace erro_manual
{
namespace base_xi_eta
{
namespace wal_invariant
{

tf::real_type
I1_value(tf::impedance const& imp)
{
  tf::real_type xi1 = ip::xi1_value(imp);
  tf::real_type xi4 = ip::xi4_value(imp);

  return sqrt( pow(xi1, 2) + pow(xi4, 2) );
} // I1_value

tf::real_type
I1_variance(tf::impedance const& imp)
{
  tf::real_type xi1     = ip::xi1_value(imp);
  tf::real_type xi1_var = ip::xi1_variance(imp);
  tf::real_type xi4     = ip::xi4_value(imp);
  tf::real_type xi4_var = ip::xi4_variance(imp);

  tf::real_type d_xi1 =
    xi1/sqrt(pow(xi4,2) + pow(xi1,2));
  tf::real_type d_xi4 =
    xi4/sqrt(pow(xi4,2) + pow(xi1,2));

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi4, 2 ) * xi4_var;
} // I1_variance

tf::real_measurement_type
I1(tf::impedance const& imp)
{
  return tf::real_measurement_type( I1_value(imp), I1_variance(imp) );
}//I1

tf::real_type
I2_value(tf::impedance const& imp)
{
  tf::real_type eta1 = ip::eta1_value(imp);
  tf::real_type eta4 = ip::eta4_value(imp);

  return sqrt( pow(eta1, 2) + pow(eta4, 2) );
} // I2_value

tf::real_type
I2_variance(tf::impedance const& imp)
{
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_eta1 =
    eta1/sqrt(pow(eta4,2) + pow(eta1,2));
  tf::real_type d_eta4 =
    eta4/sqrt(pow(eta4,2) + pow(eta1,2));

  return
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta4, 2 ) * eta4_var;
} // I2_variance

tf::real_measurement_type
I2(tf::impedance const& imp)
{
  return tf::real_measurement_type( I2_value(imp), I2_variance(imp) );
}//I2

tf::real_type
I3_value(tf::impedance const& imp)
{
  tf::real_type xi1 = ip::xi1_value(imp);
  tf::real_type xi2 = ip::xi2_value(imp);
  tf::real_type xi3 = ip::xi3_value(imp);
  tf::real_type xi4 = ip::xi4_value(imp);

  return sqrt( ( pow(xi2, 2) + pow(xi3, 2) )
	       /
	       ( pow(xi1, 2) + pow(xi4, 2) ) );
} // I3_value

tf::real_type
I3_variance(tf::impedance const& imp)
{
  tf::real_type xi1 = ip::xi1_value(imp);
  tf::real_type xi2 = ip::xi2_value(imp);
  tf::real_type xi3 = ip::xi3_value(imp);
  tf::real_type xi4 = ip::xi4_value(imp);
  tf::real_type xi1_var = ip::xi1_variance(imp);
  tf::real_type xi2_var = ip::xi2_variance(imp);
  tf::real_type xi3_var = ip::xi3_variance(imp);
  tf::real_type xi4_var = ip::xi4_variance(imp);

  tf::real_type d_xi1 =
    -xi1*sqrt(pow(xi3,2) + pow(xi2,2))/pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0);
  tf::real_type d_xi2 =
    xi2/(sqrt(pow(xi3,2) + pow(xi2,2))*sqrt(pow(xi4,2) + pow(xi1,2)));
  tf::real_type d_xi3 =
    xi3/(sqrt(pow(xi3,2) + pow(xi2,2))*sqrt(pow(xi4,2) + pow(xi1,2)));
  tf::real_type d_xi4 =
    -sqrt(pow(xi3,2) + pow(xi2,2))*xi4/pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0);

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi2, 2 ) * xi2_var
    +
    pow( d_xi3, 2 ) * xi3_var + pow( d_xi4, 2 ) * xi4_var;
} // I3_variance

tf::real_measurement_type
I3(tf::impedance const& imp)
{
  return tf::real_measurement_type( I3_value(imp), I3_variance(imp) );
}//I3

tf::real_type
I4_value(tf::impedance const& imp)
{
  tf::real_type eta1 = ip::eta1_value(imp);
  tf::real_type eta2 = ip::eta2_value(imp);
  tf::real_type eta3 = ip::eta3_value(imp);
  tf::real_type eta4 = ip::eta4_value(imp);

  return sqrt( ( pow(eta2, 2) + pow(eta3, 2) )
	       /
	       ( pow(eta1, 2) + pow(eta4, 2) ) );
} // I4_value

tf::real_type
I4_variance(tf::impedance const& imp)
{
  tf::real_type eta1 = ip::eta1_value(imp);
  tf::real_type eta2 = ip::eta2_value(imp);
  tf::real_type eta3 = ip::eta3_value(imp);
  tf::real_type eta4 = ip::eta4_value(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta2_var = ip::eta2_variance(imp);
  tf::real_type eta3_var = ip::eta3_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_eta1 =
    -eta1*sqrt(pow(eta3,2) + pow(eta2,2))/pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0);
  tf::real_type d_eta2 =
    eta2/(sqrt(pow(eta3,2) + pow(eta2,2))*sqrt(pow(eta4,2) + pow(eta1,2)));
  tf::real_type d_eta3 =
    eta3/(sqrt(pow(eta3,2) + pow(eta2,2))*sqrt(pow(eta4,2) + pow(eta1,2)));
  tf::real_type d_eta4 =
    -sqrt(pow(eta3,2) + pow(eta2,2))*eta4/pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0);

  return
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta2, 2 ) * eta2_var
    +
    pow( d_eta3, 2 ) * eta3_var + pow( d_eta4, 2 ) * eta4_var;
} // I4_variance

tf::real_measurement_type
I4(tf::impedance const& imp)
{
  return tf::real_measurement_type( I4_value(imp), I4_variance(imp) );
}//I4

tf::real_type
I5_value(tf::impedance const& imp)
{
  tf::real_type xi1   = ip::xi1_value(imp);
  tf::real_type xi4   = ip::xi4_value(imp);
  tf::real_type eta1  = ip::eta1_value(imp);
  tf::real_type eta4  = ip::eta4_value(imp);
  tf::real_type I1    = I1_value(imp);
  tf::real_type I2    = I2_value(imp);

  return ( xi4 * eta1 + xi1 * eta4 ) / (I1 * I2);
} // I5_value

tf::real_type
I5_variance(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type xi1_var  = ip::xi1_variance(imp);
  tf::real_type xi4_var  = ip::xi4_variance(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_xi1 =
    eta4/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - xi1*(eta1*xi4 + eta4*xi1)/(sqrt(pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0));
  tf::real_type d_xi4 =
    eta1/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - xi4*(eta1*xi4 + eta4*xi1)/(sqrt(pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0));
  tf::real_type d_eta1 =
    xi4/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - eta1*(eta1*xi4 + eta4*xi1)/(pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0)*sqrt(pow(xi4,2) + pow(xi1,2)));
  tf::real_type d_eta4 =
    xi1/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - eta4*(eta1*xi4 + eta4*xi1)/(pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0)*sqrt(pow(xi4,2) + pow(xi1,2)));

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi4, 2 ) * xi4_var
    +
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta4, 2) * eta4_var;
} // I5_variance

tf::real_measurement_type
I5(tf::impedance const& imp)
{
  return tf::real_measurement_type( I5_value(imp), I5_variance(imp) );
}//I5

tf::real_type
I6_value(tf::impedance const& imp)
{
  tf::real_type xi1   = ip::xi1_value(imp);
  tf::real_type xi4   = ip::xi4_value(imp);
  tf::real_type eta1  = ip::eta1_value(imp);
  tf::real_type eta4  = ip::eta4_value(imp);
  tf::real_type I1    = I1_value(imp);
  tf::real_type I2    = I2_value(imp);

  return ( xi4 * eta1 - xi1 * eta4 ) / (I1 * I2);
} // I6_value

tf::real_type
I6_variance(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type xi1_var  = ip::xi1_variance(imp);
  tf::real_type xi4_var  = ip::xi4_variance(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_xi1 =
    -eta4/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - xi1*(eta1*xi4 - eta4*xi1)/(sqrt(pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0));
  tf::real_type d_xi4 =
    eta1/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - xi4*(eta1*xi4 - eta4*xi1)/(sqrt(pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0));
  tf::real_type d_eta1 =
    xi4/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2))) - eta1*(eta1*xi4 - eta4*xi1)/(pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0)*sqrt(pow(xi4,2) + pow(xi1,2)));
  tf::real_type d_eta4 =
    -eta4*(eta1*xi4 - eta4*xi1)/(pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0)*sqrt(pow(xi4,2) + pow(xi1,2))) - xi1/(sqrt(pow(eta4,2) + pow(eta1,2))*sqrt(pow(xi4,2) + pow(xi1,2)));

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi4, 2 ) * xi4_var
    +
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta4, 2) * eta4_var;

} // I6_variance

tf::real_measurement_type
I6(tf::impedance const& imp)
{
  return tf::real_measurement_type( I6_value(imp), I6_variance(imp) );
}//I6

tf::real_type
I7_value(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);

  return
    (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2)));
} // I7_value

tf::real_type
I7_variance(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type xi1_var  = ip::xi1_variance(imp);
  tf::real_type xi2_var  = ip::xi2_variance(imp);
  tf::real_type xi3_var  = ip::xi3_variance(imp);
  tf::real_type xi4_var  = ip::xi4_variance(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta2_var = ip::eta2_variance(imp);
  tf::real_type eta3_var = ip::eta3_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_xi1 =
    -eta4/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (2*(pow(eta3,2) + pow(eta2,2))*xi1 - 2*(eta1*(eta3*xi3 + eta2*xi2) - eta4*(eta3*xi2 - eta2*xi3)))*(eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_xi2 =
    -eta3/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*(pow(eta4,2) + pow(eta1,2))*xi2 - 2*(eta2*(eta4*xi4 + eta1*xi1) + eta3*(eta1*xi4 - eta4*xi1)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_xi3 =
    eta2/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*(pow(eta4,2) + pow(eta1,2))*xi3 - 2*(eta3*(eta4*xi4 + eta1*xi1) - eta2*(eta1*xi4 - eta4*xi1)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_xi4 =
    eta1/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*(pow(eta3,2) + pow(eta2,2))*xi4 - 2*(eta4*(eta3*xi3 + eta2*xi2) + eta1*(eta3*xi2 - eta2*xi3)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_eta1 =
    xi4/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*eta1*(pow(xi3,2) + pow(xi2,2)) - 2*((eta3*xi2 - eta2*xi3)*xi4 + xi1*(eta3*xi3 + eta2*xi2)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_eta2 =
    xi3/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*eta2*(pow(xi4,2) + pow(xi1,2)) - 2*(xi2*(eta4*xi4 + eta1*xi1) - xi3*(eta1*xi4 - eta4*xi1)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_eta3 =
    -xi2/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*eta3*(pow(xi4,2) + pow(xi1,2)) - 2*(xi3*(eta4*xi4 + eta1*xi1) + xi2*(eta1*xi4 - eta4*xi1)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));
  tf::real_type d_eta4 =
    -xi1/sqrt((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))) - (eta1*xi4 + eta2*xi3 - eta3*xi2 - eta4*xi1)*(2*eta4*(pow(xi3,2) + pow(xi2,2)) - 2*((eta3*xi3 + eta2*xi2)*xi4 - xi1*(eta3*xi2 - eta2*xi3)))/(2*pow(((pow(eta3,2) + pow(eta2,2))*(pow(xi4,2) + pow(xi1,2)) - 2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1)) + (pow(eta4,2) + pow(eta1,2))*(pow(xi3,2) + pow(xi2,2))),3.0/2.0));

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi2, 2 ) * xi2_var
    +
    pow( d_xi3, 2 ) * xi3_var + pow( d_xi4, 2 ) * xi4_var
    +
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta2, 2) * eta2_var
    +
    pow( d_eta3, 2 ) * eta3_var + pow( d_eta4, 2) * eta4_var;

} // I7_variance

tf::real_measurement_type
I7(tf::impedance const& imp)
{
  return tf::real_measurement_type( I7_value(imp), I7_variance(imp) );
}//I7

tf::real_type
Q_value(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);

  return
    sqrt(-2.0*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2)));
} // Q_value

tf::real_type
Q_variance(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type xi1_var  = ip::xi1_variance(imp);
  tf::real_type xi2_var  = ip::xi2_variance(imp);
  tf::real_type xi3_var  = ip::xi3_variance(imp);
  tf::real_type xi4_var  = ip::xi4_variance(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta2_var = ip::eta2_variance(imp);
  tf::real_type eta3_var = ip::eta3_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_xi1 =
    (-2*(eta1*(eta3*xi3 + eta2*xi2) - eta4*(eta3*xi2 - eta2*xi3))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + 4*xi1*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),2)) - 2*xi1*(pow(xi3,2) + pow(xi2,2))/pow((pow(xi4,2) + pow(xi1,2)),2))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_xi2 =
    (2*xi2/(pow(xi4,2) + pow(xi1,2)) - 2*(eta2*(eta4*xi4 + eta1*xi1) + eta3*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_xi3 =
    (2*xi3/(pow(xi4,2) + pow(xi1,2)) - 2*(eta3*(eta4*xi4 + eta1*xi1) - eta2*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_xi4 =
    (-2*(eta4*(eta3*xi3 + eta2*xi2) + eta1*(eta3*xi2 - eta2*xi3))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + 4*xi4*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*pow((pow(xi4,2) + pow(xi1,2)),2)) - 2*(pow(xi3,2) + pow(xi2,2))*xi4/pow((pow(xi4,2) + pow(xi1,2)),2))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_eta1 =
    (4*eta1*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/(pow((pow(eta4,2) + pow(eta1,2)),2)*(pow(xi4,2) + pow(xi1,2))) - 2*((eta3*xi2 - eta2*xi3)*xi4 + xi1*(eta3*xi3 + eta2*xi2))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) - 2*eta1*(pow(eta3,2) + pow(eta2,2))/pow((pow(eta4,2) + pow(eta1,2)),2))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_eta2 =
    (2*eta2/(pow(eta4,2) + pow(eta1,2)) - 2*(xi2*(eta4*xi4 + eta1*xi1) - xi3*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_eta3 =
    (2*eta3/(pow(eta4,2) + pow(eta1,2)) - 2*(xi3*(eta4*xi4 + eta1*xi1) + xi2*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));
  tf::real_type d_eta4 =
    (4*eta4*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/(pow((pow(eta4,2) + pow(eta1,2)),2)*(pow(xi4,2) + pow(xi1,2))) - 2*((eta3*xi3 + eta2*xi2)*xi4 - xi1*(eta3*xi2 - eta2*xi3))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) - 2*(pow(eta3,2) + pow(eta2,2))*eta4/pow((pow(eta4,2) + pow(eta1,2)),2))/(2*sqrt(-2*((eta3*xi3 + eta2*xi2)*(eta4*xi4 + eta1*xi1) + (eta3*xi2 - eta2*xi3)*(eta1*xi4 - eta4*xi1))/((pow(eta4,2) + pow(eta1,2))*(pow(xi4,2) + pow(xi1,2))) + (pow(xi3,2) + pow(xi2,2))/(pow(xi4,2) + pow(xi1,2)) + (pow(eta3,2) + pow(eta2,2))/(pow(eta4,2) + pow(eta1,2))));

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi2, 2 ) * xi2_var
    +
    pow( d_xi3, 2 ) * xi3_var + pow( d_xi4, 2 ) * xi4_var
    +
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta2, 2) * eta2_var
    +
    pow( d_eta3, 2 ) * eta3_var + pow( d_eta4, 2) * eta4_var;
} // Q_variance

tf::real_measurement_type
Q(tf::impedance const& imp)
{
  return tf::real_measurement_type( Q_value(imp), Q_variance(imp) );
}//Q

tf::real_type
cos_alpha_value(tf::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);

  return
    xi4/sqrt(pow(xi4,2) + pow(xi1,2));
} // cos_alpha_value

transfer_function::real_type
cos_alpha_variance(transfer_function::impedance const& imp)
{
  tf::real_type xi1     = ip::xi1_value(imp);
  tf::real_type xi1_var = ip::xi1_variance(imp);
  tf::real_type xi4     = ip::xi4_value(imp);
  tf::real_type xi4_var = ip::xi4_variance(imp);

  tf::real_type d_xi1 =
    -xi1*xi4/pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0);
  tf::real_type d_xi4 =
    1/sqrt(pow(xi4,2) + pow(xi1,2)) - pow(xi4,2)/pow((pow(xi4,2) + pow(xi1,2)),3.0/2.0);

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi4, 2 ) * xi4_var;
} // cos_alpha_variance

transfer_function::real_measurement_type
cos_alpha(transfer_function::impedance const& imp)
{
  return tf::real_measurement_type( cos_alpha_value(imp), cos_alpha_variance(imp) );
}//cos_alpha

transfer_function::real_type
cos_beta_value(transfer_function::impedance const& imp)
{
  tf::real_type eta1      = ip::eta1_value(imp);
  tf::real_type eta4      = ip::eta4_value(imp);

  return
    eta4/sqrt(pow(eta4,2) + pow(eta1,2));
} // cos_beta_value

transfer_function::real_type
cos_beta_variance(transfer_function::impedance const& imp)
{
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_eta1 =
    -eta1*eta4/pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0);
  tf::real_type d_eta4 =
    1/sqrt(pow(eta4,2) + pow(eta1,2)) - pow(eta4,2)/pow((pow(eta4,2) + pow(eta1,2)),3.0/2.0);

  return
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta4, 2 ) * eta4_var;
} // cos_beta_variance

transfer_function::real_measurement_type
cos_beta(transfer_function::impedance const& imp)
{
  return tf::real_measurement_type( cos_beta_value(imp), cos_beta_variance(imp) );
}//cos_beta

transfer_function::real_type
N_value(transfer_function::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);

  //     return sqrt(
  //                 sqrt( pow( (eta1*eta1+xi2*xi2+xi3*xi3)- (xi1*xi1+eta2*eta2+eta3*eta3) ,2) + 4.0*pow(xi2*eta2+xi3*eta3-xi1*eta1,2) )
  //                 /
  //                 (xi4*xi4+eta4*eta4) );
  return pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),1.0/4.0)/sqrt(pow(xi4,2) + pow(eta4,2));

} // N_value

transfer_function::real_type
N_variance(transfer_function::impedance const& imp)
{
  tf::real_type xi1      = ip::xi1_value(imp);
  tf::real_type xi2      = ip::xi2_value(imp);
  tf::real_type xi3      = ip::xi3_value(imp);
  tf::real_type xi4      = ip::xi4_value(imp);
  tf::real_type eta1     = ip::eta1_value(imp);
  tf::real_type eta2     = ip::eta2_value(imp);
  tf::real_type eta3     = ip::eta3_value(imp);
  tf::real_type eta4     = ip::eta4_value(imp);
  tf::real_type xi1_var  = ip::xi1_variance(imp);
  tf::real_type xi2_var  = ip::xi2_variance(imp);
  tf::real_type xi3_var  = ip::xi3_variance(imp);
  tf::real_type xi4_var  = ip::xi4_variance(imp);
  tf::real_type eta1_var = ip::eta1_variance(imp);
  tf::real_type eta2_var = ip::eta2_variance(imp);
  tf::real_type eta3_var = ip::eta3_variance(imp);
  tf::real_type eta4_var = ip::eta4_variance(imp);

  tf::real_type d_xi1 =
    (-4*xi1*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)) - 8.0*eta1*(eta3*xi3 + eta2*xi2 - eta1*xi1))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_xi2 =
    (4*xi2*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)) + 8.0*eta2*(eta3*xi3 + eta2*xi2 - eta1*xi1))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_xi3 =
    (4*xi3*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)) + 8.0*eta3*(eta3*xi3 + eta2*xi2 - eta1*xi1))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_xi4 =
    -pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),1.0/4.0)*xi4/pow((pow(xi4,2) + pow(eta4,2)),3.0/2.0);
  tf::real_type d_eta1 =
    (4*eta1*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)) - 8.0*xi1*(eta3*xi3 + eta2*xi2 - eta1*xi1))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_eta2 =
    (8.0*xi2*(eta3*xi3 + eta2*xi2 - eta1*xi1) - 4*eta2*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_eta3 =
    (8.0*xi3*(eta3*xi3 + eta2*xi2 - eta1*xi1) - 4*eta3*(pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)))/(4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),3.0/4.0)*sqrt(pow(xi4,2) + pow(eta4,2)));
  tf::real_type d_eta4 =
    -eta4*pow((pow((pow(xi3,2) + pow(xi2,2) - pow(xi1,2) - pow(eta3,2) - pow(eta2,2) + pow(eta1,2)),2) + 4.0*pow((eta3*xi3 + eta2*xi2 - eta1*xi1),2)),1.0/4.0)/pow((pow(xi4,2) + pow(eta4,2)),3.0/2.0);

  return
    pow( d_xi1, 2 ) * xi1_var + pow( d_xi2, 2 ) * xi2_var
    +
    pow( d_xi3, 2 ) * xi3_var + pow( d_xi4, 2 ) * xi4_var
    +
    pow( d_eta1, 2 ) * eta1_var + pow( d_eta2, 2) * eta2_var
    +
    pow( d_eta3, 2 ) * eta3_var + pow( d_eta4, 2) * eta4_var;
} // N_variance

transfer_function::real_measurement_type
N(transfer_function::impedance const& imp)
{
  return tf::real_measurement_type( N_value(imp), N_variance(imp) );
}//N

} //namespace wal_invariant

} //namespace base_xi_eta

} //namespace erro_manual
