#include "bahr_invariants.hpp"
#include "bahr_invariants_impl.hpp"
#include "invariants_parameters.hpp"

namespace ip = erro_manual::base_xi_eta::invariant_parameter;
namespace tf = transfer_function;

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::kappa_squared( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi1  = ip::xi1_value(imp);
    tf::real_type xi4  = ip::xi4_value(imp);
    tf::real_type eta1 = ip::eta1_value(imp);
    tf::real_type eta4 = ip::eta4_value(imp);

    return ( pow(xi1, 2) + pow(eta1, 2) ) / ( pow(xi4, 2) + pow(eta4, 2) );
}// kappa_squared

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::kappa_value( 
    transfer_function::impedance const& imp)
{
    return sqrt( kappa_squared(imp) );
}// kappa_value

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::kappa_variance( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi1       = ip::xi1_value(imp);
    tf::real_type xi1_var   = ip::xi1_variance(imp);
    tf::real_type xi4       = ip::xi4_value(imp);
    tf::real_type xi4_var   = ip::xi4_variance(imp);
    tf::real_type eta1      = ip::eta1_value(imp);
    tf::real_type eta1_var  = ip::eta1_variance(imp);
    tf::real_type eta4      = ip::eta4_value(imp);
    tf::real_type eta4_var  = ip::eta4_variance(imp);
    tf::real_type kappa_sq  = kappa_squared(imp);

    return
	( pow(xi1, 2) * xi1_var + pow(eta1, 2) * eta1_var
	  +
	  kappa_sq * ( pow(xi4, 2) * xi4_var + pow(eta4, 2) * eta4_var ) )
	/
	( kappa_sq * pow( pow(xi4, 2) + pow(eta4, 2), 2 ) );
}// kappa_variance

transfer_function::real_measurement_type
erro_manual::base_xi_eta::bahr_invariant::kappa(
    transfer_function::impedance const& imp)
{
    return tf::real_measurement_type( kappa_value(imp), kappa_variance(imp) );
}//kappa

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::mu_squared( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi1  = ip::xi1_value(imp);
    tf::real_type xi2  = ip::xi2_value(imp);
    tf::real_type xi3  = ip::xi3_value(imp);
    tf::real_type xi4  = ip::xi4_value(imp);
    tf::real_type eta1 = ip::eta1_value(imp);
    tf::real_type eta2 = ip::eta2_value(imp);
    tf::real_type eta3 = ip::eta3_value(imp);
    tf::real_type eta4 = ip::eta4_value(imp);

    return 
	( fabs(xi3 * eta2 - xi2 * eta3) + fabs( xi1 * eta4 - xi4 * eta1 ) )
	/
	( pow(xi4, 2) + pow(eta4, 2) );
}// mu_squared

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::mu_value( 
    transfer_function::impedance const& imp)
{
    return sqrt( mu_squared(imp) );
}// mu_value

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::mu_variance( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi4      = ip::xi4_value(imp);
    tf::real_type xi4_var  = ip::xi4_variance(imp);
    tf::real_type eta4     = ip::eta4_value(imp);
    tf::real_type eta4_var = ip::eta4_variance(imp);
    tf::real_type mu_sq    = mu_squared(imp);
    tf::real_type d32_var  = d_board_variance(3, 2, imp);
    tf::real_type d14_var  = d_board_variance(1, 4, imp);
    return 
	( d32_var + d14_var 
	  +
	  4 * mu_sq * ( pow(xi4, 2) * xi4_var + pow(eta4, 2) * eta4_var ) )
	/
	( 4 * mu_sq * pow( pow(xi4, 2) + pow(eta4, 2), 2 ) );
}// mu_variance

transfer_function::real_measurement_type
erro_manual::base_xi_eta::bahr_invariant::mu(
    transfer_function::impedance const& imp)
{
    return tf::real_measurement_type( mu_value(imp), mu_variance(imp) );
}//mu

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::eta_squared( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi1  = ip::xi1_value(imp);
    tf::real_type xi2  = ip::xi2_value(imp);
    tf::real_type xi3  = ip::xi3_value(imp);
    tf::real_type xi4  = ip::xi4_value(imp);
    tf::real_type eta1 = ip::eta1_value(imp);
    tf::real_type eta2 = ip::eta2_value(imp);
    tf::real_type eta3 = ip::eta3_value(imp);
    tf::real_type eta4 = ip::eta4_value(imp);

    return 
	fabs(xi3 * eta2 - xi2 * eta3 - xi1 * eta4 + xi4 * eta1 )
	/
	( pow(xi4, 2) + pow(eta4, 2) );
}// eta_squared

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::eta_value( 
    transfer_function::impedance const& imp)
{
    return sqrt( eta_squared(imp) );
}// eta_value

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::eta_variance( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi4      = ip::xi4_value(imp);
    tf::real_type xi4_var  = ip::xi4_variance(imp);
    tf::real_type eta4     = ip::eta4_value(imp);
    tf::real_type eta4_var = ip::eta4_variance(imp);
    tf::real_type eta_sq    = eta_squared(imp);
    tf::real_type d32_var  = d_board_variance(3, 2, imp);
    tf::real_type d14_var  = d_board_variance(1, 4, imp);
    return 
	( d32_var + d14_var 
	  +
	  4 * eta_sq * ( pow(xi4, 2) * xi4_var + pow(eta4, 2) * eta4_var ) )
	/
	( 4 * eta_sq * pow( pow(xi4, 2) + pow(eta4, 2), 2 ) );
}// eta_variance

transfer_function::real_measurement_type
erro_manual::base_xi_eta::bahr_invariant::eta(
    transfer_function::impedance const& imp)
{
    return tf::real_measurement_type( eta_value(imp), eta_variance(imp) );
}//eta

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::Sigma_value( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi2_sq  = pow(ip::xi2_value(imp), 2);
    tf::real_type xi3_sq  = pow(ip::xi3_value(imp), 2);
    tf::real_type xi4_sq  = pow(ip::xi4_value(imp), 2);
    tf::real_type eta2_sq = pow(ip::eta2_value(imp), 2);
    tf::real_type eta3_sq = pow(ip::eta3_value(imp), 2);
    tf::real_type eta4_sq = pow(ip::eta4_value(imp), 2);
    
    return ( xi3_sq + eta3_sq + xi2_sq + eta2_sq ) / (xi4_sq + eta4_sq );
}// Sigma_value

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::Sigma_variance( 
    transfer_function::impedance const& imp)
{
    tf::real_type xi2_sq   = pow(ip::xi2_value(imp), 2);
    tf::real_type xi3_sq   = pow(ip::xi3_value(imp), 2);
    tf::real_type xi4_sq   = pow(ip::xi4_value(imp), 2);
    tf::real_type eta2_sq  = pow(ip::eta2_value(imp), 2);
    tf::real_type eta3_sq  = pow(ip::eta3_value(imp), 2);
    tf::real_type eta4_sq  = pow(ip::eta4_value(imp), 2);
    tf::real_type Sigma_sq = pow(Sigma_value(imp), 2);
    tf::real_type xi2_var  = ip::xi4_variance(imp);
    tf::real_type xi3_var  = ip::xi4_variance(imp);
    tf::real_type xi4_var  = ip::xi4_variance(imp);
    tf::real_type eta2_var = ip::eta4_variance(imp);
    tf::real_type eta3_var = ip::eta4_variance(imp);
    tf::real_type eta4_var = ip::eta4_variance(imp);

    return 
	( xi3_sq * xi3_var + eta3_sq * eta3_var 
	  + 
	  xi2_sq * xi2_var + eta2_sq * eta2_var 
	  + 
	  Sigma_sq * ( xi4_sq * xi4_var + eta4_sq * eta4_var ) ) * 4
	/
	pow( xi4_sq + eta4_sq, 2);	
}// Sigma_variance

transfer_function::real_measurement_type
erro_manual::base_xi_eta::bahr_invariant::Sigma(
    transfer_function::impedance const& imp)
{
    return tf::real_measurement_type( Sigma_value(imp), Sigma_variance(imp) );
}//Sigma

transfer_function::real_type
erro_manual::base_xi_eta::bahr_invariant::d_board_variance( 
    int i, int j, transfer_function::impedance const& imp)
{
    tf::real_type xi_i      = ip::xi_value(i, imp);
    tf::real_type xi_i_var  = ip::xi_variance(i, imp);
    tf::real_type xi_j      = ip::xi_value(j, imp);
    tf::real_type xi_j_var  = ip::xi_variance(j, imp);
    tf::real_type eta_i     = ip::eta_value(i, imp);
    tf::real_type eta_i_var = ip::eta_variance(i, imp);
    tf::real_type eta_j     = ip::eta_value(j, imp);
    tf::real_type eta_j_var = ip::eta_variance(j, imp);

    return 
	pow(xi_i, 2) * eta_j_var + pow(eta_j, 2) * xi_i_var
	+
	pow(eta_i, 2) * xi_j_var + pow(xi_j, 2) * eta_i_var;
}// d_board_variance

