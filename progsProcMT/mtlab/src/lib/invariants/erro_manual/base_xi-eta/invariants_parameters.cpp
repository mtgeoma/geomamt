#include "invariants_parameters.hpp"

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi1_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().real().value + imp.Zyy().real().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi1_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().real().variance + imp.Zyy().real().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi2_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().real().value + imp.Zyx().real().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi2_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().real().variance + imp.Zyx().real().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi3_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().real().value - imp.Zyy().real().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi3_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().real().variance + imp.Zyy().real().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi4_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().real().value - imp.Zyx().real().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi4_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().real().variance + imp.Zyx().real().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta1_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().imag().value + imp.Zyy().imag().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta1_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().imag().variance + imp.Zyy().imag().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta2_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().imag().value + imp.Zyx().imag().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta2_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().imag().variance + imp.Zyx().imag().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta3_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().imag().value - imp.Zyy().imag().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta3_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxx().imag().variance + imp.Zyy().imag().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta4_value( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().imag().value - imp.Zyx().imag().value ) / 2;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta4_variance( 
    transfer_function::impedance const& imp )
{
    return ( imp.Zxy().imag().variance + imp.Zyx().imag().variance ) / 4.0;
}

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi_value(
    unsigned short int i, transfer_function::impedance const& imp)
{
    switch (i)
    {
	case 1:
	    return xi1_value( imp );
	case 2:
	    return xi2_value( imp );
	case 3:
	    return xi3_value( imp );
	case 4:
	    return xi4_value( imp );
	default:
	    throw std::out_of_range( "invariant_parameter::xi vai de 1 a 4" );
    }
} // xi(i, imp)

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::xi_variance(
    unsigned short int i, transfer_function::impedance const& imp)
{
    switch (i)
    {
	case 1:
	    return xi1_variance( imp );
	case 2:
	    return xi2_variance( imp );
	case 3:
	    return xi3_variance( imp );
	case 4:
	    return xi4_variance( imp );
	default:
	    throw std::out_of_range( "invariant_parameter::xi vai de 1 a 4" );
    }
} // xi(i, imp)

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta_value(
    unsigned short int i, transfer_function::impedance const& imp)
{
    switch (i)
    {
	case 1:
	    return eta1_value( imp );
	case 2:
	    return eta2_value( imp );
	case 3:
	    return eta3_value( imp );
	case 4:
	    return eta4_value( imp );
	default:
	    throw std::out_of_range( "invariant_parameter::eta vai de 1 a 4" );
    }
} // eta(i, imp)

transfer_function::real_type
erro_manual::base_xi_eta::invariant_parameter::eta_variance(
    unsigned short int i, transfer_function::impedance const& imp)
{
    switch (i)
    {
	case 1:
	    return eta1_variance( imp );
	case 2:
	    return eta2_variance( imp );
	case 3:
	    return eta3_variance( imp );
	case 4:
	    return eta4_variance( imp );
	default:
	    throw std::out_of_range( "invariant_parameter::eta vai de 1 a 4" );
    }
} // eta(i, imp)
