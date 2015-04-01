#ifndef ERRO_MANUAL_BASE_XI_ETA_INVARIANTS_PARAMETERS_HPP
#define ERRO_MANUAL_BASE_XI_ETA_INVARIANTS_PARAMETERS_HPP

#include "transfer_functions/transfer_functions.hpp"

namespace erro_manual
{
    namespace base_xi_eta
    {

//! Parâmetros de base para o cálculo dos invariantes
	namespace invariant_parameter
	{
//! @name Valores dos parâmetros base para o cálculos dos invariantes
/*@{*/
	    transfer_function::real_type
	    xi1_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi2_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi3_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi4_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta1_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta2_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta3_value( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta4_value( transfer_function::impedance const& imp );
/*@}*/

//! @name Variâncas dos parâmetros de base para o cálculo dos invariantes
/*@{*/	    transfer_function::real_type
	    xi1_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi2_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi3_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    xi4_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta1_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta2_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta3_variance( transfer_function::impedance const& imp );

	    transfer_function::real_type
	    eta4_variance( transfer_function::impedance const& imp );
/*@}*/

//! 
/*! @name Funções auxiliares para o cálculo dos invariantes
  
Essas funções simplesmente devolvem o valor (eta/xi)_value e as variânvias
(eta/xi)_variance dos paremetros \a xi(i) ou \a eta(i) onde \p i é o número do
parâmetro desejado. Quando possível use as funções diretas ????_value e
????_variance
*/
/*@{*/
	    transfer_function::real_type
	    xi_value(unsigned short int i, 
		     transfer_function::impedance const& imp);

	    transfer_function::real_type
	    xi_variance(unsigned short int i, 
			transfer_function::impedance const& imp);

	    transfer_function::real_type
	    eta_value(unsigned short int i, 
		      transfer_function::impedance const& imp);

	    transfer_function::real_type
	    eta_variance(unsigned short int i, 
			 transfer_function::impedance const& imp);
/*@}*/

	}// namespace invariant_parameter
    }// namespace base_xi_eta
}// namespace erro_manual

#endif
