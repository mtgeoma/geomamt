#ifndef ERRO_MANUAL_BASE_XI_ETA_BAHR_INVARIANTS_HPP
#define ERRO_MANUAL_BASE_XI_ETA_BAHR_INVARIANTS_HPP

#include "transfer_functions/transfer_functions.hpp"

namespace erro_manual
{
    namespace base_xi_eta
    {
//! 
/*! \brief Funções para calcular os invariantes de Bahr
*/
	namespace bahr_invariant
	{
//! @name Valores de cada invariente
/*@{*/
	    transfer_function::real_type
	    kappa_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    mu_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    eta_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Sigma_value( transfer_function::impedance const& imp);
/*@}*/

//! @name Variança de cada Invariante
/*@{*/
	    transfer_function::real_type
	    kappa_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    mu_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    eta_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Sigma_variance( transfer_function::impedance const& imp);
/*@}*/

//! @name Medida (valor +/- incerteza) de cada invariante
/*@{*/
	    transfer_function::real_measurement_type
	    kappa( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    mu( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    eta( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    Sigma( transfer_function::impedance const& imp);
/*@}*/
	}// namespace bahr_invariant
    }// namespace base_xi_eta
}// namespace erro_manual

#endif
