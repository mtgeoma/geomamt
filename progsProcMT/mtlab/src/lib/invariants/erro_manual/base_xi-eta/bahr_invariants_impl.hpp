#ifndef ERRO_MANUAL_BASE_XI_ETA_BAHR_INVARIANTS_IMPL_HPP
#define ERRO_MANUAL_BASE_XI_ETA_BAHR_INVARIANTS_IMPL_HPP

namespace erro_manual
{
    namespace base_xi_eta
    {
	namespace bahr_invariant
	{

//! @name Funções auxilíares para o cálculos dos invariants de Bahr
/*@{*/
	    transfer_function::real_type
	    kappa_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    mu_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    eta_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    d_board_variance( int i, int j, 
			      transfer_function::impedance const& imp);
/*@}*/
	}// namespace bahr_invariant
    }// namespace base_xi_eta
}// namespace erro_manual


#endif
