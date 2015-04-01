#ifndef ERRO_MANUAL_BASE_XI_ETA_WAL_INVARIANTS_IMPL_HPP
#define ERRO_MANUAL_BASE_XI_ETA_WAL_INVARIANTS_IMPL_HPP

namespace erro_manual
{
    namespace base_xi_eta
    {
	namespace wal_invariant
	{

//! 
/*! @name Quadrados dos Invariantes
Nas expressões de cálculo das varianças e dos valores dos invariantes
  várias vezes aparecem termos que envolvem os quadrados dos invariantes que
  seguem, daí a sua necessidade  
*/
/*@{*/
	    transfer_function::real_type
	    I1_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I2_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Rre_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Rim_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I5_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I6_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I7_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Q_squared( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    N_squared( transfer_function::impedance const& imp);
/*@}*/

//! 
/*! @name Funções auxiliares
  Essas funções auxiliam no cálculo de vários invariantes, sobretudo I7 e Q
*/
/*@{*/
	    transfer_function::real_type 
	    d_value( unsigned short int i, unsigned short int j, 
		     transfer_function::impedance const& imp);

	    transfer_function::real_type 
	    d_variance( unsigned short int i, unsigned short int j, 
			transfer_function::impedance const& imp);
/*@}*/
	}// namespace wal_invariant
    }// namespace base_xi_eta
}// namespace erro_manual


#endif
