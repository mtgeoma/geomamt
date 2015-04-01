/*!
  \file   wal_invariants.hpp
  \author Marcos Banik de Padua
  \date   Mon Jun  4 17:19:37 2007
  
  \brief  WAL Invariants
  
  
*/
#ifndef ERRO_MANUAL_BASE_XI_ETA_WAL_INVARIANTS_HPP
#define ERRO_MANUAL_BASE_XI_ETA_WAL_INVARIANTS_HPP

#include "transfer_functions/transfer_functions.hpp"

namespace erro_manual
{
    namespace base_xi_eta
    {

//! Funções para calcular os invariantes WAL
	namespace wal_invariant
	{
//! @name Valores de cada invariente (sem as incertezas)
/*@{*/
	    transfer_function::real_type
	    I1_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I2_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I3_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I4_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I5_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I6_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I7_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Q_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    cos_alpha_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    cos_beta_value( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    N_value( transfer_function::impedance const& imp);
/*@}*/

//! @name Variança de cada Invariante
/*@{*/
	    transfer_function::real_type
	    I1_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I2_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I3_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I4_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I5_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I6_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    I7_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    Q_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    cos_alpha_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    cos_beta_variance( transfer_function::impedance const& imp);

	    transfer_function::real_type
	    N_variance( transfer_function::impedance const& imp);
/*@}*/

//! @name Medida (valor +/- incerteza) de cada invariante
/*@{*/
	    transfer_function::real_measurement_type
	    I1( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    I2( transfer_function::impedance const& imp);
	    
	    transfer_function::real_measurement_type
	    I3( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    I4( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    I5( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    I6( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    I7( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    Q( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    cos_alpha( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    cos_beta( transfer_function::impedance const& imp);

	    transfer_function::real_measurement_type
	    N( transfer_function::impedance const& imp);
/*@}*/

	}// namespace wal_invariant
    }// namespace base_xi_eta
}// namespace erro_manual
#endif
