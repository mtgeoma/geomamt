#ifndef IMPORT_EGBERT_IMPL_HPP
#define IMPORT_EGBERT_IMPL_HPP

#include "import_egbert.hpp"
#include "../transfer_functions/transfer_functions.hpp"

namespace import_egbert
{
namespace tf = transfer_function;

//!
/*! \brief Le o valor do próximo tipper a partir da posição atual de \a
 *  egbert_file

 \param tv
 \param egbert_file
*/void
import_tipper_vector_value(std::vector< tf::complex_type >& tv,
			   std::ifstream& egbert_file );

//!
/*! \brief Le o valor do próximo tensor de impedância a partir da posição
 *  atual de \a egbert_file

 \param imp_tensor tensor de impedância
 \param egbert_file arquivo do tipo egbert
*/void
import_impedance_tensor_value(tf::complex_matrix_type& imp_tensor,
			      std::ifstream& egbert_file );

//!
/*! \brief Le o valor da próxima "inverse coherent signal power matrix" a
 *  partir da posição de \a egbert_file

 \param matrix inverse coherent signal power matrix
 \param egbert_file arquivo do tipo egbert
*/
void
import_inverse_coherent_signal_power_matrix(tf::complex_matrix_type& matrix,
					    std::ifstream& egbert_file );

//!
/*! \brief Le o valor da próxima matriz de covariância residual a partir da
 *  posição de \a egbert_file

 \param[out] matrix a matriz de covariância residual
 \param egbert_file arquivo do tipo egbert
*/
void
import_residual_covariance_matrix(tf::complex_matrix_type& matrix,
				  std::ifstream& egbert_file );

//!
/*! \brief Le, a partir da posição de \a egbert_file o valor do período
 *  tipper e da impedância
 */

/*!

  \param period
  \param tipper
  \param impedance
  \param egbert_file
*/void
import_transfer_function(int const number_of_channels,
			 tf::real_type& period,
			 tf::tipper_vector_type& tipper,
			 tf::impedance_tensor_type& impedance,
			 std::ifstream& egbert_file );

//! Constrói o tipper (valor e incerteza) através do seu valor e das matrizes de
/*! covariância

  \param[out] tipper
  \param[in] tipper_value
  \param[in] s_matrix
  \param[in] n_matrix
*/
void
build_tipper(tf::tipper_vector_type& tipper,
	     std::vector< tf::complex_type > const& tipper_value,
	     tf::complex_matrix_type const& s_matrix,
	     tf::complex_matrix_type const& n_matrix );

//! Constrói o tensor de impedância (valor e incerteza)
/*!

  \note O egbert não justifica porque usa apenas a parte real das
  matrizes S e N para o cálculo das varianças

  \param[out] imp
  \param[in] impedance_value
  \param[in] s_matrix
  \param[in] n_matrix
*/
void
build_impedance(tf::impedance_tensor_type& imp,
		tf::complex_matrix_type const& impedance_value,
		tf::complex_matrix_type const& s_matrix,
		tf::complex_matrix_type const& n_matrix );

////////////////////////////////////////////////////////////////////////////////
// Funções documentadas em import_egbert.hpp
int
number_of_channels( std::ifstream& egbert_file );

int
number_of_transfer_functions( std::ifstream& egbert_file );

void
import_impedances(tf::impedance_collection_type& tf,
		  std::ifstream& egbert_file );

} // namespace import_egbert

#endif
