#ifndef IMPORT_IMPL_HPP
#define IMPORT_IMPL_HPP

#include "import.hpp"
#include "../transfer_functions/transfer_functions.hpp"
#include <fstream>

namespace import
{
//! \brief verifica se \a file é do tipo egbert
/*! Um arquivo é identificado como do tipo egbert se conter uma linha que
  inicie com o texto " Transfer Functions" (o espaço no início é importante

  \param[in] file arquivo

  \return true se \a file for do tipo egbert
*/
bool
is_egbert( std::ifstream& file );


////////////////////////////////////////////////////////////////////////////////
//! \brief verifica se \a file é do tipo jones
/*! Um arquivo é identificado como do tipo jones se conter uma linha que
 *  inicie com o texto ">AZIMUTH ="
  
 \param[in] file arquivo

 \return true se \p file for do tipo jones
*/
bool
is_jones( std::ifstream& file );

// funções documentadas em import.hpp -----------------------------------
transfer_function::impedance_collection_type::size_type 
number_of_transfer_functions( std::ifstream& file );

void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& file );
    
} //namespace import

#endif
