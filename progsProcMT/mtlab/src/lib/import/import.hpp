#ifndef IMPORT_HPP
#define IMPORT_HPP

#include "import_impl.hpp"
#include "../transfer_functions/transfer_functions.hpp"
#include <fstream>

//! \brief funções para importar os valores de impedâncias de arquivos
/*! 
  
 */

namespace import{

////////////////////////////////////////////////////////////////////////////////
//! \brief Retorna o número de funções de transferência em \a file
/*!
  \pre \a file deve ser um arquivo válido
  \invariant a posição de leitura de \a file

  \param[in] file stream de um arquivo válido no formato jones

  \return o número de funções de transferência em \a file

  \todo implementar exception
*/
transfer_function::impedance_collection_type::size_type 
number_of_transfer_functions( std::ifstream& file );



////////////////////////////////////////////////////////////////////////////////
//! \brief função para adquirir as funções de transferência de um arquivo
//
/*!
  
  \param[out] tf container com as funções de transferência
  \param[in] file arquivo em um formato conhecido

*/
void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& file );

} //namespace import


#endif
