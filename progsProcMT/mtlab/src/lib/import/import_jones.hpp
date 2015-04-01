#ifndef IMPORT_JONES_HPP
#define IMPORT_JONES_HPP

#include <fstream>
#include "../transfer_functions/transfer_functions.hpp"

//! \brief funções para ler um arquivo do tipo jones
/*!

 */
namespace import_jones
{

////////////////////////////////////////////////////////////////////////////////
//! \brief Retorna o número de funções de transferência em \a file
/*!
  \pre \a file deve ser um arquivo válido no formato jones
  \invariant a posição de leitura de \a file

  \param[in] file stream de um arquivo válido no formato jones

  \return o número de funções de transferência em \a file
*/transfer_function::impedance_collection_type::size_type
number_of_transfer_functions( std::ifstream& file );



////////////////////////////////////////////////////////////////////////////////
//! \brief função para adquirir as funções de transferência de um arquivo
//jones
/*!

  \param[out] tf container com as funções de transferência
  \param[in] file arquivo no formato jones
*/
void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& file );

} //namespace import_jones

#endif
