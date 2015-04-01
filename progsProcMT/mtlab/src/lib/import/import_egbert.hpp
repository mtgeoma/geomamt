#ifndef IMPORT_EGBERT_HPP
#define IMPORT_EGBERT_HPP
#include <fstream>
#include "../transfer_functions/transfer_functions.hpp"

//! \brief funções para ler um arquivo do tipo egbert
namespace import_egbert 
{
//! \brief Retorna o número de funções de transferência em \p egbert_file
/*!
 * \pre \p egbert_file deve ser um arquivo válido no formato egbert
 * \invariant a posição de leitura de \p egbert_file

 \param egbert_file stream de um arquivo válido no formato egbert

 \return o número de funções de transferência em \p egbert_file
*/int 
number_of_transfer_functions( std::ifstream& egbert_file );
    
//! \brief função para adquirir as funções de transferência de um arquivo
//egbert
/*!
  
  \param[out] tf container com as funções de transferência
  \param[in] egbert_file arquivo egbert
*/
void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& egbert_file );

} //namespace import_egbert

#endif
