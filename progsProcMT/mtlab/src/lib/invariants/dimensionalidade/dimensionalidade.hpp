/*!
  \file   dimensionalidade.hpp
  \author Marcelo Banik de Padua
  \date   Mon Set  2 15:10:00 2008
  
  \brief  Dimensionality
  
  
*/
#ifndef DIMENSIONALITY_HPP
#define DIMENSIONALITY_HPP

#include "transfer_functions/transfer_functions.hpp"
#include <string>

//#include "invariants/erro_manual/base_xi-eta/wal_invariants.hpp"

//! Funções para calcular a dimensionalidade conforme WAL ou Bahr-Q
namespace dimensionalidade
{
//! @name Define a estrutura dimensao (por enquanto é apenas um comentário num string)
/*@{*/
  typedef std::string dimensao;
/*@}*/

//! @name Define a estrutura thresholds com os limites para os invariantes
/*@{*/
   typedef struct
   {
     transfer_function::real_type I3;
     transfer_function::real_type I4;
     transfer_function::real_type I5;
     transfer_function::real_type I6;
     transfer_function::real_type I7;
     transfer_function::real_type Q;
     transfer_function::real_type cos_alpha;
     transfer_function::real_type cos_beta;
     transfer_function::real_type kappa;
     transfer_function::real_type mu;
     transfer_function::real_type eta;
     transfer_function::real_type Sigma;
   }
   thresholds;
/*@}*/

//! @name função que inicializa os valores dos thresholds
/*@{*/
thresholds inicializa_thresholds( std::string const& file_thresholds );
/*@}*/
  
//! @name métodos de classificar a dimensionalidade
/*@{*/
  dimensao
  wal( transfer_function::impedance const& imp, thresholds threshold);
  dimensao
  bahr_Q( transfer_function::impedance const& imp, thresholds threshold);
/*@}*/

}// namespace dimensionalidade
#endif
