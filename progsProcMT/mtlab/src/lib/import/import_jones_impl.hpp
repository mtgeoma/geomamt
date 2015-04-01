#ifndef IMPORT_JONES_IMPL_HPP
#define IMPORT_JONES_IMPL_HPP

#include "import_jones.hpp"
#include <map> //component_map, components_collection

namespace import_jones
{
//! \brief Conjunto de todos os valores de uma determinada componente
/*! A chave do mapa é período é o valor mapeado é o valor da componente
 *  daquele período
 */
typedef std::map<
  transfer_function::real_type,
  transfer_function::complex_measurement_type > component_map;



//! \brief Coleção de todas as componentes de um arquivo Jones
/*!  A chave do mapa é um \p std::string que identifica a componente e o valor
  mapeado é um \p std::map com os todos os valores da componente identificada
  pela string
*/
typedef std::map< std::string, component_map > components_collection;

////////////////////////////////////////////////////////////////////////////////
//! \brief importa a componente \a component identificada por \a comp_string do arquivo \a file
/*! Os arquivos do tipo jones identificam cada componente por um
 *  string. Quando essa string é encontrada os valores são inseridos na
 *  componente fornecida

 \param[out] component mapa com os valores do componente
 \param[in] comp_string identificação do componente (case sensitive)
 \param[in] file arquivo tipo jones
*/void
import_component( component_map& component, const std::string comp_string,
		  std::ifstream& file);



////////////////////////////////////////////////////////////////////////////////
//! \brief confirma se os periodos do arquivo estão validos
/*! Se uma componente (Tzx, Tzy, Zxx, Zxy, etc.) tiver um certo período todas
 *  as demais componentes devem ter o mesmo período
 \todo não é claro o que o jones calcula os erros
 \param[in] collection vetor com as componentes a serem verificadas

 \return true se todas as componentes de \a collection tiverem os mesmos
 períodos
*/bool periods_are_valid( components_collection const& collection );



////////////////////////////////////////////////////////////////////////////////
//!!\brief transfere todos os componentes de \a file em \a collection
/*! collection contém todas as 42 componentes possíveis, as que não estão em
 *  \a file tem tamanho igual a zero.

 \param[out] collection
 \param[in] file
*/
void
read_components( components_collection& collection, std::ifstream& file);



////////////////////////////////////////////////////////////////////////////////
//! \brief Verifica se todos os periodos de \a first_it estão em \a second_it
/*!

  \note Não verifica se os containers de \a first_it e \a second_it
  tem o mesmo tamanho

  \param first_it
  \param second_it

  \return
*/bool
periods_are_equal( components_collection::const_iterator const first_it,
		   components_collection::const_iterator const second_it );



////////////////////////////////////////////////////////////////////////////////
//! \brief salta o bloco de comentário de \a file
/*!
  \pre \a file deve ser um arquivo válido no jones (file == true)
  \post \a file fica posicionado após o bloco de comentário
  \param file arquivo no formato jones
*/void
skip_comment_block( std::ifstream& file );



////////////////////////////////////////////////////////////////////////////////
//! \brief salta o bloco de informação de \a file
/*!

  \pre \a file deve estar posicionado no início do bloco de
  informação. Essa condição pode ser consiguida com uma chamada a \ref
  skip_comment_block

  \post \a file fica posicionado após o bloco de informação

  \param file arquivo no formato jones
*/
void
skip_information_block( std::ifstream& file );

// funções documentadas em import_jones.hpp -----------------------------------
transfer_function::impedance_collection_type::size_type
number_of_transfer_functions( std::ifstream& file );

void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& file );

} // namespace import_jones

#endif
