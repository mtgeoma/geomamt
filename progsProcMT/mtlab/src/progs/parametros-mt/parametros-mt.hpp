#ifndef INVARIANTES_HPP
#define INVARIANTES_HPP

#include <iostream>
#include <boost/program_options.hpp>
#include <string>
#include "transfer_functions/transfer_functions.hpp"

//! Funções usadas por main()
namespace program
{
//! 
/*! \brief Verifica se a linha de comando está ok.
  
Se houver algum problema com a linha de comando, ou o programa tenha sido
executado com a opção \p --help exibe a mensagem de ajuda e retorna \a
false. Caso contrário retorna \a true e vm contem as infomações do arquivo em
\a vm["file"] e do invariante em \a vm["invariant"]

\param[out] vm mapa com as opções de linha de comando
\param[out] visible descrições das opções da linha de comando
\param[in] argc número de argumentos da linha de comando
\param[in] argv os argumentos da linha de comando
*/
    void command_line( boost::program_options::variables_map& vm, 
			  boost::program_options::options_description& visible,
			  int argc, char* argv[] );

    void print_help( std::string program_name, 
		     boost::program_options::options_description& visible );

    void lista_parametro( 
	std::string const& par,  
	long double azimute,
	transfer_function::impedance_collection_type& imp );
}//namespace program
#endif
