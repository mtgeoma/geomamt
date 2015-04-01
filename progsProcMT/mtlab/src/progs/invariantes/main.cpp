#include "import/import.hpp"
#include "invariantes.hpp"
#include "invariants/dimensionalidade/dimensionalidade.hpp"
#include <cstdlib>

int
main(int argc, char* argv[])
{
  namespace po = boost::program_options;
  namespace tf = transfer_function;

  po::variables_map vm;  //mapa com os argumentos da linha de comando

  //Descrição das opções
  po::options_description options("Opções");

  program::command_line( vm, options, argc, argv );

  if ( argc == 2 && vm.count("help") )
  {
    program::print_help(argv[0], options);
    return EXIT_SUCCESS;
  }
  else if ( argc == 1 || argc > 5 )
  {
    program::print_help(argv[0], options);
    return EXIT_FAILURE;
  }

  try
  {
    dimensionalidade::thresholds threshold =
      dimensionalidade::inicializa_thresholds( vm[ "thresholds" ].
					       as< std::string >() );

    std::ifstream in( vm[ "file" ].as< std::string >().c_str() );
    if (!in) throw std::runtime_error("Não foi possível abrir o arquivo "
				      +
				      vm["file"].as< std::string>() );

    int n = import::number_of_transfer_functions( in );
    tf::impedance_collection_type impedances( n );
    import::import_impedances( impedances, in );
    program::lista_invariante( vm["invariant"].as< std::string >(),
			       impedances, threshold );
    return EXIT_SUCCESS;
  }
  catch( std::exception& e )
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

}// main(int argc, char*argv[])
