#include "import/import.hpp"
#include "parametros-mt.hpp"
#include <cstdlib>

int
main(int argc, char* argv[])
{
  namespace po = boost::program_options;
  namespace tf = transfer_function;

  po::variables_map vm;  //mapa com os argumentos da linha de comando

  //Descrição das opções
  po::options_description options("Opções"); 

  try
  {
  program::command_line( vm, options, argc, argv );
  }
  catch(const std::exception& e)
  {
    std::cout << e.what() << "\n";
    return 1;
  }
  catch(...)
  {
    std::cout << "Unknown error\n";
    return 1;
  }

  if (vm.count("help"))
  {
    program::print_help(argv[0], options);
    return EXIT_SUCCESS;
  }

  if (!vm.count("file") || !vm.count("parameter"))
  {
    program::print_help(argv[0], options);
    return EXIT_FAILURE;
  }

  try
  {
    std::ifstream in( vm[ "file" ].as< std::string >().c_str() );
    if (!in) throw std::runtime_error("Não foi possível abrir o arquivo " 
				      + 
				      vm["file"].as< std::string>() );

    int n = import::number_of_transfer_functions( in );
    tf::impedance_collection_type impedances( n );
    import::import_impedances( impedances, in );

    program::lista_parametro( vm["parameter"].as< std::string >(), 
			      vm["azimute"].as< long double >(),
			      impedances );
    return EXIT_SUCCESS;
  }
  catch( std::exception& e )
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

}// main(int argc, char*argv[])
