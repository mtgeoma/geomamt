#include "invariantes.hpp"
#include "invariants/erro_manual/base_xi-eta/wal_invariants.hpp"
#include "invariants/erro_manual/base_xi-eta/bahr_invariants.hpp"
#include "invariants/dimensionalidade/dimensionalidade.hpp"
#include <iostream> //cout

void program::command_line(boost::program_options::variables_map& vm,
			   boost::program_options::options_description& visible,
			   int argc, char *argv[])
{
  namespace po = boost::program_options;

  // Declara as opções genericas
  visible.add_options()
    ("help", "exibe essa ajuda")
    ("thresholds", po::value< std::string >()->default_value(""),
     "arquivo com valores dos thresholds");

  // Declara as opções posicionais. As opções que dependem apenas da posição
  //devem ficar ocultas do usuário.
  po::options_description hidden("Opções Ocultas");
  hidden.add_options()
    ("file", po::value< std::string >(), "")
    ("invariant", po::value< std::string >(), "");

  // Estabelece a ordem das opções posicionais.
  po::positional_options_description pos_order;
  pos_order.add("file", 1);
  pos_order.add("invariant", 1);

  // As opções de linha de comando ( ocultas + genéricas )
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  //Descrição da linha de comando para o programa principal

  // interpreta as opções
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).
	    positional(pos_order).run(), vm);

}// command_line_ok

void
program::print_help(std::string program_name,
		    boost::program_options::options_description& visible)
{
  std::cerr << "Uso: " << program_name << " ARQUIVO INVARIANTE\n";
  std::cerr << "ARQUIVO é o arquivo com os tensores de impedância, e\n"
	    << "INVARIANTE é um invariante de WAL: I[1 a 7], Q, cos_alpha e "
	    << "cos_beta ou\n"
	    << "um dos invariantes de Bahr: kappa, mu, eta ou Sigma\n"
	    << "ou N (invariante de inomogeneidade)\n"
	    << "ou dim-wal, dim-bahr-q (verifica a dimensionalidade)\n";
  std::cerr << visible;
}

void program::lista_invariante(
    std::string const& inv,
    transfer_function::impedance_collection_type const& imp,
    dimensionalidade::thresholds threshold)
{
  namespace tf   = transfer_function;
  namespace wal  = erro_manual::base_xi_eta::wal_invariant;
  namespace bahr = erro_manual::base_xi_eta::bahr_invariant;

  for( tf::impedance_collection_type::size_type i = 0;
       i != imp.size(); ++i )
  {
    if ( inv == "I1" )
    {
      std::cout << imp[i].period << "\t" << wal::I1(imp[i]).value << "\t"
	   << wal::I1(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I2" )
    {
      std::cout << imp[i].period << "\t" << wal::I2(imp[i]).value << "\t"
	   << wal::I2(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I3" )
    {
      std::cout << imp[i].period << "\t" << wal::I3(imp[i]).value << "\t"
	   << wal::I3(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I4" )
    {
      std::cout << imp[i].period << "\t" << wal::I4(imp[i]).value << "\t"
	   << wal::I4(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I5" )
    {
      std::cout << imp[i].period << "\t" << wal::I5(imp[i]).value << "\t"
	   << wal::I5(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I6" )
    {
      std::cout << imp[i].period << "\t" << wal::I6(imp[i]).value << "\t"
	   << wal::I6(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "I7" )
    {
      std::cout << imp[i].period << "\t" << wal::I7(imp[i]).value << "\t"
	   << wal::I7(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "Q" )
    {
      std::cout << imp[i].period << "\t" << wal::Q(imp[i]).value << "\t"
	   << wal::Q(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "cos_alpha" )
    {
      std::cout << imp[i].period << "\t" << wal::cos_alpha(imp[i]).value << "\t"
	   << wal::cos_alpha(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "cos_beta" )
    {
      std::cout << imp[i].period << "\t" << wal::cos_beta(imp[i]).value << "\t"
	   << wal::cos_beta(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "N" )
    {
      std::cout << imp[i].period << "\t" << wal::N(imp[i]).value << "\t"
	   << wal::N(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "kappa" )
    {
      std::cout << imp[i].period << "\t" << bahr::kappa(imp[i]).value << "\t"
	   << bahr::kappa(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "mu" )
    {
      std::cout << imp[i].period << "\t" << bahr::mu(imp[i]).value << "\t"
	   << bahr::mu(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "eta" )
    {
      std::cout << imp[i].period << "\t" << bahr::eta(imp[i]).value << "\t"
	   << bahr::eta(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "Sigma" )
    {
      std::cout << imp[i].period << "\t" << bahr::Sigma(imp[i]).value << "\t"
	   << bahr::Sigma(imp[i]).uncertainty() << "\n";
    }
    else if ( inv == "dim-wal" )
    {
      std::cout << imp[i].period << "\t" << dimensionalidade::wal(imp[i], threshold)
	   << "\n";
    }
    else if ( inv == "dim-bahr-q" )
    {
      std::cout << imp[i].period << "\t"
	   << dimensionalidade::bahr_Q(imp[i], threshold) << "\n";
    }
    else
    {
      throw std::runtime_error("Invariante desconhecido: " + inv );
    }
  }
} // lista_invariante
