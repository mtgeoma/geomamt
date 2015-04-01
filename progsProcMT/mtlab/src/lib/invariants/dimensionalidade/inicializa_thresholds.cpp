#include "dimensionalidade.hpp"
#include "transfer_functions/transfer_functions.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

dimensionalidade::thresholds dimensionalidade::inicializa_thresholds( std::string const& file_thresholds )
{
  dimensionalidade::thresholds threshold;

  threshold.I3 = 0.1;
  threshold.I4 = 0.1;
  threshold.I5 = 0.1;
  threshold.I6 = 0.1;
  threshold.I7 = 0.1;
  threshold.Q  = 0.1;
  threshold.cos_alpha = 0.1;
  threshold.cos_beta = 0.1;
  threshold.kappa = 0.06;
  threshold.mu = 0.34;
  threshold.eta = 0.12;
  threshold.Sigma = 0.01;

  if( file_thresholds.size() != 0 )
    {
      try
	{
	  std::ifstream in( file_thresholds.c_str() );
	  if (!in) throw std::runtime_error(
					    "Não foi possível abrir o arquivo "
					    +
					    file_thresholds
					    +
					    "\nusando valores defaults");
	  std::string line;
	  while (getline(in,line))
	    {
	      // elimina todos os espaços em branco
	      boost::algorithm::erase_all(line, " ");
	      std::vector<std::string> str_vec;
	      boost::algorithm::split (str_vec,line,boost::algorithm::is_any_of("="));
	      if(str_vec[0]=="I3")
		{
		  threshold.I3 = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="I4")
		{
		  threshold.I4 = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="I5")
		{
		  threshold.I5 = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="I6")
		{
		  threshold.I6 = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="I7")
		{
		  threshold.I7 = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="Q")
		{
		  threshold.Q = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="cos_alpha")
		{
		  threshold.cos_alpha = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="cos_beta")
		{
		  threshold.cos_beta = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="kappa")
		{
		  threshold.kappa = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="mu")
		{
		  threshold.mu = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="eta")
		{
		  threshold.eta = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else if(str_vec[0]=="Sigma")
		{
		  threshold.Sigma = boost::lexical_cast<transfer_function::real_type>(str_vec[1]);
		}
	      else
		{
		  std::cout << "identificação de threshold desconhecida: " << str_vec[0] << "\n";
		}
	    }
	}
      catch( std::exception& e )
	{
	  std::cout << e.what() << std::endl;
	}
    }

  return  threshold;
}
