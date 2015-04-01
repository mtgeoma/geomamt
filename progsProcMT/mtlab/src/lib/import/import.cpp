#include "import.hpp"
#include "import_jones.hpp"
#include "import_egbert.hpp"

transfer_function::impedance_collection_type::size_type 
import::number_of_transfer_functions( std::ifstream& file )
{
  if ( is_egbert( file ) )
  {
    return import_egbert::number_of_transfer_functions( file );
  }
  else if ( is_jones( file ) )
  {
    return import_jones::number_of_transfer_functions( file );
  }
  throw std::runtime_error("formato de arquivo desconhecido");
} //number_of_transfer_functions( std::ifstream& file )

void
import::import_impedances(transfer_function::impedance_collection_type& tf,
			  std::ifstream& file )
{
  if ( is_egbert( file ) )
  {
    import_egbert::import_impedances( tf, file );
  }
  else if ( is_jones( file ) )
  {
    import_jones::import_impedances( tf, file );
  }
  //deveria lançar uma exceção se chegasse aqui
    
} // import_impedances

bool
import::is_egbert( std::ifstream& file )
{
  std::string s;
  file.seekg( 0 );
  while(getline( file,  s ))
  {
    if ( s == " Transfer Functions" )
    {
      file.seekg(0);
      return true;
    }
  }

  // se o encontrar eof() então todo o arquivo foi lido e está tudo bem.
  // Sempre que eof() é encontrado o stream é marcado como fail(), por isso é
  // necessário o clear() antes de recolocar o arquivo na posição inicial
  if ( file.eof() )
  {
    file.clear();
    file.seekg( 0 );
  }
  return false;
}

bool
import::is_jones( std::ifstream& file )
{
  file.seekg( 0 );

  //string que identifica um arquivo jones
  std::string const tag = ">AZIMUTH   =";

  std::string s;
  while(getline( file,  s ) )
  {
    if ( s.size() >= tag.size() )
    {
      if ( s.substr(0, 12) == tag )
      {
	file.seekg(0);
	return true;
      }
    }
  }

  // se o encontrar eof() então todo o arquivo foi lido e está tudo bem.
  // Sempre que eof() é encontrado o stream é marcado como fail(), por isso é
  // necessário o clear() antes de recolocar o arquivo na posição inicial
  if ( file.eof() )
  {
    file.clear();
    file.seekg( 0 );
  }
  return false;
}
