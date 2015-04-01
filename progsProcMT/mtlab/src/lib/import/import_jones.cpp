#include "import_jones_impl.hpp"
#include "../transfer_functions/transfer_functions.hpp"
#include <fstream>
#include <ios> // std::pos_type
#include <map>
#include <cmath> //M_PIl

using namespace std;
using namespace transfer_function;

namespace import_jones
{

void
import_component( component_map& component, const std::string comp_string,
		  std::ifstream& file )
{
  file.seekg( 0 );

  string s;
  real_type period, real, imag, error, weight;
  component_map::size_type num_of_records;
  while(getline( file,  s ))
  {
    if (s.find(comp_string.c_str()) < s.length())
    {
      real_type escala = 1.0;
      if (s.find("field") < s.length())
      {
	escala = 4.0e-4 * M_PIl;
      }

      file >> num_of_records;
      for (component_map::size_type i = 0; i != num_of_records; ++i)
      {
	file >> period >> real >> imag >> error >> weight;
	if (period < 0.0)
	{
	  period = abs( 1.0 / period );
	}
	real *= escala;
	imag *= escala;
	error *= escala;
	component[period] =
	  complex_measurement_type(measurement(real, error * error / 2.0),
				   measurement(imag, error * error / 2.0) );
      }
    }
  }

  // se o encontrar eof() então todo o arquivo foi lido e está tudo
  // bem.  Sempre que eof() é encontrado o stream é marcado como
  // fail(), por isso é necessário o clear() antes de recolocar o
  // arquivo na posição inicial
  if ( file.eof() )
  {
    file.clear();
    file.seekg( 0 );
  }
} // import_component

void
read_components(components_collection& collection, std::ifstream& file)
{
  const int num_of_tags = 42;
  string const tag[num_of_tags] = {
    "RXX", "RXY", "RYX", "RYY", "RTE", "RTM", "RAV", "RDE",
    "SXX", "SXY", "SYX", "SYY", "STE", "STM", "SAV", "SDE",
    "ZXX", "ZXY", "ZYX", "ZYY", "ZTE", "ZTM", "ZAV", "ZDE",
    "QXX", "QXY", "QYX", "QYY", "QTE", "QTM", "QAV", "QDE",
    "CXX", "CXY", "CYX", "CYY", "CTE", "CTM", "CAV", "CDE",
    "TZX", "TZY"
  };

  for( int i = 0; i != num_of_tags; ++i )
  {
    import_component( collection[tag[i]], tag[i], file );
  }
}

bool
periods_are_valid( components_collection const& collection)
{
  components_collection::const_iterator const end = collection.end();
  components_collection::const_iterator it = collection.begin();

  // encontra a primeira coleção não nula
  while( it != end && it->second.size() == 0 )
  {
    ++it;
  }
  component_map::size_type const num_of_periods = it->second.size();
  components_collection::const_iterator const first_component = it;

  bool file_is_valid = true;
  if ( it != end ) ++it;

  // verifica se todos as coleções tem o mesmo tamanho da primeira e
  // se os períodos são iguais
  while( it != end && file_is_valid )
  {
    if ( it->second.size() != 0 )
    {
      if (it->second.size() == num_of_periods )
      {
	file_is_valid = periods_are_equal( first_component, it );
      }
      else
      {
	file_is_valid = false;
      }
    }
    ++it;
  }

  return file_is_valid;
} // periods_are_valid

bool
periods_are_equal(components_collection::const_iterator const first_it,
		  components_collection::const_iterator const second_it )
{
  component_map::const_iterator first_map_it = first_it->second.begin();
  component_map::const_iterator second_map_it = second_it->second.begin();
  component_map::const_iterator const end = first_it->second.end();

  bool is_valid = true;
  while ( ( first_map_it != end ) && is_valid )
  {
    is_valid = ( first_map_it->first == second_map_it->first );
    ++first_map_it;
    ++second_map_it;
  }
  return is_valid;
} // periods_are_equal

void
skip_comment_block( std::ifstream& file )
{
  file.seekg( 0 );
  string s;

  while ( file.peek() == '#' || file.peek() == '\n' )
  {
    getline( file, s );
  }
} //skip_comment_block

void
skip_information_block( std::ifstream& file )
{
  string s;

  while ( file.peek() == '>' )
  {
    getline( file, s );
  }
  // ignore first line after last '>'
  getline( file, s );
} //skip_information_block

transfer_function::impedance_collection_type::size_type
number_of_transfer_functions( std::ifstream& file )
{
  ifstream::pos_type initial_position = file.tellg();

  skip_comment_block( file );
  skip_information_block( file );

  string s;
  transfer_function::impedance_collection_type::size_type x = 0;
  getline( file, s );
  file >> x;

  file.seekg( initial_position );
  return x;
} // number_of_transfer_functions

void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& file )
{
  components_collection cc;
  read_components( cc, file );

  impedance_tensor_type impedanse_tensor(2, 2);
  real_type period = 0.0;
  int i = 0;
  component_map::const_iterator const end = cc["ZXX"].end();
  for ( component_map::const_iterator it = cc["ZXX"].begin();
	it != end; ++it )
  {
    period = it->first;
    impedanse_tensor(0, 0) = cc["ZXX"][period];
    impedanse_tensor(0, 1) = cc["ZXY"][period];
    impedanse_tensor(1, 0) = cc["ZYX"][period];
    impedanse_tensor(1, 1) = cc["ZYY"][period];
    tf[i] = impedance(period, impedanse_tensor);
    ++i;
  }
} // import_impedances

} // namespace import_jones
