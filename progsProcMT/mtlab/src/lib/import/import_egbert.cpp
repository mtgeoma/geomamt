#include "import_egbert_impl.hpp"
#include "../transfer_functions/transfer_functions.hpp"
#include <cmath> //M_PIl
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace transfer_function;

namespace import_egbert
{

int
number_of_channels(std::ifstream& egbert_file )
{
  int number_of_channels_;
  ifstream::pos_type initial_position;
  string s;
  initial_position = egbert_file.tellg();

  while ( getline( egbert_file,  s) )
  {
    if (s.find("number of channels") < s.length())
    {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep(" \t");
      tokenizer tok(s, sep);
      for(tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
      {
	if(*it == "channels")
	{
	  ++it;
	  number_of_channels_ = boost::lexical_cast<int>(*it);
	  break;
	}
      }
      break;
    }
  }
  // se o encontrar eof() então todo o arquivo foi lido e está tudo bem.
  // Sempre que eof() é encontrado o stream é marcado como fail(), por isso é
  // necessário o clear() antes de recolocar o arquivo na posição inicial
  if ( egbert_file.eof() ) {
    egbert_file.clear();
    egbert_file.seekg( initial_position );
  }

  egbert_file.seekg( initial_position, ios::beg );
  return number_of_channels_;
} //number_of_channels

int
number_of_transfer_functions(std::ifstream& egbert_file )
{
  int number_of_transfer_functions_ = 0;
  ifstream::pos_type initial_position;
  string s;
  initial_position = egbert_file.tellg();

  while(getline( egbert_file,  s ))
  {
    if (s.find("period :") < s.length()) {
      ++number_of_transfer_functions_;
    }
  }

  // se o encontrar eof() então todo o arquivo foi lido e está tudo bem.
  // Sempre que eof() é encontrado o stream é marcado como fail(), por isso é
  // necessário o clear() antes de recolocar o arquivo na posição inicial
  if ( egbert_file.eof() ) {
    egbert_file.clear();
    egbert_file.seekg( initial_position );
  }

  egbert_file.seekg( initial_position, ios::beg );
  return number_of_transfer_functions_;
} //number_of_transfer_functions

void
import_impedances(transfer_function::impedance_collection_type& tf,
		  std::ifstream& egbert_file )
{
  int Nchannels = number_of_channels( egbert_file );
  int end = number_of_transfer_functions( egbert_file );
  real_type period;
  tipper_vector_type tipper(2);
  impedance_tensor_type impedance_tensor(2, 2);

  for (int i = 0; i != end; ++i) 
  {
    import_transfer_function(Nchannels, period, tipper, impedance_tensor, egbert_file);
    tf[i] = impedance(period, impedance_tensor);
  }
}

void
import_tipper_vector_value(std::vector< transfer_function::complex_type >& tv,
			   std::ifstream& egbert_file )
{
  real_type real, imag;
  egbert_file >> real >> imag;

  real_type escala = 4.0e-4 * M_PIl;
  tv.push_back(complex_type( real * escala, imag * escala));
  egbert_file >> real >> imag;
  tv.push_back(complex_type( real * escala, imag * escala));
}

void
import_impedance_tensor_value(complex_matrix_type& imp_tensor,
			      std::ifstream& egbert_file )
{
  real_type real, imag;
  egbert_file >> real >> imag;

  real_type escala = 4.0e-4 * M_PIl;
  imp_tensor( 0, 0 ) = complex_type( real * escala, imag * escala);
  egbert_file >> real >> imag;
  imp_tensor( 0, 1 ) = complex_type( real * escala, imag * escala);
  egbert_file >> real >> imag;
  imp_tensor( 1, 0 ) = complex_type( real * escala, imag * escala);
  egbert_file >> real >> imag;
  imp_tensor( 1, 1 ) = complex_type( real * escala, imag * escala);
}

void
import_inverse_coherent_signal_power_matrix(complex_matrix_type& matrix,
					    std::ifstream& egbert_file )
{
  real_type real, imag;
  egbert_file >> real >> imag;

  real_type escala = 4.0e-4 * M_PIl;
  matrix(0, 0) = complex_type(real * escala, imag * escala);
  egbert_file >> real >> imag;
  matrix(1, 0) = complex_type(real * escala, imag * escala);
  egbert_file >> real >> imag;
  matrix(1, 1) = complex_type(real * escala, imag * escala);

  //completa a matrix hermetiana
  matrix(0, 1) = conj(matrix(1, 0));
}

void
import_residual_covariance_matrix(complex_matrix_type& matrix,
				  std::ifstream& egbert_file )
{
  real_type real, imag;

  real_type escala = 4.0e-4 * M_PIl;

  if( gmm::mat_nrows(matrix)==3 )
  {
    egbert_file >> real >> imag;
    matrix(0, 0) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(1, 0) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(1, 1) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(2, 0) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(2, 1) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(2, 2) = complex_type(real * escala, imag * escala);

    //completa a matriz hermetiana
    matrix(0, 1) = conj(matrix(1, 0));
    matrix(0, 2) = conj(matrix(2, 0));
    matrix(1, 2) = conj(matrix(2, 1));
  }
  else
  {
    egbert_file >> real >> imag;
    matrix(0, 0) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(1, 0) = complex_type(real * escala, imag * escala);
    egbert_file >> real >> imag;
    matrix(1, 1) = complex_type(real * escala, imag * escala);

    //completa a matriz hermetiana
    matrix(0, 1) = conj(matrix(1, 0));
  }
}

void
build_tipper(transfer_function::tipper_vector_type& tipper,
	     std::vector< transfer_function::complex_type > const& tipper_value,
	     complex_matrix_type const& s_matrix,
	     complex_matrix_type const& n_matrix )
{
  long double realVal = tipper_value[0].real();
  long double realVar = n_matrix(0, 0).real() * s_matrix(0, 0).real() / 2.0;

  long double imagVal = tipper_value[0].imag();
  long double imagVar = n_matrix(0, 0).real() * s_matrix(0, 0).real() / 2.0;

  tipper[0] = 
    complex_measurement_type(real_measurement_type(realVal, realVar),
			     real_measurement_type(imagVal, imagVar));

  realVal = tipper_value[1].real();
  realVar = n_matrix(0, 0).real() * s_matrix(1, 1).real() / 2.0;

  imagVal = tipper_value[1].imag();
  imagVar = n_matrix(0, 0).real() * s_matrix(1, 1).real() / 2.0;

  tipper[1] = 
    complex_measurement_type(real_measurement_type(realVal, realVar),
			     real_measurement_type(imagVal, imagVar));
}//build_tipper

void
build_impedance(transfer_function::impedance_tensor_type& imp,
		complex_matrix_type const& impedance_value,
		complex_matrix_type const& s_matrix,
		complex_matrix_type const& n_matrix )
{
  int n;
  if( gmm::mat_nrows(n_matrix)==3 )
  {
    n=1;
  }
  else
  {
    n=0;
  }
  real_type real_value = impedance_value(0, 0).real();
  real_type real_variance = n_matrix(0+n, 0+n).real() * s_matrix(0, 0).real() / 2.0;
  real_type imag_value = impedance_value(0, 0).imag();
  real_type imag_variance = real_variance;
  real_measurement_type real_measure(real_value, real_variance);
  real_measurement_type imag_measure(imag_value, imag_variance);
  imp(0, 0) = complex_measurement_type(real_measure, imag_measure);

  real_value = impedance_value(0, 1).real();
  real_variance = n_matrix(0+n, 0+n).real() * s_matrix(1, 1).real() / 2.0;
  imag_value = impedance_value(0, 1).imag();
  imag_variance = real_variance;
  real_measure = real_measurement_type(real_value, real_variance);
  imag_measure = real_measurement_type(imag_value, imag_variance);
  imp(0, 1) = complex_measurement_type(real_measure, imag_measure);

  real_value = impedance_value(1, 0).real();
  real_variance = n_matrix(1+n, 1+n).real() * s_matrix(0, 0).real() / 2.0;
  imag_value = impedance_value(1, 0).imag();
  imag_variance = real_variance;
  real_measure = real_measurement_type(real_value, real_variance);
  imag_measure = real_measurement_type(imag_value, imag_variance);
  imp(1, 0) = complex_measurement_type(real_measure, imag_measure);

  real_value = impedance_value(1, 1).real();
  real_variance = n_matrix(1+n, 1+n).real() * s_matrix(1, 1).real() / 2.0;
  imag_value = impedance_value(1, 1).imag();
  imag_variance = real_variance;
  real_measure = real_measurement_type(real_value, real_variance);
  imag_measure = real_measurement_type(imag_value, imag_variance);
  imp(1, 1) = complex_measurement_type(real_measure, imag_measure);

}//build_impedance

void
import_transfer_function(int const number_of_channels,
			 transfer_function::real_type& period,
			 transfer_function::tipper_vector_type& tipper,
			 transfer_function::impedance_tensor_type& impedance,
			 std::ifstream& egbert_file )
{

  string s;

  //s_matrix e n_matrix são hermetianas mas é melhor construir essas
  //matrizes na marra até aprender a acessar um elemento de
  //ublas::hermitic_matrix facilmente

  // para usar os elementos de uma matriz hermetiana diretamente
  // seria necessário um const_cast como nas linhas abaixo
    
  // typedef const boost::numeric::ublas::hermitian_matrix<
  //     transfer_functions::complex, 
  //     boost::numeric::ublas::lower>& access_hermetian_type;
  // const_cast< access_hermetian_type >(n_matrix)(0, 0)

  std::vector< complex_type > tipper_value;
  gmm::dense_matrix< complex_type > impedance_value(2, 2);
  gmm::dense_matrix< complex_type > s_matrix(2, 2);
  gmm::dense_matrix< complex_type > n_matrix(2, 2);
  if(number_of_channels==5)
  {
    gmm::resize(n_matrix,3, 3);
  }
  else if (number_of_channels!=4)
  {
    std::cerr << "importação do formato egbert não está preparado para "
	      << number_of_channels
	      << " canais\n";
    exit(1);
  }

  while ( getline(egbert_file, s) )
  {
    if (s.find("period :") < s.length()) 
    {
      period = atof( string(s, 8, 17).c_str() );
    }
    else if (s.find("Transfer Functions") < s.length())
    {
      if(number_of_channels==5)
      {
	import_tipper_vector_value(tipper_value, egbert_file);
      }
      import_impedance_tensor_value(impedance_value, egbert_file);
    }
    else if (s.find("Inverse Coherent Signal Power Matrix") < s.length())
    {
      import_inverse_coherent_signal_power_matrix(s_matrix, egbert_file);
    }
    else if (s.find("Residual Covar") < s.length())
    {
      import_residual_covariance_matrix(n_matrix, egbert_file);
      break;
    }
  }

  if(number_of_channels==5)
  {
    build_tipper( tipper, tipper_value, s_matrix, n_matrix );
  }
  build_impedance( impedance, impedance_value, s_matrix, n_matrix );
} // import_transfer_function

} // namespace import_egbert
