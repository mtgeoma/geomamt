/*!
  \file   ats2asc.hpp
  \author marcos banik
  \date   Sat Aug 12 12:04:22 2006

  \brief

*/

#ifndef E87D5DD8E244D4D8F6AD4ADCF1546DC1
#define E87D5DD8E244D4D8F6AD4ADCF1546DC1

#include <set>
#include <string>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include <boost/ptr_container/ptr_vector.hpp> //ptr_vector std::vector
#include <boost/program_options.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

//! Elementos necessários para a conversão dos arquivos ats
namespace ats2asc
{
//! Functor para a criação dos conjuntos \p ats_set
/*!
  \sa ats_set
*/
struct compare_ats_filenames
  : public std::binary_function<fs::path, fs::path, bool>
{
  //! compara 2 arquivos ats com relação à ordem das colunas no arquivo //asc
  /*!
    \param path_1
    \param path_2

    \return \p true se os dados de \p path_1 devem aparecer em uma coluna
    anterior aos dados de \p path_2

    \sa compare_ats_channels, \ref est_asc_sec
  */
  bool operator() (const fs::path& path_1, const fs::path& path_2) const;
};

//! Typedef do conjunto que contém os nomes-dos-arquivos \p ats
typedef std::set<fs::path, compare_ats_filenames> ats_set;


//---------------

//! Compara 2 arquivos ats com relação à ordem das colunas no arquivo asc.
/*! Compara os nomes dos arquivos apenas com relação aos canais. A função \ref
 *  compare_ats_filenames compara o nome dos arquivos com relação ao número do
 *  equipamento, canal e rodada.

\param file_1
\param file_2

\return \p true se os dados de \p file_1 devem aparecer em uma coluna anterior
aos dados de \p file_2

\sa compare_ats_filenames, \ref est_asc_sec
*/
bool compare_ats_channels(const std::string& file_1, const std::string& file_2);




//! Preenche \p arquivos com os nomes dos arquivos \p ats em \p diretorio
/*\param arquivos conjunto onde os nomes-de-arquivos serão armazenados
\param diretorio diretório que contém os arquivos \p ats
*/
void preenche(ats_set& arquivos, const fs::path& diretorio);



struct nao_e_do_mesmo_grupo
  : public std::binary_function<fs::path, fs::path, bool>
{
  bool operator()(const fs::path& path_1, const fs::path& path_2) const;
};

void checkOptions(po::variables_map const& vm);

bool le_ats( boost::ptr_vector< fs::ifstream >& arq,
	     std::vector< int >& dados );

void pula_cabecalho( boost::ptr_vector< fs::ifstream >& arq_ats );

void registra(std::string const& ats_dir, std::string const& site_name,
	      std::string const& model_file);

void registra_asc( ats_set::iterator inicio, ats_set::iterator fim,
		   const std::string& prefix, const fs::path& data_dir );

void registra_clk( ats_set::iterator inicio, const std::string& prefix,
		   const fs::path& data_dir );

void registra_sp(ats_set::iterator inicio, ats_set::iterator fim,
		 const std::string& model_file, const std::string& prefix,
		 const fs::path& data_dir );

std::string nome_do_grupo( const fs::path& path );

void le_linha_de_comando(po::variables_map& vm,
			 po::options_description& generic,
			 int argc, char* argv[]);

//! Cabeçalho do arquivo ats.
/*!  Nem todos os atributos de um arquivo ats constam da struct. Apenas os que
  são importantes para a conversão para sp.
*/

struct ats_header {
  ats_header(const fs::path& path);
  float sampleFreq;
  int dateTime;
  double lsbval;
  char channelType[2];
  char sensorType[6];
  short sensorSerNum;
  float eFieldDipoleLength;
  float angle;
  int latitude;		/*!< latitude */
  int longitude;	/*!< longitude */
  int elevation;
};

po::options_description genericOptions();

po::options_description configurationOptions();

po::options_description positionalOptions();

po::options_description commandLineOptions();

void helpMessage(char* program);

po::positional_options_description positionalOrder();

void parseOptions(int const argc, char* argv[], po::variables_map& vm);

} // end namespace ats2asc

#endif
