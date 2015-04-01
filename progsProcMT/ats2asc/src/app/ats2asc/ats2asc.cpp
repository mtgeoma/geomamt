#include "ats2asc.hpp"
#include "declination.h"
#include "system_dependent.hpp"
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/ptr_container/ptr_vector.hpp> //ptr_vector std::vector
#include <cmath> // atan2 angle ats_header
#include <iomanip> //setw setprecision
#include <iostream>

namespace ats2asc
{
bool
compare_ats_filenames::operator()(const fs::path& path_1,
				  const fs::path& path_2) const
{
  std::string file_1_numero_do_equipamento;
  std::string file_1_run;
  std::string file_1_canal;
  std::string file_1_banda;
  std::string file_2_numero_do_equipamento;
  std::string file_2_run;
  std::string file_2_canal;
  std::string file_2_banda;
  std::string const file_1_name = path_1.filename().string();
  std::string const file_2_name = path_2.filename().string();

  if (file_1_name.size() == 12 && file_2_name.size() == 12)
  {
    file_1_numero_do_equipamento = file_1_name.substr(0, 3);
    file_1_run = file_1_name.substr(4, 2);
    file_1_canal = file_1_name.substr(6, 1);
    file_1_banda = file_1_name.substr(7, 1);
    file_2_numero_do_equipamento = file_2_name.substr(0, 3);
    file_2_run = file_2_name.substr(4, 2);
    file_2_canal = file_2_name.substr(6, 1);
    file_2_banda = file_2_name.substr(7, 1);
  }
  else
  {
    file_1_numero_do_equipamento = file_1_name.substr(0, 3);
    file_1_run = file_1_name.substr(13, 3);
    file_1_canal = file_1_name.substr(18, 2);
    file_1_banda = file_1_name.substr(24, file_1_name.size() - 24 - 4);
    file_2_numero_do_equipamento = file_2_name.substr(0, 3);
    file_2_run = file_2_name.substr(13, 3);
    file_2_canal = file_2_name.substr(18, 2);
    file_2_banda = file_2_name.substr(24, file_2_name.size() - 24 - 4);
  }

  boost::algorithm::to_lower(file_1_canal);
  boost::algorithm::to_lower(file_2_canal);

  if (file_1_numero_do_equipamento == file_2_numero_do_equipamento)
  {
    if (file_1_run == file_2_run)
    {
      if (file_1_banda == file_2_banda)
      {
	return compare_ats_channels(file_1_canal, file_2_canal);
      }
      else
      {
	return (file_1_banda < file_2_banda);
      }
    }
    else
    {
      return (file_1_run < file_2_run);
    }
  }
  else
  {
    return (file_1_numero_do_equipamento < file_2_numero_do_equipamento);
  }
} //ats2asc::compare_ats_filenames::operator()

void
checkOptions(po::variables_map const& vm)
{
  //string para testar a existencia dos arquivos e diretórios
  std::string test;

  //condições para ats-dir
  test = vm["ats-dir"].as<std::string>();
  if (!fs::is_directory(test))
  {
    throw std::runtime_error("ats-dir not found: " + test);
  }

  //condições para model-file
  if (!vm.count("model-file"))
  {
    throw std::runtime_error("model-file not found");
  }

  test = vm["model-file"].as<std::string>();
  if (!fs::exists(test))
  {
    throw std::runtime_error("model-file not found: " + test);
  }
}

bool
compare_ats_channels(const std::string& canal_1, const std::string& canal_2)
{
  //inicializa com um valor inválido que será conferido antes do retorno
  int file_1_col = -1, file_2_col = -1; //colunas das variáveis x e y

  if (canal_1 == "a" || canal_1 == "ex")
  {
    file_1_col = 3;
  }
  else if (canal_1 == "b" || canal_1 == "ey")
  {
    file_1_col = 4;
  }
  else if (canal_1 == "x" || canal_1 == "hx")
  {
    file_1_col = 0;
  }
  else if (canal_1 == "y" || canal_1 == "hy")
  {
    file_1_col = 1;
  }
  else if (canal_1 == "z" || canal_1 == "hz")
  {
    file_1_col = 2;
  }
  if (file_1_col == -1)
  {
    std::cerr << "unknow channel 1" << canal_1 << "\n";
    exit(1);
  }

  if (canal_2 == "a" || canal_2 == "ex")
  {
    file_2_col = 3;
  }
  else if (canal_2 == "b" || canal_2 == "ey")
  {
    file_2_col = 4;
  }
  else if (canal_2 == "x" || canal_2 == "hx")
  {
    file_2_col = 0;
  }
  else if (canal_2 == "y" || canal_2 == "hy")
  {
    file_2_col = 1;
  }
  else if (canal_2 == "z" || canal_2 == "hz")
  {
    file_2_col = 2;
  }
  if (file_2_col == -1)
  {
    std::cerr << "unknow channel 2" << canal_2 << "\n";
    exit(1);
  }

  return (file_1_col < file_2_col);

} //ats2asc::compare_ats_channels

void
preenche(ats_set& arquivos, const fs::path& diretorio)
{
  fs::directory_iterator end_itr;
  for (fs::directory_iterator itr(diretorio); itr != end_itr; ++itr)
  {
    if (itr->path().extension() == ".ats")
    {
      arquivos.insert(*itr);
    }
  }
} //preenche(std::set<fs::path>& arquivos, fs::path diretorio);

bool
nao_e_do_mesmo_grupo::operator()(const fs::path& path_1,
				 const fs::path& path_2) const
{
  std::string arq_1_numero_do_equipamento;
  std::string arq_1_run;
  std::string arq_1_banda;
  std::string arq_2_numero_do_equipamento;
  std::string arq_2_run;
  std::string arq_2_banda;
  std::string const file_1_name = path_1.filename().string();
  std::string const file_2_name = path_2.filename().string();

  if (file_1_name.size() == 12 && file_2_name.size() == 12)
  {
    arq_1_numero_do_equipamento = file_1_name.substr(0, 3);
    arq_1_run = file_1_name.substr(4, 2);
    arq_1_banda = file_1_name.substr(7, 1);
    arq_2_numero_do_equipamento = file_2_name.substr(0, 3);
    arq_2_run = file_2_name.substr(4, 2);
    arq_2_banda = file_2_name.substr(7, 1);
  }
  else
  {
    arq_1_numero_do_equipamento = file_1_name.substr(0, 3);
    arq_1_run = file_1_name.substr(13, 3);
    arq_1_banda = file_1_name.substr(24, file_1_name.size() - 24 - 4);
    arq_2_numero_do_equipamento = file_2_name.substr(0, 3);
    arq_2_run = file_2_name.substr(13, 3);
    arq_2_banda = file_2_name.substr(24, file_2_name.size() - 24 - 4);
  }

  return !((arq_1_numero_do_equipamento == arq_2_numero_do_equipamento) &&
	   (arq_1_run == arq_2_run) &&
	   (arq_1_banda == arq_2_banda));
}

bool
le_ats(boost::ptr_vector<fs::ifstream>& arq, std::vector<int>& data)
{
  bool r = true;
  std::vector<int>::size_type size = data.size();

  for(std::vector<int>::size_type i = 0; i != size; ++i)
  {
    r = r && arq[i].read(reinterpret_cast<char *>(&data[i]), sizeof(data[i]));
  } /* end of for(unsigned int i = 0; i != size; ++i) */

  return r;
}

void
pula_cabecalho(boost::ptr_vector<fs::ifstream>& arq_ats)
{
  short header_length;

  for(boost::ptr_vector<fs::ifstream>::iterator i = arq_ats.begin();
      i != arq_ats.end();
      ++i)
  {
    i->read(reinterpret_cast<char *>(&header_length), sizeof(header_length));
    i->ignore(static_cast<int>(header_length - i->tellg()));
  } /* end of for(boost::ptr_vector<fs::ifstream>::iterator i =
       arq_ats.begin(); i != end; ++i) */
}

void
registra_asc(ats_set::iterator inicio, ats_set::iterator fim,
	     const std::string& prefix, const fs::path& data_dir)
{
  //container com os arquivos ats
  boost::ptr_vector<fs::ifstream> arquivos_ats;
  for(ats_set::iterator itr = inicio; itr != fim; ++itr)
  {
    arquivos_ats.push_back(new fs::ifstream(*itr, std::ios_base::binary));
  } /* end of for(itr = inicio; itr != fim; ++itr) */

  //arquivo asc
  std::string nome_arq_asc = prefix + nome_do_grupo(*inicio) + ".asc";
  fs::ofstream arquivo_asc(data_dir / nome_arq_asc);

  //dados
  std::vector<int> dados(arquivos_ats.size());
  pula_cabecalho(arquivos_ats);
  while(le_ats(arquivos_ats, dados))
  {
    for(std::vector<int>::iterator i = dados.begin(); i != dados.end(); ++i)
    {
      arquivo_asc <<std::setw(12) <<*i;
    } /* end of for(std::vector<int>::iterator i = dados.begin();
	 i != end; ++i) */

    arquivo_asc << "\n";
  } /* end of while(le_ats(arquivos_ats, dados)) */
}

std::string
nome_do_grupo(const fs::path& path)
{
  std::string nome;
  std::string const filename = path.filename().string();
  if (filename.size() == 12)
  {
    nome = filename.substr(0, 3) + "_" + filename.substr(4, 2) +
      filename.substr(7, 1);
  }
  else
  {
    nome = filename.substr(0, 3) + "_" +
      filename.substr(13, 3) + "_" +
      filename.substr(24, filename.size() - 24 - 4);
  }

  return nome;
}

void
registra(std::string const& ats_dir, std::string const& site_name,
	 std::string const& model_file)
{
  using namespace ats2asc;
  namespace fs = boost::filesystem;

  //! Preenche o \p set \p arquivos_ats com os arquivos ats em ats_dir
  ats_set arquivos_ats;
  preenche(arquivos_ats, ats_dir);

  ats_set::iterator inicio_do_grupo = arquivos_ats.begin();
  ats_set::iterator fim_do_grupo = arquivos_ats.begin();
  ats_set::iterator fim_dos_arquivos = arquivos_ats.end();

  fs::create_directory( "DATA" );
  fs::path asc_dir = "DATA";

  fs::create_directory( "SP" );
  fs::path sp_dir = "SP";

  while(inicio_do_grupo != fim_dos_arquivos)
  {
    fim_do_grupo = std::find_if(inicio_do_grupo, fim_dos_arquivos,
				std::bind1st( nao_e_do_mesmo_grupo(),
					      *inicio_do_grupo ) );

    registra_asc(inicio_do_grupo, fim_do_grupo, site_name, asc_dir);
    registra_clk(inicio_do_grupo, site_name, asc_dir);
    registra_sp(inicio_do_grupo, fim_do_grupo, model_file, site_name, sp_dir);
    inicio_do_grupo = fim_do_grupo;
  } /* end of while(inicio_do_grupo != fim_dos_arquivos) */
}

void
registra_clk(ats_set::iterator inicio, const std::string& prefix,
	     const fs::path& data_dir)
{
  fs::ifstream arq_ats(*inicio, std::ios_base::binary);

  std::string nome_arq_clk = prefix + nome_do_grupo(*inicio) + ".clk";
  fs::ofstream arq_clk(data_dir / nome_arq_clk);

  float sample_frequency;
  int   start_time;
  arq_ats.ignore(8);
  arq_ats.read(reinterpret_cast<char *>(&sample_frequency),
	       sizeof(sample_frequency));
  arq_ats.read(reinterpret_cast<char *>(&start_time), sizeof(start_time));
  arq_clk.setf(std::ios::scientific);
  arq_clk << std::setprecision(10) << (1.0f / sample_frequency) << "\n";

  const time_t seconds(start_time);
  struct tm* ptm = gmtime(&seconds);

  //a mesma linha 2 vezes inicio
  arq_clk << " " << std::setw(2) << std::setfill('0')
	  << (ptm->tm_year % 100)
	  << " "  << std::setw(2) << std::setfill('0')
	  << (ptm->tm_mon + 1)
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_mday
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_hour
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_min
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_sec << "\n";

  arq_clk << " " << std::setw(2) << std::setfill('0')
	  << (ptm->tm_year % 100)
	  << " "  << std::setw(2) << std::setfill('0')
	  << (ptm->tm_mon + 1)
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_mday
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_hour
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_min
	  << " "  << std::setw(2) << std::setfill('0')
	  << ptm->tm_sec << "\n";
  //a mesma linha 2 vezes fim
}

void
registra_sp(ats_set::iterator inicio, ats_set::iterator fim,
	    const std::string& model_file, const std::string& prefix,
	    const fs::path& data_dir)
{
  boost::ptr_vector<ats_header> headers;
  unsigned numero_de_canais = 0;
  std::map<std::string, unsigned int> mapa_canais;
  for(std::set<fs::path, compare_ats_filenames>::iterator i = inicio;
      i != fim;
      ++i)
  {
    headers.push_back(new ats_header(*i));
    mapa_canais[std::string(headers[numero_de_canais].channelType,2)] =
      numero_de_canais;
    ++numero_de_canais;
  } /* end of for(std::set<fs::path, compare_ats_filenames>::iterator
       i = inicio; i != fim; ++i) */

  if (mapa_canais.find("Hx") != mapa_canais.end())
  {
    if (mapa_canais.find("Ex") != mapa_canais.end())
    {
      headers[mapa_canais["Hx"]].angle=headers[mapa_canais["Ex"]].angle;
    }
  }
  if (mapa_canais.find("Hy") != mapa_canais.end())
  {
    if (mapa_canais.find("Ey") != mapa_canais.end())
    {
      headers[mapa_canais["Hy"]].angle=headers[mapa_canais["Ey"]].angle;
    }
  }

  std::string nome_do_arquivo_sp = prefix + nome_do_grupo(*inicio) + ".sp";
  fs::ofstream arquivo_sp(data_dir / nome_do_arquivo_sp);
  arquivo_sp << prefix + nome_do_grupo(*inicio) << "\n";

  double latitude_deg = headers[0].latitude / 3600000.0;
  double longitude_deg = headers[0].longitude / 3600000.0;

  const time_t seconds(headers[0].dateTime);
  struct tm* ptm = gmtime(&seconds);

  double declina = declination(latitude_deg, longitude_deg,
			       (headers[0].elevation / 100000.0),
			       (1900 + ptm->tm_year), (ptm->tm_mon + 1),
			       ptm->tm_mday, model_file);

  arquivo_sp << latitude_deg << " " << longitude_deg << "\n";

  const std::ios_base::fmtflags old_opt = arquivo_sp.flags();
  int prec = static_cast<int>(arquivo_sp.precision());
  arquivo_sp << std::fixed << std::setprecision(0) << declina
	     << std::setprecision(prec) << "\n";
  arquivo_sp.flags(old_opt);

  arquivo_sp << numero_de_canais << "\n";
  arquivo_sp.setf(std::ios::scientific);
  arquivo_sp << std::setprecision(10) << (1.0f / headers[0].sampleFreq) << "\n";
  arquivo_sp.flags(old_opt);

  arquivo_sp << "0. 0." << "\n";

  for(unsigned i = 0; i != numero_de_canais; ++i)
  {
    if (headers[i].channelType[0] == 'H')
    {
      std::string sensorType = std::string(headers[i].sensorType,6);
      boost::algorithm::trim_if(sensorType,boost::algorithm::is_cntrl());
      boost::algorithm::to_upper(sensorType);
      std::string cal_parameters;
      if ( sensorType == "MFS07")
      {
	cal_parameters = "0.64 32.0 20.0e3";
      }
      else
      {
	cal_parameters = "0.80 4.00 8192.0";
      }

      arquivo_sp << headers[i].channelType[0]
		 << headers[i].channelType[1] << "\n";
      arquivo_sp << headers[i].angle << " " << "0.\n";
      arquivo_sp.setf(std::ios::scientific);
      arquivo_sp << std::setprecision(5) << headers[i].lsbval << " "
		 << "1\n";
      arquivo_sp.flags(old_opt);
      arquivo_sp << "MS\n";
      arquivo_sp << sensorType
		 << std::setw(3) << std::setfill('0')
		 << headers[i].sensorSerNum;
      if(headers[i].sampleFreq<=512.0f)
      {
	arquivo_sp << ".TXT on " << cal_parameters
		   << " # change 'on' to 'ttf' to use theoretical "
		   << "MFS transfer function\n";
      }
      else
      {
	arquivo_sp << ".TXT off " << cal_parameters
		   << " # change 'off' to 'ttf' to use theoretical "
		   << "MFS transfer function\n";
      }
    }
    else if (headers[i].channelType[0] == 'E')
    {
      arquivo_sp << headers[i].channelType[0]
		 << headers[i].channelType[1] << "\n";
      arquivo_sp.setf(std::ios::fixed);
      arquivo_sp << std::setprecision(4)
		 << (headers[i].eFieldDipoleLength / 1000.0f)
		 << " " << headers[i].angle << " 0. 1.\n";
      arquivo_sp.setf(std::ios::scientific);
      arquivo_sp << std::setprecision(6) << (-1.0 * headers[i].lsbval)
		 << " 0\n";
      arquivo_sp.flags(old_opt);
    }
  } /* end of for(unsigned i = 0; i != numero_de_canais; ++i) */
}

ats_header::ats_header(const fs::path& path)
{
  fs::ifstream arquivo_ats(path);
  arquivo_ats.seekg(0x8);
  arquivo_ats.read(reinterpret_cast<char *>(&sampleFreq), sizeof(sampleFreq));

  arquivo_ats.read(reinterpret_cast<char *>(&dateTime), sizeof(dateTime));

  arquivo_ats.read(reinterpret_cast<char *>(&lsbval), sizeof(lsbval));

  arquivo_ats.seekg(0x26);
  arquivo_ats.read(reinterpret_cast<char *>(&channelType), sizeof(channelType));

  arquivo_ats.seekg(0x28);
  arquivo_ats.read(reinterpret_cast<char *>(&sensorType), sizeof(sensorType));

  arquivo_ats.read(reinterpret_cast<char *>(&sensorSerNum),
		   sizeof(sensorSerNum));

  float x1, y1, z1, x2, y2, z2;
  arquivo_ats.seekg(0x30);
  arquivo_ats.read(reinterpret_cast<char *>(&x1), sizeof(x1));
  arquivo_ats.read(reinterpret_cast<char *>(&y1), sizeof(y1));
  arquivo_ats.read(reinterpret_cast<char *>(&z1), sizeof(z1));
  arquivo_ats.read(reinterpret_cast<char *>(&x2), sizeof(x2));
  arquivo_ats.read(reinterpret_cast<char *>(&y2), sizeof(y2));
  arquivo_ats.read(reinterpret_cast<char *>(&z2), sizeof(z2));
  eFieldDipoleLength = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
				 + (z2-z1)*(z2-z1));

  float const pi = boost::math::constants::pi<float>();
  angle = std::atan2(y2-y1,x2-x1)*(180.0f/pi);

  arquivo_ats.seekg(0x60);
  arquivo_ats.read(reinterpret_cast<char *>(&latitude), sizeof(latitude));

  arquivo_ats.read(reinterpret_cast<char *>(&longitude), sizeof(longitude));

  arquivo_ats.read(reinterpret_cast<char *>(&elevation), sizeof(elevation));
}

//! Opções usadas exclusivamente na linha de comando
/*!
*/
po::options_description
genericOptions()
{
  po::options_description opt("Command line");
  opt.add_options()
    ("help", "show this message")
    ("site-name", po::value<std::string>()->default_value(""),
     "label to site name");
  return opt;
}

//! Opções usadas no arquivo de configuração
/*!
 */
po::options_description
configurationOptions()
{
  po::options_description opt("Configuration");
  opt.add_options()("model-file", po::value<std::string>(),
		    "IGRF11 model file");
  return opt;
}

//! Opções Posicionadas
/*! Opções que dependem da posição da linha de comando,
  os seus nomes não deveriam constar da mensagem de ajuda.
*/
po::options_description
positionalOptions()
{
  po::options_description opt;
  opt.add_options()("ats-dir", po::value<std::string>(), "ats dir");
  return opt;
}

po::options_description
commandLineOptions()
{
  po::options_description generic = genericOptions();
  po::options_description config = configurationOptions();
  po::options_description pos = positionalOptions();

  po::options_description cmd;
  cmd.add(generic).add(config).add(pos);
  return cmd;
}

//! Texto que descreve o uso das opções
/*!
 */
void
helpMessage(char* program)
{
  po::options_description generic = genericOptions();
  po::options_description config = configurationOptions();
  po::options_description visible("Options");
  visible.add(generic).add(config);

  std::cout << program << " convert ats files to EMTF "
	    << "input files: asc clk and sp\n"
	    << "Usage: "
	    << program << " [Options] ats-dir\n"
	    << "  ats-dir\t\tdirectory with ats files\n\n"
	    << visible << "\n";
}

//! //relaciona o nome da opção à sua posição na linha de comando
/*!
*/
po::positional_options_description
positionalOrder()
{
  po::positional_options_description opt;
  opt.add("ats-dir", -1);

  return opt;
}

// analisa as opções da linha de comando e do arquivo de configuração
void
parseOptions(int const argc, char* argv[], po::variables_map& vm)
{
  po::options_description cmdline_options = commandLineOptions();
  po::positional_options_description p = positionalOrder();

  store(po::command_line_parser(argc, argv).options(cmdline_options)
	.positional(p).run(), vm);

  fs::path teste(diretorio_do_usuario() + "/.ats2asc");
  fs::ifstream ifs(teste);

  po::options_description config = configurationOptions();
  store(parse_config_file(ifs, config), vm);
  notify(vm);
}

}
