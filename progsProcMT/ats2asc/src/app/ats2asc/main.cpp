/*! \mainpage

  \section intro_sec Introdução

  Converte os arquivos ats em arquivos asc, sp e clk.

  Os arquivos \p asc e \p clk serão criados no sub-diretório \p DATA
  do diretório atual e os arquivos \p sp serão criados no
  sub-diretório \p SP. Ampos os diretórios serão criados casos não
  existam.

  \b Obs. O nome dos arquivos \p ats descrevem os agrupamentos e os
  canais. Como essas informações não constam do cabeçalhos dos
  arquivos os seus nomes originais não devem ser alterados.

  \section sinopse_sec Sinopse:
  <tt> ats2asc [Opções] DIR </tt>

  - \p DIR é o diretório que contém os arquivos ats.
  - \p Opções
  - \p --help exibe a mensagem de ajuda
  - \p --site-name \e nome prefixo para os arquivos de saída
  - \p --model-file \e file nome do arquivo de modelagem.

  \section config_sec Configuração

  Se existir, o arquivo <tt> .ats2asc </tt> no diretório raíz do
  usuário será usado como arquivo de configuração.

  A única opção permitida nesse arquivo é \p model-file que deve estar
  em uma linha da forma:

  \code
  model-file = caminho_completo_para_o_arquivo_de_modelagem
  \endcode

  \b Obs. A opção \p model-file da linha de comando sobrescreve a do
  arquivo de configuração.

  \section est_asc_sec Estrutura dos arquivos asc

  Os arquivos asc são arquivos ascii que contém os dados de um
  agrupamento separados em colunas

  A ordem das 5 colunas iniciais é Hx Hy Hz Ex Ey, a partir dessas
  colunas não há nenhuma ordem definida */

#include <iostream>      // std::cout
#include "ats2asc.hpp"

int main(int argc, char* argv[]) {
  using namespace ats2asc;

  try
  {
    po::variables_map vm;

    parseOptions(argc, argv, vm);
    if (vm.count("help") || !vm.count("ats-dir"))
    {
      helpMessage(argv[0]);
      return(1);
    }

    checkOptions(vm);

    registra(vm["ats-dir"].as<std::string>(),
	     vm["site-name"].as<std::string>(),
	     vm["model-file"].as<std::string>());
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

} //end main()
