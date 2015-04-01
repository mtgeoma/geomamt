#ifndef TRANSFER_FUNCTIONS_HPP
#define TRANSFER_FUNCTIONS_HPP

#include "udouble/measure.hpp"

#include <gmm/gmm.h> //namespace gmm
#include <utility> //std::pair
#include <vector>

//! \brief estrutura e funções para manipular as funções de transferência
/*!  \todo verificar outras bibliotecas para matriz como
  <a href="http://www.boost.org/libs/numeric/ublas/doc/overview.html">uBLAS</a>,
  <a href="http://home.gna.org/getfem/gmm_intro.html">GMM++</a>,  
  <a href="http://flens.sf.net/">FLENS</a>,
  <a href="http://www.lsc.nd.edu/research/mtl/">mtl</a>,
  <a href="http://www.oonumerics.org/blitz/">blitz++</a>,
  <a href="http://tvmet.sourceforge.net">tvmet</a>, ou
  <a href="http://acts.nersc.gov/pooma/">POOMA</a>
*/
namespace transfer_function
{
/*@{*/ 
//!\brief Tipos usados ao longo do programa

typedef long double real_type;
typedef std::complex<real_type> complex_type;
typedef measurement real_measurement_type;
typedef gmm::dense_matrix< complex_type > complex_matrix_type;
typedef std::complex< real_measurement_type > complex_measurement_type;
typedef gmm::dense_matrix< complex_measurement_type > impedance_tensor_type;

typedef std::vector< complex_measurement_type > tipper_vector_type;
/*@}*/

//! \brief A impedância
/*! A estrutura associa um tensor de impedância a um período
  
 */
struct impedance
{
  real_type period;	/*!< Período */
  impedance_tensor_type Z; /*!< Tensor de impedância */
	
  //! @name Construtores

  //! Construtor
  /*!
    \param period_ período do tensor de impedância
    \param Z_ tensor de impedância

    \return Um objeto do tipo impedance
  */
  impedance(real_type period_ = 0.0,
	    impedance_tensor_type Z_ = impedance_tensor_type(2, 2)) :
    period(period_), Z(Z_) {};
	
  //! @name Componentes do Tensor de Impedância
  //! Retorna a componente de mesmo nome do tensor de impedância
  /*@{*/
  complex_measurement_type Zxx() const {return Z(0, 0);}
  complex_measurement_type Zxy() const {return Z(0, 1);}
  complex_measurement_type Zyx() const {return Z(1, 0);}
  complex_measurement_type Zyy() const {return Z(1, 1);}
  /*@}*/

  //!
  /*! \brief gira a matriz
  
    A matriz é girada em torno do eixo Z no sentido Norte-Leste
    \param angle angulo de giro em graus

    \todo verificar se a matriz precisar ser rotacionada ou não com
boost/fp_comparison.hpp que pode ser encontrado em <a
href="http://www.boost-consulting.com/vault/">Boost Vault/Math -
Numerics</a>,

  */
  void rotate(real_type angle);
};
    
//! Conjunto de todas as impedâncias
typedef std::vector< impedance > impedance_collection_type;

} //namespace transfer_function

#endif
