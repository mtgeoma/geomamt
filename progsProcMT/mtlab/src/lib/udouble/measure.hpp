#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <cmath> // sqrt

//!
/*! \brief Representação de uma medida (valor +/- incerteza)

 */
struct measurement
{
  long double value;		/*!< valor da medida */
  long double variance;	/*!< variança da medida */

  //! @name Construtor
  /*@{*/
  //!
  /*! \brief Construtor default

    \param val valor da medida
    \param var variança da medida

    \return um objeto do tipo \a measurement
  */
  measurement ( long double val = 0.0, long double var = 0.0 )
    : value(val), variance(var) {};
  /*@}*/

  //!
  /*! \brief Calcula o valor da incereteza
    \return O valor da Incerteza
  */
  long double uncertainty() { return sqrt( variance ); }
};

#endif
