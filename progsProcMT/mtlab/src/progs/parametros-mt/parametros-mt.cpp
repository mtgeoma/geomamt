#include "parametros-mt.hpp"
#include <cmath> //M_PIl
#include <iostream> //cout

void
program::command_line(
    boost::program_options::variables_map& vm,
    boost::program_options::options_description& visible,
    int argc, char *argv[] )
{
    namespace po = boost::program_options;

    // Declara as opções genericas
    visible.add_options()
	("help", "exibe essa ajuda")
	("azimute", po::value<long double>()->default_value(0.0), "ângulo de rotação do tensor em relação ao sistema de coordenada da medida");

    // Declara as opções posicionais. As opções que dependem apenas da posição
    //devem ficar ocultas do usuário.
    po::options_description hidden("Opções Ocultas");
    hidden.add_options()
	("file", po::value< std::string >(), "")
	("parameter", po::value< std::string >(), "");

    // Estabelece a ordem das opções posicionais.
    po::positional_options_description pos_order;
    pos_order.add("file", 1);
    pos_order.add("parameter", 1);

    // As opções de linha de comando ( ocultas + genéricas )
    po::options_description cmdline_options;
    cmdline_options.add(visible).add(hidden);

    //Descrição da linha de comando para o programa principal

    // interpreta as opções
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos_order).run(), vm);

}// command_line_ok

void
program::print_help( std::string program_name,
		     boost::program_options::options_description& visible )
{
    std::cerr << "Uso: " << program_name << " ARQUIVO PARAMETRO-COMPONENTE\n";
    std::cerr << "ARQUIVO é o arquivo com os tensores de impedância\n"
	      << "PARAMETRO: rho, phi, ZReal, ZImag, rhoReal, rhoImag\n"
	      << "COMPONENTE: xx, xy, yx, yy, mean, eff[ective]\n";
    std::cerr << visible;
}

void program::lista_parametro(
    std::string const& par,
    long double azimute,
    transfer_function::impedance_collection_type& imp )
{
    namespace tf   = transfer_function;
    long double pi = M_PIl;
    long double mu = 4.0e-7 * pi;

    for( tf::impedance_collection_type::size_type i = 0;
	 i != imp.size(); ++i )
    {
	imp[i].rotate(azimute);
	if ( par == "ZReal-xx" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zxx().real().value;
	    long double uncertainty = imp[i].Zxx().real().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZReal-xy" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zxy().real().value;
	    long double uncertainty = imp[i].Zxy().real().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZReal-yx" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zyx().real().value;
	    long double uncertainty = imp[i].Zyx().real().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZReal-yy" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zyy().real().value;
	    long double uncertainty = imp[i].Zyy().real().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZImag-xx" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zxx().imag().value;
	    long double uncertainty = imp[i].Zxx().imag().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZImag-xy" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zxy().imag().value;
	    long double uncertainty = imp[i].Zxy().imag().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZImag-yx" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zyx().imag().value;
	    long double uncertainty = imp[i].Zyx().imag().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "ZImag-yy" )
	{
	    long double      period = imp[i].period;
	    long double       value = imp[i].Zyy().imag().value;
	    long double uncertainty = imp[i].Zyy().imag().uncertainty();
 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoReal-xx" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zxx().real().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zxx().real().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoReal-xy" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zxy().real().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zxy().real().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoReal-yx" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zyx().real().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zyx().real().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoReal-yy" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zyy().real().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zyy().real().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoImag-xx" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zxx().imag().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zxx().imag().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoImag-xy" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zxy().imag().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zxy().imag().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoImag-yx" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zyx().imag().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zyx().imag().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rhoImag-yy" )
	{
	    long double      period = imp[i].period;
	    long double       value =
		imp[i].Zyy().imag().value / sqrt(2.0 * pi * mu / period);
	    long double uncertainty =
		imp[i].Zyy().imag().uncertainty() / sqrt(2.0 * pi * mu
							 /
							 period);

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rho-xx" )
	{
	    long double period = imp[i].period;
	    long double Zxx_r  = imp[i].Zxx().real().value;
	    long double Zxx_i  = imp[i].Zxx().imag().value;
	    long double Zxx_rv = imp[i].Zxx().real().variance;
	    long double Zxx_iv = imp[i].Zxx().imag().variance;
	    long double mod2_Z  = pow(Zxx_r, 2) + pow(Zxx_i, 2);

	    long double       value =
		( mod2_Z ) / (2.0 * pi * mu / period);

	    long double uncertainty =
		2.0 / ( mod2_Z * log( 10.0L ) )
		*
		sqrt( pow(Zxx_r, 2) * Zxx_rv + pow(Zxx_i, 2) * Zxx_iv );

 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
	}
	else if ( par == "rho-xy" )
	{
	    long double period = imp[i].period;
	    long double Zxy_r  = imp[i].Zxy().real().value;
	    long double Zxy_i  = imp[i].Zxy().imag().value;
	    long double Zxy_rv = imp[i].Zxy().real().variance;
	    long double Zxy_iv = imp[i].Zxy().imag().variance;
	    long double mod2_Z  = pow(Zxy_r, 2) + pow(Zxy_i, 2);

	    long double       value =
		( mod2_Z ) / (2.0 * pi * mu / period);

	    long double uncertainty =
		2.0 / ( mod2_Z * log( 10.0L ) )
		*
		sqrt( pow(Zxy_r, 2) * Zxy_rv + pow(Zxy_i, 2) * Zxy_iv );

 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
	}
	else if ( par == "rho-yx" )
	{
	    long double period = imp[i].period;
	    long double Zyx_r  = imp[i].Zyx().real().value;
	    long double Zyx_i  = imp[i].Zyx().imag().value;
	    long double Zyx_rv = imp[i].Zyx().real().variance;
	    long double Zyx_iv = imp[i].Zyx().imag().variance;
	    long double mod2_Z  = pow(Zyx_r, 2) + pow(Zyx_i, 2);

	    long double       value =
		( mod2_Z ) / (2.0 * pi * mu / period);

	    long double uncertainty =
		2.0 / ( mod2_Z * log( 10.0L ) )
		*
		sqrt( pow(Zyx_r, 2) * Zyx_rv + pow(Zyx_i, 2) * Zyx_iv );

 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
	}
	else if ( par == "rho-yy" )
	{
	    long double period = imp[i].period;
	    long double Zyy_r  = imp[i].Zyy().real().value;
	    long double Zyy_i  = imp[i].Zyy().imag().value;
	    long double Zyy_rv = imp[i].Zyy().real().variance;
	    long double Zyy_iv = imp[i].Zyy().imag().variance;
	    long double mod2_Z  = pow(Zyy_r, 2) + pow(Zyy_i, 2);

	    long double       value =
		( mod2_Z ) / (2.0 * pi * mu / period);

	    long double uncertainty =
		2.0 / ( mod2_Z * log( 10.0L ) )
		*
		sqrt( pow(Zyy_r, 2) * Zyy_rv + pow(Zyy_i, 2) * Zyy_iv );

 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
	}
	else if ( par == "phi-xx" )
	{
	    long double period = imp[i].period;
	    long double      y = imp[i].Zxx().imag().value;
	    long double      x = imp[i].Zxx().real().value;
	    long double    y_v = 2.0L * imp[i].Zxx().real().variance;

	    long double value = atan2( y, x ) * 180.0L / pi;

	    long double uncertainty =
		sqrt ( y_v / ( 2.0L * ( pow( x, 2 ) + pow( y, 2 ) ) ) )
		*
		180.0L / pi;

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "phi-xy" )
	{
	    long double period = imp[i].period;
	    long double      y = imp[i].Zxy().imag().value;
	    long double      x = imp[i].Zxy().real().value;
	    long double    y_v = 2.0L * imp[i].Zxy().real().variance;

	    long double value = atan2( y, x ) * 180.0L / pi;

	    long double uncertainty =
		sqrt ( y_v / ( 2.0L * ( pow( x, 2 ) + pow( y, 2 ) ) ) )
		*
		180.0L / pi;

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "phi-yx" )
	{
	    long double period = imp[i].period;
	    long double      y = imp[i].Zyx().imag().value;
	    long double      x = imp[i].Zyx().real().value;
	    long double    y_v = 2.0L * imp[i].Zyx().real().variance;

	    long double value = atan2( y, x ) * 180.0L / pi;

	    long double uncertainty =
		sqrt ( y_v / ( 2.0L * ( pow( x, 2 ) + pow( y, 2 ) ) ) )
		*
		180.0L / pi;

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "phi-yy" )
	{
	    long double period = imp[i].period;
	    long double      y = imp[i].Zyy().imag().value;
	    long double      x = imp[i].Zyy().real().value;
	    long double    y_v = 2.0L * imp[i].Zyy().real().variance;

	    long double value = atan2( y, x ) * 180.0L / pi;

	    long double uncertainty =
		sqrt ( y_v / ( 2.0L * ( pow( x, 2 ) + pow( y, 2 ) ) ) )
		*
		180.0L / pi;

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rho-mean" )
	{
	    tf::real_type period  = imp[i].period;
	    tf::real_type Rmean   =
 		( imp[i].Zxy().real().value - imp[i].Zyx().real().value )
		/
		tf::real_type(2);

	    tf::real_type Imean   =
		( imp[i].Zxy().imag().value - imp[i].Zyx().imag().value )
		/
		tf::real_type(2);

	    long double Rmean_v =
		( imp[i].Zxy().real().variance + imp[i].Zyx().real().variance )
		/
		tf::real_type(4);

	    long double Imean_v =
		( imp[i].Zxy().imag().variance + imp[i].Zyx().imag().variance )
		/
		tf::real_type(4);

	    long double norm_Z  = pow(Rmean, 2) + pow(Imean, 2);

	    long double  value =
		( norm_Z ) / (2.0 * pi * mu / period);

	    long double uncertainty =
		2.0 / ( norm_Z * log( 10.0L ) )
		*
		sqrt( pow(Rmean, 2) * Rmean_v + pow(Imean, 2) * Imean_v );

 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
	}
	else if ( par == "phi-mean" )
	{
	    tf::real_type period = imp[i].period;

	    tf::real_type y =
		( imp[i].Zxy().imag().value - imp[i].Zyx().imag().value )
		/
		tf::real_type(2);

	    tf::real_type x =
		( imp[i].Zxy().real().value - imp[i].Zyx().real().value )
		/
		tf::real_type(2);

	    tf::real_type y_v =
		tf::real_type(2)
		*
		( imp[i].Zxy().imag().variance + imp[i].Zyx().imag().variance );

	    long double value = atan2( y, x ) * 180.0L / pi;

	    long double uncertainty =
		sqrt ( y_v / ( 2.0L * ( pow( x, 2 ) + pow( y, 2 ) ) ) )
		*
		180.0L / pi;

 	    std::cout << period << "\t" << value << "\t" << uncertainty << "\n";
	}
	else if ( par == "rho-eff" || par == "rho-effective"  )
	{
	    tf::real_type period  = imp[i].period;
	    tf::complex_type Zxx = tf::complex_type(imp[i].Zxx().real().value,
						    imp[i].Zxx().imag().value);
	    tf::complex_type Zxy = tf::complex_type(imp[i].Zxy().real().value,
						    imp[i].Zxy().imag().value);
	    tf::complex_type Zyx = tf::complex_type(imp[i].Zyx().real().value,
						    imp[i].Zyx().imag().value);
	    tf::complex_type Zyy = tf::complex_type(imp[i].Zyy().real().value,
						    imp[i].Zyy().imag().value);

	    tf::complex_type rho_eff   =
		sqrt( Zxx * Zyy - Zxy * Zyx );

 	    tf::real_type norm_rho  = norm(rho_eff);

	    tf::real_type value =
		( norm_rho ) / (tf::real_type(2) * pi * mu / period);

	    tf::real_type Rxx = imp[i].Zxx().real().value;
	    tf::real_type Rxy = imp[i].Zxy().real().value;
	    tf::real_type Ryx = imp[i].Zyx().real().value;
	    tf::real_type Ryy = imp[i].Zyy().real().value;
	    tf::real_type Ixx = imp[i].Zxx().imag().value;
	    tf::real_type Ixy = imp[i].Zxy().imag().value;
	    tf::real_type Iyx = imp[i].Zyx().imag().value;
	    tf::real_type Iyy = imp[i].Zyy().imag().value;

	    tf::real_type d_Rxx =
		( tf::real_type(2) * Ryy
		  *
		  ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy *Iyx )
		  +
		  tf::real_type(2) * Iyy
		  *
		  ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx ))
		/
		( log(tf::real_type(10))
		  *
		  ( pow( ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx ), 2 )
		    +
		    pow( ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Rxy	=
		(tf::real_type(-2) * Ryx
		 *
		 (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		 -
		 tf::real_type(2) * Iyx
		 *
		 ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx))
		/
		( log(tf::real_type(10))
		  *
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Ryx =
		(tf::real_type(-2) * Rxy
		 *
		 ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx )
		 -
		 tf::real_type(2) * Ixy
		 *
		 (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx))
		/
		( log(tf::real_type(10))
		  *
		  ( pow( (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow( (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Ryy =
		(tf::real_type(2) * Rxx
		 *
		 ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		 +
		 tf::real_type(2) * Ixx
		 *
		 ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx))
		/
		( log(tf::real_type(10) )
		  *
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Ixx =
		( tf::real_type(2) * Ryy
		  *
		  ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		  -
		  tf::real_type(2) * Iyy
		  *
		  ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx))
		/
		( log(tf::real_type(10))
		  *
		  ( pow( (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow( (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Ixy =
		( tf::real_type(2) * Iyx
		  *
		  ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx )
		  -
		  tf::real_type(2) * Ryx
		  *
		  ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx))
		/
		( log(tf::real_type(10))
		  *
		  ( pow( (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow( (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Iyx =
		( tf::real_type(2) * Ixy
		  *
		  ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  -
		  tf::real_type(2) * Rxy
		  *
		  ( Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx))
		/
		( log(tf::real_type(10))
		  *
		  ( pow( (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow( (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type d_Iyy =
		(tf::real_type(2) * Rxx
		 *
		 (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		 -
		 tf::real_type(2) * Ixx
		 *
		 ( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx))
		/
		( log(tf::real_type(10) )
		  *
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)
		    +
		    pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)));

	    tf::real_type Rxx_v = imp[i].Zxx().real().variance;
	    tf::real_type Rxy_v = imp[i].Zxy().real().variance;
	    tf::real_type Ryx_v = imp[i].Zyx().real().variance;
	    tf::real_type Ryy_v = imp[i].Zyy().real().variance;
	    tf::real_type Ixx_v = imp[i].Zxx().imag().variance;
	    tf::real_type Ixy_v = imp[i].Zxy().imag().variance;
	    tf::real_type Iyx_v = imp[i].Zyx().imag().variance;
	    tf::real_type Iyy_v = imp[i].Zyy().imag().variance;

	    tf::real_type uncertainty =
		sqrt( pow(d_Rxx, 2) * Rxx_v
		      +
		      pow(d_Rxy, 2) * Rxy_v
		      +
		      pow(d_Ryx, 2) * Ryx_v
		      +
		      pow(d_Ryy, 2) * Ryy_v
		      +
		      pow(d_Ixx, 2) * Ixx_v
		      +
		      pow(d_Ixy, 2) * Ixy_v
		      +
		      pow(d_Iyx, 2) * Iyx_v
		      +
		      pow(d_Iyy, 2) * Iyy_v );


 	    std::cout << period << "\t" << value << "\t"
		 << value / pow(10.0L, uncertainty )
		 << "\t" << value << "\t" << value << "\t"
		 << value * pow(10.0L, uncertainty ) << "\n";
    }
	else if ( par == "phi-eff" || par == "phi-effective" )
	{
	    tf::real_type period  = imp[i].period;
	    tf::complex_type Zxx = tf::complex_type(imp[i].Zxx().real().value,
						    imp[i].Zxx().imag().value);
	    tf::complex_type Zxy = tf::complex_type(imp[i].Zxy().real().value,
						    imp[i].Zxy().imag().value);
	    tf::complex_type Zyx = tf::complex_type(imp[i].Zyx().real().value,
						    imp[i].Zyx().imag().value);
	    tf::complex_type Zyy = tf::complex_type(imp[i].Zyy().real().value,
						    imp[i].Zyy().imag().value);

	    tf::complex_type phi_eff = sqrt( Zxx * Zyy - Zxy * Zyx );
	    tf::real_type       Reff = phi_eff.real();
	    tf::real_type       Ieff = phi_eff.imag();
	    tf::real_type      value =
		atan2(Ieff, Reff) * tf::real_type(180) / pi;

	    tf::real_type Rxx = imp[i].Zxx().real().value;
	    tf::real_type Rxy = imp[i].Zxy().real().value;
	    tf::real_type Ryx = imp[i].Zyx().real().value;
	    tf::real_type Ryy = imp[i].Zyy().real().value;
	    tf::real_type Ixx = imp[i].Zxx().imag().value;
	    tf::real_type Ixy = imp[i].Zxy().imag().value;
	    tf::real_type Iyx = imp[i].Zyx().imag().value;
	    tf::real_type Iyy = imp[i].Zyy().imag().value;

	    tf::real_type d_Rxx =
		( Iyy / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx )
		  -
		  ( Ryy * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx )
		    /
		    (pow(( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2 )
		  /
		  ( pow(( Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  +
		  tf::real_type(1));

	    tf::real_type d_Rxy =
		(Ryx * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		 /
		 (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		 -
		 (Iyx / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)))
		/
		(pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		 /
		 (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx),2))
		 +
		 tf::real_type(1));

	    tf::real_type d_Ryx =
		( Rxy * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		  /
		  (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  -
		  (Ixy / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx),2) )
		  +
		  tf::real_type(1) );

	    tf::real_type d_Ryy =
		( Ixx / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  -
		  ( Rxx * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		    /
		    ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  +
		  tf::real_type(1));

	    tf::real_type d_Ixx =
		( Ryy / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  +
		  Iyy * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		  /
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  +
		  tf::real_type(1));

	    tf::real_type d_Ixy =
		( -Ryx / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  -
		  ( Iyx * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		    /
		    ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx),2))))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  +
		  tf::real_type(1));

	    tf::real_type d_Iyx
		=
		( -Rxy / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  -
		  ( Ixy * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		    /
		    ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx),2))))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2))
		  +
		  tf::real_type(1));

	    tf::real_type d_Iyy =
		( Rxx / (Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx)
		  +
		  Ixx * (Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx)
		  /
		  ( pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx), 2)))
		/
		( pow((Ixx * Ryy - Ixy * Ryx - Iyx * Rxy + Iyy * Rxx), 2)
		  /
		  (pow((Rxx * Ryy - Rxy * Ryx - Ixx * Iyy + Ixy * Iyx),2))
		  +
		  tf::real_type(1));

	    tf::real_type Rxx_v = imp[i].Zxx().real().variance;
	    tf::real_type Rxy_v = imp[i].Zxy().real().variance;
	    tf::real_type Ryx_v = imp[i].Zyx().real().variance;
	    tf::real_type Ryy_v = imp[i].Zyy().real().variance;
	    tf::real_type Ixx_v = imp[i].Zxx().imag().variance;
	    tf::real_type Ixy_v = imp[i].Zxy().imag().variance;
	    tf::real_type Iyx_v = imp[i].Zyx().imag().variance;
	    tf::real_type Iyy_v = imp[i].Zyy().imag().variance;

	    tf::real_type uncertainty =
		sqrt( pow(d_Rxx, 2) * Rxx_v
		      +
		      pow(d_Rxy, 2) * Rxy_v
		      +
		      pow(d_Ryx, 2) * Ryx_v
		      +
		      pow(d_Ryy, 2) * Ryy_v
		      +
		      pow(d_Ixx, 2) * Ixx_v
		      +
		      pow(d_Ixy, 2) * Ixy_v
		      +
		      pow(d_Iyx, 2) * Iyx_v
		      +
		      pow(d_Iyy, 2) * Iyy_v )
		*
		tf::real_type(180) / pi;


 	    std::cout << period << "\t" << value << "\t"
		 << uncertainty << "\n";
    }
	else
	{
	    throw std::runtime_error("Parâmetro desconhecido: " + par );
	}
    }
} // lista_invariante
