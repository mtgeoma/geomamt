#include <cmath>
#include <boost/format.hpp>
#include "dimensionalidade.hpp"
#include "transfer_functions/transfer_functions.hpp"
#include "invariants/erro_manual/base_xi-eta/invariants_parameters.hpp"
#include "invariants/erro_manual/base_xi-eta/wal_invariants.hpp"

dimensionalidade::dimensao
dimensionalidade::wal(
    transfer_function::impedance const& imp,
    dimensionalidade::thresholds t )
{
  namespace ip = erro_manual::base_xi_eta::invariant_parameter;
  namespace tf = transfer_function;
  namespace wal  = erro_manual::base_xi_eta::wal_invariant;

  dimensionalidade::dimensao classe;
  tf::real_type I3 = wal::I3_value(imp);
  tf::real_type I4 = wal::I4_value(imp);
  tf::real_type I5 = wal::I5_value(imp);
  tf::real_type I6 = wal::I6_value(imp);
  tf::real_type I7 = wal::I7_value(imp);
  tf::real_type Q  = wal::Q_value(imp);
  tf::real_type cos_alpha = wal::cos_alpha_value(imp);
  tf::real_type cos_beta = wal::cos_beta_value(imp);
  tf::real_type xi1 = ip::xi1_value(imp);
  tf::real_type xi2 = ip::xi2_value(imp);
  tf::real_type xi3 = ip::xi3_value(imp);
  tf::real_type xi4 = ip::xi4_value(imp);
  tf::real_type eta1 = ip::eta1_value(imp);
  tf::real_type eta2 = ip::eta2_value(imp);
  tf::real_type eta3 = ip::eta3_value(imp);
  tf::real_type eta4 = ip::eta4_value(imp);

  // verifica se é 3D
  if( (fabs(I7) >= t.I7) && (fabs(Q) >= t.Q) )
    {
      classe="7  3D";
    }
  else if( fabs(Q) <= t.Q )
    {
      if( fabs(I6) <= t.I6 )
	{
	  if( fabs(I5) <= t.I5 )
	    {
	      if( (fabs(I4) <= t.I4) && (fabs(I3) <= t.I3) )
		{
		  classe="1  1D";
		}
	      else
		{
		  if( ( fabs(cos_alpha) <= t.cos_alpha ) &&
		      ( fabs(cos_beta) <= t.cos_beta ) )
		    {
		      classe="0  classe_WAL_indefinida[Q==0&&zeta4==0]";
		    }
		  else
		    {
		      tf::real_type strikeRe = (atan(-xi3/xi2)/2.0)*(180.0/M_PI);
		      tf::real_type strikeIm = (atan(-eta3/eta2)/2.0)*(180.0/M_PI);
		      classe="2  2D[PhiTE==PhiTM]" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
		    }
		}
	    }
	  else
	    {
	      classe="6  3D/1D_ou_2D[PhiTE==PhiTM]";
	    }
	}
      else // classe indefinida; parece-me que WAL assume Q==0 => d41==0 (I6==0)
	{  // mas não consegui estabelecer essa relação 
	  classe="0  classe_WAL_indefinida[Q==0&&I6!=0]";
	}
    }
  else // I7==0
    {
      if( fabs(I6) <= t.I6 )
	{
	  if( fabs(I5) <= t.I5 )
	    {
	      if( (fabs(I4) <= t.I4) && (fabs(I3) <= t.I3) )
		{
		  classe="1  1D[Q!=0(?)]";
		}
	      else
		{
		  if( ( fabs(cos_alpha) <= t.cos_alpha ) &&
		      ( fabs(cos_beta) <= t.cos_beta ) )
		    {
		      tf::real_type strikeRe = (atan(xi2/xi3)/2.0)*(180.0/M_PI);
		      tf::real_type strikeIm = (atan(eta2/eta3)/2.0)*(180.0/M_PI);
		      classe="5  3D/1D2Ddiag" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
		    }
		  else
		    {
		      tf::real_type strikeRe = (atan(-xi3/xi2)/2.0)*(180.0/M_PI);
		      tf::real_type strikeIm = (atan(-eta3/eta2)/2.0)*(180.0/M_PI);
		      classe="2  2D" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
		    }
		}
	    }
	  else
	    {
	      tf::real_type strike = (atan( ( (xi1*eta2-xi2*eta1)-(xi3*eta4-xi4*eta3) ) / ( (xi1*eta3-xi3*eta1)+(xi2*eta4-xi4*eta2) ) )/2.0)*(180.0/M_PI);
	      classe="3  3D/2Dtwist" + str( boost::format(" %f") % strike );
	    }
	}
      else
	{
	  tf::real_type strike = (atan( ( (xi1*eta2-xi2*eta1)-(xi3*eta4-xi4*eta3) ) / ( (xi1*eta3-xi3*eta1)+(xi2*eta4-xi4*eta2) ) )/2.0)*(180.0/M_PI);
	  classe="4  3D/2D" + str( boost::format(" %f") % strike );
	}
    }

  return classe;
} // wal
