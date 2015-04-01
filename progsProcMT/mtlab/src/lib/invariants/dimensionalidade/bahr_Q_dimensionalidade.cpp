#include <cmath>
#include <boost/format.hpp>
#include "dimensionalidade.hpp"
#include "transfer_functions/transfer_functions.hpp"
#include "invariants/erro_manual/base_xi-eta/invariants_parameters.hpp"
#include "invariants/erro_manual/base_xi-eta/wal_invariants.hpp"
#include "invariants/erro_manual/base_xi-eta/bahr_invariants.hpp"

dimensionalidade::dimensao
dimensionalidade::bahr_Q(
    transfer_function::impedance const& imp,
    dimensionalidade::thresholds t )
{
  namespace ip = erro_manual::base_xi_eta::invariant_parameter;
  namespace tf = transfer_function;
  namespace wal  = erro_manual::base_xi_eta::wal_invariant;
  namespace bahr = erro_manual::base_xi_eta::bahr_invariant;

  dimensionalidade::dimensao classe;
  tf::real_type kappa = bahr::kappa_value(imp);
  tf::real_type mu = bahr::mu_value(imp);
  tf::real_type eta = bahr::eta_value(imp);
  tf::real_type Sigma = bahr::Sigma_value(imp);
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

  // classificação segundo a tabela 2 de Marti, Queralt, Jones e Ledo
  // GJI (2005) v.163,pp 38-41
  // verifica se é 3D
  if( (fabs(eta) > t.eta) && (fabs(Q) > t.Q) )
    {
      classe="7  3D";
    }
  else if( fabs(eta) <= t.eta )
    {
      if( (fabs(kappa) <= t.kappa) && (fabs(mu) <= t.mu) )
	{
	  if(fabs(Sigma) <= t.Sigma)
	    {
	      classe="1  1D";
	    }
	  else
	    {
	      tf::real_type strikeRe = (atan(-xi3/xi2)/2.0)*(180.0/M_PI);
	      tf::real_type strikeIm = (atan(-eta3/eta2)/2.0)*(180.0/M_PI);
	      if(fabs(Q) <= t.Q)
		{
		  classe="2  2D[PhiTE==PhiTM]" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
		}
	      else
		{
		  classe="2  2D" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
		}
	    }
	}
      else if( (fabs(kappa) > t.kappa) && (fabs(Sigma) > t.Sigma) )
	{
	  if( fabs(mu) <= t.mu )
	    {
	      if( fabs(Q) > t.Q )
		{
		  tf::real_type strike = (atan( ( (xi1*eta2-xi2*eta1)-(xi3*eta4-xi4*eta3) ) / ( (xi1*eta3-xi3*eta1)+(xi2*eta4-xi4*eta2) ) )/2.0)*(180.0/M_PI);
		  classe="3  3D/2Dtwist" + str( boost::format(" %f") % strike );
		}
	      else
		{
		  classe="6  3D/1D_ou_2D[PhiTE==PhiTM]";
		}
	    }
	  else // eta==0 && (kappa&&Sigma&&mu>0)
	    {
	      tf::real_type strike = (atan( ( (xi1*eta2-xi2*eta1)-(xi3*eta4-xi4*eta3) ) / ( (xi1*eta3-xi3*eta1)+(xi2*eta4-xi4*eta2) ) )/2.0)*(180.0/M_PI);
	      classe="4  3D/2D" + str( boost::format(" %f") % strike );
	    }
	}
      else // condições para compatibilizar com WAL
	{
	  if( fabs(Q) <= t.Q )
	    {
              classe="6  3D/1D_ou_2D[PhiTE==PhiTM]";
	    }
	  else
	    {
	      classe="7  3D";
	    }
	}
    }
  else // ( eta > t.eta && Q <= t.Q )
    {
      if( (fabs(kappa) <= t.kappa) && (fabs(mu) <= t.mu) )
	{
	  if(fabs(Sigma) <= t.Sigma)
	    {
	      classe="1  1D";
	    }
	  else
	    {
	      tf::real_type strikeRe = (atan(-xi3/xi2)/2.0)*(180.0/M_PI);
	      tf::real_type strikeIm = (atan(-eta3/eta2)/2.0)*(180.0/M_PI);
	      classe="2  2D[PhiTE==PhiTM]" + str( boost::format(" %f") % strikeRe ) + str( boost::format(" %f") % strikeIm );
	    }
	}
      else if( (fabs(kappa) > t.kappa) && (fabs(mu) > t.mu) && (fabs(Sigma) > t.Sigma) )
	{
	  tf::real_type strike = (atan( ( (xi1*eta2-xi2*eta1)-(xi3*eta4-xi4*eta3) ) / ( (xi1*eta3-xi3*eta1)+(xi2*eta4-xi4*eta2) ) )/2.0)*(180.0/M_PI);
	  classe="4  3D/2D" + str( boost::format(" %f") % strike );
	}
      else
	{ // condição para compatibilizar com WAL
	  classe="6  3D/1D_ou_2D[PhiTE==PhiTM]";
	}
    }

  return classe;
} // Barh-Q
