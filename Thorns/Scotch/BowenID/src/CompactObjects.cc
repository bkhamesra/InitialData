#include "cctk_Parameters.h"
#include "CompactObjects.hh"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <OutputTable.hh>
#include "cctk.h"
#define M_TO_KM  1.47695518286//1.48
#define MSUN_DENS_TO_CGI 6.17327129839e17  //6.1694e17 // Convert rest-mass density to CGI for M=Msun

//Pure Virtual Destructors; required
Source::~Source() {}
CompactSource::~CompactSource() {}
BowenSource::~BowenSource() {}
SphericalHydroSource::~SphericalHydroSource() {}

CompactSource::CompactSource(int idx)
{
	DECLARE_CCTK_PARAMETERS;
	index = idx;
	position.x = object_rx[index];
	position.y = object_ry[index];
	position.z = object_rz[index];
}

double BowenSource::getNonFlatPsi(double xx, double yy, double zz)
{
	double rr = getRelativeDistance(xx, yy, zz);
	double psi = getNonFlatPsiOfR(smoothR(rr));
//	CCTK_VInfo(CCTK_THORNSTRING, "psi information at (%g, %g, %g) (r = %g): psi = %g", xx, yy, zz, rr, psi);
	return psi;
}

BowenSource::BowenSource(int idx) : CompactSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	momentum.x = object_Px[index];
	momentum.y = object_Py[index];
	momentum.z = object_Pz[index];
	spin.x = object_Jx[index];
	spin.y = object_Jy[index];
	spin.z = object_Jz[index];
}

VectorAlgebra::Vector BowenSource::getJhat(double xx, double yy, double zz)
{
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
	double rmag = smoothR(rr.getMagnitude());
	return (momentum*getMomentumDist(rmag) + VectorAlgebra::crossProduct(spin, rr)*getSpinDist(rmag));
}

VectorAlgebra::LinearOperator BowenSource::getAij(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
	double rmag = smoothR(rr.getMagnitude());

	VectorAlgebra::Vector direction = rr/rmag;

	VectorAlgebra::Vector lxs = VectorAlgebra::crossProduct(direction, spin);

	VectorAlgebra::LinearOperator psym(momentum.x*direction + direction.x*momentum, momentum.y*direction + direction.y*momentum, momentum.z*direction + direction.z*momentum);

	VectorAlgebra::LinearOperator project(direction.x*direction, direction.y*direction, direction.z*direction);

	VectorAlgebra::LinearOperator lxssym(lxs.x*direction + direction.x*lxs, lxs.y*direction+direction.y*lxs, lxs.z*direction+direction.z*lxs);

	
	double ldotp = VectorAlgebra::dotProduct(direction, momentum);
	

	double qq = getMomentumWeight(rmag);
	double cc = getPointSourceDeviation(rmag);
	double nn = getSpinWeight(rmag);

	VectorAlgebra::LinearOperator Aij_mom = 1.5*qq/SQR(rmag)*(psym+ ldotp*(project-VectorAlgebra::identity)) + 3.*cc/QUAD(rmag)*(psym + ldotp*(VectorAlgebra::identity - 5.*project));

	VectorAlgebra::LinearOperator Aij_spin = - 3.*nn/CUBE(rmag)*(lxssym);

	if(momentum.getMagnitude() < 1e-6)
	{
		Aij_mom = 0.*VectorAlgebra::identity;
	}
	if(spin.getMagnitude() < 1e-6)
	{
		Aij_spin = 0.*VectorAlgebra::identity;
	}

	VectorAlgebra::LinearOperator Aij = Aij_mom + Aij_spin;
	/*
	CCTK_VInfo(CCTK_THORNSTRING, "Aij info at (%g, %g, %g) (r = %g): q = %g, c = %g, n = %g, Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g", xx, yy, zz, rmag, qq, cc, nn, Aij.x.x, Aij.x.y, Aij.z.x, Aij.y.y, Aij.z.y, Aij.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "n = (%g, %g, %g), P = (%g, %g, %g), J = (%g, %g, %g)", direction.x, direction.y, direction.z, momentum.x, momentum.y, momentum.z, spin.x, spin.y, spin.z);
	CCTK_VInfo(CCTK_THORNSTRING, "n cross J = (%g, %g, %g), (n cross J) otimes n: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", lxs.x, lxs.y, lxs.z, lxssym.x.x, lxssym.y.x, lxssym.z.x, lxssym.x.y, lxssym.y.y, lxssym.z.y, lxssym.x.z, lxssym.y.z, lxssym.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "p otimes n: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", psym.x.x, psym.y.x, psym.z.x, psym.x.y, psym.y.y, psym.z.y, psym.x.z, psym.y.z, psym.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "projector: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", project.x.x, project.y.x, project.z.x, project.x.y, project.y.y, project.z.y, project.x.z, project.y.z, project.z.z);
	*/
	return Aij;
}

BlackHole::BlackHole(int idx) : BowenSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	mass = bh_bare_mass[index];
	radius = 0.;
}

double BowenSource::radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax)
{
	double fmin = (this->*dFdr)(rmin);
	double fmax = (this->*dFdr)(rmax);

	double rmid = (rmin+rmax)/2.0;
	double fmid = (this->*dFdr)(rmid);

	if(std::isnan(fmin) || std::isnan(fmid) || std::isnan(fmax))
	{
		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "radial integration failed from nans: f min, mid, max = (%g, %g, %g) at (%g, %g, %g)", fmin, fmid, fmax, rmin, rmid, rmax);
	}
	return simpsonIntegrate(dFdr, rmin, rmid, rmax, fmin, fmid, fmax, 0);
}

double BowenSource::simpsonIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmid, double rmax, double fmin, double fmid, double fmax, int depth)
{
	DECLARE_CCTK_PARAMETERS;
	double S1 = (fmin+4.*fmid+fmax)/6.*(rmax-rmin);

	double rq1 = (rmin+rmid)/2.0;
	double fq1 = (this->*dFdr)(rq1);
	double S21 = (fmin+4*fq1+fmid)/6.*(rmid-rmin);


	double rq3 = (rmax+rmid)/2.0;
	double fq3 = (this->*dFdr)(rq3);
	double S22 = (fmid+4*fq3+fmax)/6.*(rmax-rmid);

	double S2 = S21 + S22;

	double delta = (S2-S1)/15.;
	if(fabs(delta) < max_quadrature_error*pow(.5,depth) || depth > 3)
	{
		return S2 + delta;
	}
	else
	{
		return simpsonIntegrate(dFdr, rmin, rq1, rmid, fmin, fq1, fmid, depth+1) + simpsonIntegrate(dFdr, rmid, rq3, rmax, fmid, fq3, fmax, depth+1);
	}
}

/*
double BowenSource::radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax)
{
	DECLARE_CCTK_PARAMETERS;
	//adaptive simpson's rule (tentative; can implement other quadratures and use this function as a wrapper)	
	//TODO: a more precise integrator (e.g. Gauss-Kronrod?)
	double ffmin = (this->*dFdr)(rmin);
	double ffmax = (this->*dFdr)(rmax);

	double rmid = (rmin+rmax)/2.0;
	double ffmid = (this->*dFdr)(rmid);
	
	double rqtr1 = (rmin+rmid)/2.0;
	double ffqtr1 = (this->*dFdr)(rqtr1);

	double rqtr3 = (rmax+rmid)/2.0;
	double ffqtr3 = (this->*dFdr)(rqtr3);

	double S1 = simpson_weighted_integral(ffmin, ffmid, ffmax)*(rmax - rmin);
	double S2 = simpson_weighted_integral(ffmin, ffqtr1, ffmid)*(rmid-rmin) + simpson_weighted_integral(ffmid, ffqtr3, ffmax)*(rmax-rmid);

	bool too_small = ((rmax-rmin) < minimum_dr);
	if(too_small)
		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Simpison integration reached minimum size; results may be invalid!");

	double delta = (S2-S1)/15.;
	if(fabs(delta) < max_quadrature_error || too_small)
	{
		return S2 - delta;
	}
	else
	{
		return radialIntegrate(dFdr, rmin, rmid) + radialIntegrate(dFdr, rmid, rmax);
	}	
}
*/

double SphericalHydroSource::getRhohat(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx, yy, zz);
	double rho_H = std::pow(1.+getNonFlatPsiOfR(rr), rhohat_exponent)*((getRho(rr)+getPress(rr))*getWW2(xx, yy, zz) - getPress(rr));
//	double rhohat = std::pow(1.+getNonFlatPsiOfR(rr), rhohat_exponent)*getRho(rr);
	if(rho_H < 0.)
	{
		double ww2 = getWW2(xx, yy, zz);
		double press = getPress(rr);
		double rho = getRho(rr);
		double psi = 1.+getNonFlatPsiOfR(rr);

		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "rhohat negative at (%g, %g, %g): psi = %g, rho = %g, press = %g, ww2 = %g, exponent = %g, rhohat = %g", xx, yy, zz, psi, rho, press, ww2, rhohat_exponent, rho_H);
	}

	return rho_H;

}

double SphericalHydroSource::getConfDensity(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx, yy, zz);
	double rhohat = std::pow(1.+getNonFlatPsiOfR(rr), rhohat_exponent)*getRho(rr);
	if(rhohat < 0.)
	{
		double ww2 = getWW2(xx, yy, zz);
		double press = getPress(rr);
		double rho = getRho(rr);
		double psi = 1.+getNonFlatPsiOfR(rr);

		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "rhohat negative at (%g, %g, %g): psi = %g, rho = %g, press = %g, ww2 = %g, exponent = %g, rhohat = %g", xx, yy, zz, psi, rho, press, ww2, rhohat_exponent, rhohat);
	}

	return rhohat;

}

double SphericalHydroSource::getConfPressure(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx, yy, zz);
	double presshat = std::pow(1.+getNonFlatPsiOfR(rr), rhohat_exponent)*getPress(rr);
	if(presshat < 0.)
	{
		double ww2 = getWW2(xx, yy, zz);
		double press = getPress(rr);
		double rho = getRho(rr);
		double psi = 1.+getNonFlatPsiOfR(rr);

		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Conformal Pressure negative at (%g, %h, %g): psi = %g, rho = %g, press = %g, ww2 = %g, exponent = %g, presshat = %g", xx, yy, zz, psi, rho, press, ww2, rhohat_exponent, presshat);
	}

	return presshat;

}

double TOVStar::getConfK(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx, yy, zz);
	double Khat = std::pow(1.+getNonFlatPsiOfR(rr), 8.*(1-eos_gamma))*eos_k;

	return Khat;

}
TOVStar::TOVStar(int idx) : SphericalHydroSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	eos_k = tov_K[index];
	eos_gamma = tov_Gamma[index];
	double nn = 1./(eos_gamma-1.);
	double kfactor = std::pow(eos_k, -nn);
	dr = tov_dr[index];

	std::stringstream ofile;
	ofile << out_dir << "/tov_internaldata_" << index << ".asc";

	OutputTable otable(ofile, 15, std::ios_base::out);
	if(testsuite)
	{
	//BK Modification - rho0, Lorentz factor added
		otable << "#1:riso" << "2:rho0" << "3:press" << "4:psi" << "5:theta" << "6:Lorentz Fac" << "7:qq" << "8:cc" << "9::nn";
		otable << std::endl;
	}

	//now using dimensionless densities
	double rho_c = tov_rho_central[index]*kfactor;
	double rho_atmo = tov_rho_atmosphere[index]*kfactor;
	double press_c = getPressFromEOS(rho_c);
	double press_atmo = getPressFromEOS(rho_atmo);

	tovdata_t temp;
	temp.press = press_c;
    temp.psi0 = 1.;  // will be adjusted later on
    temp.psi1 = 0.;                           
	temp.mass0 = 0.;
    temp.theta0 = 1.; // will be adjusted later on
    temp.theta1 = 0.;                           
    temp.riso = 0.;                          

	rk4data.push_back(temp);
    int i = 0;
    do
    {
	if(i == tov_numpoints)
	    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Integration of hydrodynamic equations failed for TOV initial data: domain too small, still inside star after %d points", tov_numpoints);

	rk4data.push_back(tov_rkstep(i*dr, dr, rk4data[i]));
    } while(rk4data[++i].press > press_atmo);

    inside_points = i;
    int last_point = i-1;
	rk4data[last_point].press = press_atmo;

    double rad_unscaled =rk4data[last_point].riso;
    double psi_rescale_factor = ((rk4data[last_point]).psi0+rk4data[last_point].psi1*rad_unscaled);
    double theta_rescale_factor = ((rk4data[last_point]).theta0+rk4data[last_point].theta1*rad_unscaled);

	mbar = rk4data[last_point].qq/SQR(psi_rescale_factor);
	ibar = rk4data[last_point].nn*SQR(psi_rescale_factor);

	otable << std::scientific;
    for(int l = 0; l < inside_points; l++)
    {
		rk4data[l].theta0 /= theta_rescale_factor;
		rk4data[l].psi0 /= psi_rescale_factor;
		rk4data[l].riso *= SQR(psi_rescale_factor);
		rk4data[l].qq /= (mbar*SQR(psi_rescale_factor));
		rk4data[l].cc *= SQR(psi_rescale_factor)/mbar;
		rk4data[l].nn *= SQR(psi_rescale_factor)/ibar;

    }
    //BK Modification - For loop with density and lorentz factor output added, original output moved here
    for(int l = 0; l < inside_points; l++)
    {
	double ww2_linear = getWW2(position.x+rk4data[l].riso, position.y, position.z);
	double ww_linear = sqrt(ww2_linear);
	double  rho0 = getRho0FromEOS(rk4data[l].press);
	if(testsuite)
	{
		otable << rk4data[l].riso << rho0 << rk4data[l].press << rk4data[l].psi0 << rk4data[l].theta0 << ww_linear << rk4data[l].qq << rk4data[l].cc << rk4data[l].nn;
		otable << std::endl;
	}
    }

    dr *= SQR(psi_rescale_factor);
    
    radius = rk4data[last_point].riso;
    mass = 2.*radius*(rk4data[last_point].psi0-1.);
    mass0 = rk4data[last_point].mass0;
/*
	mbar = 1.;
	ibar = 1.;

	mbar = getMomentumWeight(radius);
	ibar = getSpinWeight(radius);
	*/

    if(verbose >= 1)
    {
		double ww2_linear = getWW2(position.x, position.y, position.z);
		CCTK_INFO("Information on solved TOV star (where Ms are Msys in solar masses):"); 
		CCTK_VInfo(CCTK_THORNSTRING, "TOV #%d info: Mass = %g*Ms M_sun, Rest mass = %g, Mbar = %g, Ibar = %g, Radius = %g*Ms (%g*Ms km,Schwarzschild), rho_central = %g*Ms^(-2) g/cm^3, v*psi^2 = (%g, %g, %g), omega*psi^2 = (%g, %g, %g), ww2_linear = %g", 
				index, mass, mass0,
				mbar, ibar,
				radius, 
				M_TO_KM*radius*SQR(rk4data[last_point].psi0),
				rho_c*MSUN_DENS_TO_CGI, momentum.x/mbar, momentum.y/mbar, momentum.z/mbar, 
				spin.x/ibar, spin.y/ibar, spin.z/ibar, ww2_linear);
    }
}

tovdata_t TOVStar::getRHS(double riso, tovdata_t x)
{
	tovdata_t v;

	double rho = getRhoFromEOS(x.press);
	double rho0 = getRho0FromEOS(x.press);


	v.psi1 = -2.0*M_PI*pow(x.psi0,5.0)*rho;
	v.theta1 = 2.0*M_PI*x.theta0*pow(x.psi0,4.0)*(rho+6.0*x.press);
	v.psi0 = x.psi1;
	v.theta0 = x.theta1;
	v.riso = 1.0;
	v.qq = 4.*M_PI*(rho+x.press)*pow(x.psi0,8.)*SQR(riso);
	v.cc = 2.*M_PI/3*(rho+x.press)*pow(x.psi0,8.)*QUAD(riso);
	v.nn = 8.*M_PI/3*(rho+x.press)*pow(x.psi0,8.)*QUAD(riso);
	v.mass0 = 4.0*M_PI*pow(x.psi0,6.)*SQR(riso)*rho0;
	
	if(riso == 0.0)
	{
		v.press = 0.;
		v.psi1 += 4.0*M_PI*rho/3.0;
		v.theta1 += -4.0*M_PI*(rho+6.*x.press)/3.0;
		//there probably ought to be some conditions here for v.qq, v.cc; see conditions for v.psi1, v.theta1
	}
	else
	{
		v.press = -(rho + x.press)*(x.theta1/x.theta0 - x.psi1/x.psi0);
		v.theta1 += -2.0*x.theta1/riso;
		v.psi1 += -2.0*x.psi1/riso;
	}
	return v;
}

tovdata_t TOVStar::tov_rkstep(double x, double h, tovdata_t y_in)
{
    tovdata_t k1, k2, k3, k4, y_out;
    tovdata_t tmp;

    // k1
    k1 = getRHS(x, y_in)*h;

    // k2
    tmp = y_in + k1/2.;
    k2 = getRHS(x + h/2, tmp)*h;

    // k3
    tmp = y_in + k2/2.;
    k3 = getRHS(x + h/2, tmp)*h;

    // k4
    tmp = y_in + k3;
    k4 = getRHS(x + h, tmp)*h;

    // combine them
    y_out = y_in + 1/6. * (k1 + 2*k2 + 2*k3 + k4);
    return y_out;
}

tovdata_t TOVStar::interpolateData(double distance)
{
    tovdata_t interp;
    int last_point = inside_points - 1;
    int i_low1, i_high1, i_low2, i_high2;
    i_low1 = (int) (distance/dr);
    i_high1 = i_low1 + 1;
	i_low2 = i_low1-1;
	i_high2 = i_high1+1;

	if(i_low2 < 0)
	{
		i_low2 = 0;
		i_low1 = 1;
		i_high1 = 2;
		i_high2 = 3;
	}
	if(i_high2 == inside_points)
	{
		i_high2 = inside_points-1;
		i_high1 = i_high2-1;
		i_low1 = i_high1-1;
		i_low2 = i_low1-1;	
	}
    if(i_high2 <= last_point)
    {
		// box desired value such that riso[i_min] <= riso < riso[i_min+1]
		//
		/*
		double driso, delta;

		driso = distance - rk4data[i_min].riso; // distance from desired value
		delta = rk4data[i_max].riso - rk4data[i_min].riso;  // step in riso 
		double minweight, midweight, maxweight;
		maxweight = (driso/delta);
		minweight = 1. - maxweight;
		interp = minweight*rk4data[i_min]+ maxweight*rk4data[i_max];
		*/
		double c_low2, c_low1, c_high1, c_high2;
		double x_low2, x_low1, x_high1, x_high2;

		x_low2 = rk4data[i_low2].riso;
		x_low1 = rk4data[i_low1].riso;
		x_high1 = rk4data[i_high1].riso;
		x_high2 = rk4data[i_high2].riso;

		c_low2 = (distance-x_low1)*(distance-x_high1)*(distance-x_high2)/((x_low2-x_low1)*(x_low2-x_high1)*(x_low2-x_high2));
		c_low1 = (distance-x_low2)*(distance-x_high1)*(distance-x_high2)/((x_low1-x_low2)*(x_low1-x_high1)*(x_low1-x_high2));
		c_high1 = (distance-x_low2)*(distance-x_low1)*(distance-x_high2)/((x_high1-x_low2)*(x_high1-x_low1)*(x_high1-x_high2));
		c_high2 = (distance-x_low2)*(distance-x_high1)*(distance-x_low1)/((x_high2-x_low2)*(x_high2-x_high1)*(x_high2-x_low1));
		
		interp = c_low2*rk4data[i_low2]+c_low1*rk4data[i_low1]+c_high1*rk4data[i_high1]+c_high2*rk4data[i_high2];
	}
    else
	{
		interp = rk4data[last_point];
		interp.press = 0.;
		interp.psi0 = 1. + (interp.psi0 - 1.)*(radius/distance);
		interp.theta0 = 1. - (1. - interp.theta0)*(radius/distance);
		interp.riso = distance;
	}

	if(std::isnan(interp.press) || std::isnan(interp.psi0) || std::isnan(interp.theta0) || std::isnan(interp.riso))
	{
		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "interpolation of tov data returned nan: press = %g, psi = %g, theta = %g, riso = %g", interp.press, interp.psi0, interp.theta0, interp.riso);
	}
    return interp;
}

//not actually sure this is right

double SphericalHydroSource::getWW2(double xx,double yy,double zz)
{
	double dist = getRelativeDistance(xx, yy, zz);
	double psi = 1.+getNonFlatPsiOfR(dist);
	double jhatmag = getJhat(xx, yy, zz).getMagnitude();
	if(jhatmag > 0.)
	{
		double ww2 = .5+sqrt(.25+SQR(jhatmag/( (getRho(dist) + getPress(dist))*pow(1.+getNonFlatPsiOfR(dist), 8.) ) )); //that is, jhatmag/((rho+p)*Psi^8)
//		CCTK_VInfo(CCTK_THORNSTRING, "ww2 info at (%g, %g, %g): |jhat| = %g, sigma_p = %g, sigma_j = %g, Mbar = %g, ww2 = %g", xx, yy, zz, jhatmag, getMomentumDist(dist), getSpinDist(dist), mbar, ww2);
		return ww2;
	}
	else
	{
		return 1.;
	}
}

double TOVStar::getMomentumWeight(double rr)
{
	/*
	double realrad = rr > radius ? radius : rr;
	return radialIntegrate(&BowenSource::getMomentumWeightDeriv, 0., realrad);
	*/
	return interpolateData(rr).qq;
}

double TOVStar::getPointSourceDeviation(double rr)
{
	/*
	double realrad = rr > radius ? radius : rr;
	return radialIntegrate(&BowenSource::getPSDDeriv, 0., realrad);
	*/
	return interpolateData(rr).cc;
}

double TOVStar::getSpinWeight(double rr)
{
	/*
	double realrad = rr > radius ? radius : rr;
	return radialIntegrate(&BowenSource::getSpinWeightDeriv, 0., realrad);
	*/
	return interpolateData(rr).nn;
}

GaussianCloud::GaussianCloud(int idx) : SphericalHydroSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	amplitude = gaussian_amp[index];
	radius = gaussian_radius[index];
	mbar = 1.;
	ibar = 1.;

	mbar = getMomentumWeight(radius);
	ibar = getSpinWeight(radius);
}

double GaussianCloud::getRho(double rr)
{
	if(rr < radius)
	{
		double psi_nf = getNonFlatPsiOfR(rr);
		double psi = 1.+psi_nf;
		return 9.*psi_nf*(SQR(radius)-SQR(rr))/QUAD(radius)/(-2.*M_PI*std::pow(psi, 5.));
	}
	else
	{
		return 0.;
	}
}

//NOTICE: this is consistent with the RHS routine, but not used there
double SphericalHydroSource::getMomentumDist(double rr)
{
	DECLARE_CCTK_PARAMETERS;
	return (getRho(rr)+getPress(rr))*pow(1.+getNonFlatPsiOfR(rr), sigma_exponent)/mbar;
} //"sigma_p"


tovdata_t::tovdata_t() {press=psi0=theta0=psi1=theta1=mbar=mass=qq=cc=riso=nn=mass0=0.0;}

tovdata_t::tovdata_t(const tovdata_t& a) {press = a.press;psi0 = a.psi0;theta0 = a.theta0;psi1 = a.psi1;theta1 = a.theta1; riso = a.riso; mbar = a.mbar;qq = a.qq; cc = a.cc; mass = a.mass; nn=a.nn; mass0=a.mass0; }

tovdata_t& tovdata_t::operator=(const tovdata_t& a) { press = a.press; psi0 = a.psi0; theta0 = a.theta0; psi1 = a.psi1; theta1 = a.theta1; riso = a.riso; mbar = a.mbar; qq = a.qq; cc = a.cc; mass = a.mass; nn = a.nn; mass0=a.mass0; return *this; }

tovdata_t& tovdata_t::operator*=(double a) {press *= a;psi0 *= a; psi1 *= a; theta0 *= a;theta1*=a; riso *= a; mbar *= a; qq *= a; cc *= a; mass *= a; nn *= a; mass0 *= a; return *this;}

tovdata_t tovdata_t::operator*(double a) const {return tovdata_t(*this).operator*=(a);}

tovdata_t& tovdata_t::operator/=(double a) {press /= a;psi0 /= a;theta0 /= a;psi1 /= a; theta1 /= a; riso /= a; mbar /= a; qq /= a; cc /= a; mass /= a; nn /= a;  mass0 /= a; return *this;}

tovdata_t tovdata_t::operator/(double a) const {return tovdata_t(*this).operator/=(a);}

tovdata_t& tovdata_t::operator+=(const tovdata_t &a) {press += a.press;psi0 += a.psi0;theta0 += a.theta0;psi1 += a.psi1; theta1 += a.theta1; riso += a.riso; mbar += a.mbar; qq += a.qq; cc += a.cc; mass += a.mass; nn += a.nn; mass0 += a.mass0; return *this;}

tovdata_t tovdata_t::operator+(const tovdata_t &a) const {return tovdata_t(*this).operator+=(a);}

tovdata_t operator*(double a, const tovdata_t& b) {return tovdata_t(b).operator*=(a); }
