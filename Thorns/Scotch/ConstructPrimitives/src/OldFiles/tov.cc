// vim: sw=4 smarttab smartindent softtabstop=4
// tov.c
//
// solves the radial pressure equations for the TOV star as described by
// L\"offler, Hawke and Baumgarte in the documentation to the TOVSolverC thorn
// for Whisky :
// (6) dP/dr = -(mu+P) (m + 4 pi r^2 P) / (r (r - 2m))
// (7) dm/dr = 4 pi r^2 mu
// (8) dphi/dr = (m + 4 pi r^3 P) / (r (r - 2m))
// (25) riso = C r exp(int_0^r (1-sqrt(1-2m/r))/(r*sqrt(1-2m/r)))
//      C = (sqrt(R^2-2MR)+R-M)/(2R) exp(-int_0^R (1-sqrt(1-2m/r))/(r*sqrt(1-2m/r)))
//      M=m(R), R surface of star
// starting from r = 0 with rho=cloud_rho[0] up to the surface of the star where rho
// drops below cloud_rho_atmosphere
// uses EOS:
// mu = rho (1+epsilon)
// P  = (cloud_Gamma - 1) rho epsilon
// P  = cloud_K rho^cloud_Gamma
//
// Author: Roland Haas
// File created on  Sept 17, 2008
// Last Change: Thu Sep 18 09:20 2008
// 
// Copyright (C) 2008 Roland Haas
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

// NOTE: Whisky is GPL so if we use it, and distribute the software so must we.
// If we don't distribute, we don't (see the gsl docs at:
// http://www.gnu.org/software/gsl/)

#define DEBUG

#include "tov.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "matterlibs.h"

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <vector>
#include <cmath>

#define NUM_STARS 10 // max. number of stars we can generate

#define M_TO_KM 1.48
#define MSUN_DENS_TO_CGI 6.1694e17 // Convert rest-mass density to CGI for M=Msun
#define PI 3.14159265358979

static inline double SQR(const double x) {return x*x;}
static inline double QUAD(const double x) {const double x2=x*x;return x2*x2;}
static inline double CUBE(const double x) {return x*x*x;}
#define DIM(v) (sizeof(v)/sizeof((v)[0]))
static inline double MIN(const double x, const double y) {if(x < y) return x; else return y;}


extern "C"
CCTK_REAL matterlibs_getTOVMass(CCTK_INT index)
{
    DECLARE_CCTK_PARAMETERS;
    TOVController &controller = TOVController::create(num_clouds);

    return controller.getMass(index);
}

extern "C"
CCTK_REAL matterlibs_getTOVRadius(CCTK_INT index)
{
    DECLARE_CCTK_PARAMETERS;
    TOVController &controller = TOVController::create(num_clouds);

    return controller.getRadius(index);
}


extern "C"
void matterlibs_updateTOVs()
{
    DECLARE_CCTK_PARAMETERS;
    TOVController &controller = TOVController::create(num_clouds);

    controller.updateTOVs();
}

extern "C"
void matterlibs_preSolveTOVs(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_PARAMETERS;

    TOVController &controller = TOVController::create(num_clouds);
}

extern "C"
void matterlibs_calcTOVQCNRhoPsi(int star, double distance, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat)
{
    DECLARE_CCTK_PARAMETERS;
    TOVController &controller = TOVController::create(num_clouds); 
    controller.getQCNRhoPsi(star, distance, qq, cc, nn, rho, psi_nonflat);
}

double TOVController::getLapse(double x, double y, double z)
{
    DECLARE_CCTK_PARAMETERS;	
    double lapse = 1.;
    for(int i = 0; i < num_clouds; i++)
    {
	lapse *= stars.at(i).getLapse(x, y, z);
    }
    return lapse;
}

double TOVController::getMass(int index)
{
    return stars.at(index).getMass();
}

double TOVController::getRadius(int index)
{
    return stars.at(index).getRadius();
}

void TOVController::getQCNRhoPsi(int index, double distance, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat)
{
    stars.at(index).getQCNRhoPsi(distance, qq, cc, nn, rho, psi_nonflat);
}

extern "C"
CCTK_REAL calcTOVLapse(double x, double y, double z)
{
    DECLARE_CCTK_PARAMETERS;
    TOVController &controller = TOVController::create(num_clouds); 
    return controller.getLapse(x, y, z);
}

void TOVController::getConformalRho(double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz)
{
	DECLARE_CCTK_PARAMETERS;
	*dpsi = 1.;
	*drho0 = 0.;
	*dEbar = 0.;
	*dJbarx = 0.;
	*dJbary = 0.;
	*dJbarz = 0.;
	*dlorentz = 0.;

	int inside_star = -1;
	for(int i = 0; i < num_clouds; i++)
	{

	    stars.at(i).addConformalRho(x, y, z, drho0, dpsi, dEbar, dJbarx, dJbary, dJbarz, dlorentz);
	}
	
}

void TOVController::getRho(double x, double y, double z, double *dpsi, double *drho0, double *deps, double *dpress, double *dv_x, double *dv_y, double *dv_z)
{
	DECLARE_CCTK_PARAMETERS;
	*dpsi = 1.;
	*drho0 = 0.;
	*deps = 0.;
	*dpress = 0.;
	*dv_x = 0.;
	*dv_y = 0.;
	*dv_z = 0.;

	int inside_star = -1;
	for(int i = 0; i < num_clouds; i++)
	{

	    stars.at(i).addRho(x, y, z, dpsi, drho0, deps, dpress, dv_x, dv_y, dv_z);
	}
}
extern "C"
void calcTOVRho( double x, double y, double z, 
        double *dpsi, double *drho0, double *deps, double *dpress, 
        double *dv_x, double *dv_y, double *dv_z ) 
{
	DECLARE_CCTK_PARAMETERS;
	TOVController &controller = TOVController::create(num_clouds); 
	controller.getRho(x, y, z, dpsi, drho0, deps, dpress, dv_x, dv_y, dv_z);

	//CCTK_VInfo(CCTK_THORNSTRING, "calcTOVRho data at (%g, %g, %g) - rho0 = %g, press = %g, eps = %g, psi = %g", x, y, z, *drho0, *dpress, *deps, *dpsi);
}
extern "C"
void calcConformalTOVRho( double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz ) 
{
	DECLARE_CCTK_PARAMETERS;
	TOVController &controller = TOVController::create(num_clouds); 
	controller.getConformalRho(x, y, z, dpsi, drho0, dEbar, dJbarx, dJbary, dJbarz, dlorentz);

	//CCTK_VInfo(CCTK_THORNSTRING, "calcTOVRho data at (%g, %g, %g) - rho0 = %g, press = %g, eps = %g, psi = %g", x, y, z, *drho0, *dpress, *deps, *dpsi);
}

TOVController& TOVController::create(int count)
{
    static TOVController instance(count);
    return instance;
}

TOVController::~TOVController()
{
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,  "TOVController destroyed!");
}

TOVController::TOVController(int count)
{
    numstars = count;
    updateTOVs();
}

void TOVController::updateTOVs()
{
    DECLARE_CCTK_PARAMETERS;
    if(numstars > 1)
	CCTK_WARN(1, "Specifying more than one TOV star may lead to unphysical behaviour.");

    FILE *testsuitefh;
    if(testsuite)
    {
	char fn[1024];

	snprintf(fn, sizeof(fn)/sizeof(fn[0]), "%s/tov_internaldata.asc", out_dir);
	if((testsuitefh = fopen(fn, "w")) == NULL)
	    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Cannot open file %s to dump internal data", fn);

	fprintf(testsuitefh, "# 1:riso 2:press 3:psi 4:theta 5:mass 6:mbar 7:Q 8:C\n");
    }

    stars.clear();
    for(int i = 0; i < numstars; i++)
    {
	TOVStar st(i, testsuitefh);
	stars.push_back(st);
    }

    if(testsuite)
    {
	fclose(testsuitefh);
	testsuitefh = NULL;
    }
}

TOVStar::TOVStar(int index, FILE *testsuitefh)
{
    DECLARE_CCTK_PARAMETERS;
    star_index = index;
    if(testsuite)
        fprintf(testsuitefh, "# star %d\n", star_index);
    numpoints = cloud_numpoints[star_index];
    tov_k = cloud_K[star_index];
    tov_gamma = cloud_Gamma[star_index];
    dr = cloud_dr[star_index];

    position.x = cloud_x0[star_index];
    position.y = cloud_y0[star_index];
    position.z = cloud_z0[star_index];

    if (tov_centralrho[star_index] > 0.0)
	rho_c = tov_centralrho[star_index];
    else
	rho_c = cloud_rho[star_index];

    if (initial_rho_abs_min > 0.0)
	rho_atmosphere = initial_rho_abs_min;
    else if (initial_rho_rel_min > 0.0)
	rho_atmosphere = rho_c * initial_rho_rel_min;
    else if (rho_abs_min > 0.0)
	rho_atmosphere = rho_abs_min;
    else
	rho_atmosphere = rho_c * rho_rel_min;

    if (initial_atmosphere_factor > 0.0)
	rho_atmosphere *= initial_atmosphere_factor;

    double press_c = tov_k*pow(rho_c, tov_gamma);
    double press_atmosphere = tov_k*pow(rho_atmosphere, tov_gamma);
   //Solving routine
    
   
    rk4data_t temp;
    temp.press = press_c;
    temp.psi0 = 1.;  // will be adjusted later on
    temp.psi1 = 0.;                           
    temp.theta0 = 1.; // will be adjusted later on
    temp.theta1 = 0.;                           
    temp.riso = 0.;                          
    temp.qq = 0.;
    temp.cc = 0.;
    temp.mass = 0.;

    tovdata.push_back(temp);
    int i = 0;
    do
    {
	if(i == numpoints)
	    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Integration of hydrodynamic equations failed for TOV initial data: domain too small, still inside star after %d points", numpoints);

        tovdata.push_back(rkstep(i*dr, dr, tovdata[i], &TOVStar::TOV_RHS));
//	CCTK_VInfo(CCTK_THORNSTRING, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", tovdata[i+1].riso, tovdata[i+1].press, tovdata[i+1].psi0, tovdata[i+1].theta0, tovdata[i+1].mass, tovdata[i+1].mbar, tovdata[i+1].qq, tovdata[i+1].cc);
    } while(tovdata[++i].press > press_atmosphere);

    //fill in fixed data;
    inside_points = i;
    int last_point = i-1;
    tovdata[last_point].press = press_atmosphere;

    double qmax = tovdata[last_point].qq;
    
    double rad_unscaled = tovdata[last_point].riso;

    //single-pass rescaling
    double psi_rescale_factor = ((tovdata[last_point]).psi0+tovdata[last_point].psi1*rad_unscaled);
  //  CCTK_VInfo(CCTK_THORNSTRING, "Rescaling information: psi(R) = %g, dpsi(R) = %g, R(unscaled) = %g, Rescale factor = %g", tovdata[last_point].psi0, tovdata[last_point].psi1, rad_unscaled, psi_rescale_factor);
    double theta_rescale_factor = ((tovdata[last_point]).theta0+tovdata[last_point].theta1*rad_unscaled);

    mrest = tovdata[last_point].mrest;
    mass = tovdata[last_point].mass;
    mbar = tovdata[last_point].mbar/SQR(psi_rescale_factor);
    for(int l = 0; l < inside_points; l++)
    {
		tovdata[l].mbar /= SQR(psi_rescale_factor);
		tovdata[l].theta0 /= theta_rescale_factor;
		tovdata[l].psi0 /= psi_rescale_factor;
		tovdata[l].riso *= SQR(psi_rescale_factor);
		tovdata[l].cc *= QUAD(tovdata[l].riso)/CUBE(SQR(psi_rescale_factor))/mbar;
		tovdata[l].qq *= SQR(tovdata[l].riso)/CUBE(SQR(psi_rescale_factor))/mbar;
		if(testsuite)
		{
			fprintf(testsuitefh, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", tovdata[l].riso, tovdata[l].press, tovdata[l].psi0, tovdata[l].theta0, tovdata[l].mass, tovdata[l].mbar, tovdata[l].qq, tovdata[l].cc);
		}

    }
    radius = tovdata[last_point].riso;
	ibar = tovdata[last_point].cc*4.*mbar;

    dr *= SQR(psi_rescale_factor);
    //setting velocities
	if(provide_momentum)
	{
		momentum.x = cloud_px[star_index];
		momentum.y = cloud_py[star_index];
		momentum.z = cloud_pz[star_index];
		angular_momentum.x = cloud_jx[star_index];
		angular_momentum.y = cloud_jy[star_index];
		angular_momentum.z = cloud_jz[star_index];
	}
	else
	{
		double wlinear = 1/sqrt(1 - SQR(cloud_vx[0]) - SQR(cloud_vx[1]) - SQR(cloud_vx[2]));
		momentum.x = mbar*cloud_vx[0]*wlinear*wlinear;
		momentum.y = mbar*cloud_vy[1]*wlinear*wlinear;
		momentum.z = mbar*cloud_vz[2]*wlinear*wlinear;

		//implicitly assuming rotation is "small"
		angular_momentum.x = ibar*cloud_omega_x[0];
		angular_momentum.y = ibar*cloud_omega_y[0];
		angular_momentum.z = ibar*cloud_omega_z[0];
	}
	double linear_speed = momentum.getMagnitude()/mbar;
	w_linear = 1./(1.-SQR(linear_speed));

    if(verbose >= 1)
    {
	CCTK_INFO("Information on solved TOV star (where Ms are Msys in solar masses):"); 
	CCTK_VInfo(CCTK_THORNSTRING, "TOV #%d info: Mass = %g*Ms M_sun, Rest mass = %g, Mbar = %g, Ibar = %g, Radius = %g*Ms (%g*Ms km,Schwarzschild), rho_central = %g*Ms^(-2) g/cm^3, v*psi^2 = (%g, %g, %g), omega*psi^2 = (%g, %g, %g)", 
		star_index, mass,
		mrest, mbar, ibar,
		radius, 
		M_TO_KM*radius*SQR(tovdata[last_point].psi0),
		rho_c*MSUN_DENS_TO_CGI, momentum.x/mbar, momentum.y/mbar, momentum.z/mbar, 
		angular_momentum.x/ibar, angular_momentum.y/ibar, angular_momentum.z/ibar);
    }
}

rk4data_t TOVStar::rkstep(double x, double h, rk4data_t y_in,  rk4data_t (TOVStar::*f)(double x, rk4data_t y) const)
{
    rk4data_t k1, k2, k3, k4, y_out;
    rk4data_t tmp;

    // k1
    k1 = (this->*f)(x, y_in)*h;

    // k2
    tmp = y_in + k1/2.;
    k2 = (this->*f)(x + h/2, tmp)*h;

    // k3
    tmp = y_in + k2/2.;
    k3 = (this->*f)(x + h/2, tmp)*h;

    // k4
    tmp = y_in + k3;
    k4 = (this->*f)(x + h, tmp)*h;

    // combine them
    y_out = y_in + 1/6. * (k1 + 2*k2 + 2*k3 + k4);
    return y_out;
}

rk4data_t TOVStar::TOV_RHS(double riso, rk4data_t x) const
{
	rk4data_t v;

	double rho0, eps, rho;
	getRhoEpsOfPress(x.press, &rho0, &eps, &rho);    
	//    CCTK_VInfo(CCTK_THORNSTRING, "X in at riso = %g: %.15e %.15e %.15e %.15e %.15e %.15e\n", riso, x.psi0, x.psi1, x.theta0, x.theta1, x.press, x.mass);
	double rhobar = (rho + x.press)*pow(x.psi0, 8.);

	v.psi1 = -2.0*PI*pow(x.psi0,5.0)*rho;
	v.theta1 = 2.0*PI*x.theta0*pow(x.psi0,4.0)*(rho+6.0*x.press);
	v.psi0 = x.psi1;
	v.theta0 = x.theta1;
	v.riso = 1.0;
	v.qq = 4.*PI*rhobar;
	v.cc = 2./3.*PI*rhobar;

	if(riso == 0.0)
	{
		v.press = 0.;
		v.psi1 += 4.0*PI*rho/3.0;
		v.theta1 += -4.0*PI*(rho+6.*x.press)/3.0;
		v.mass = 0.;
		v.mbar = 0.;
		v.mrest = 0.;
		//there probably ought to be some conditions here for v.qq, v.cc; see conditions for v.psi1, v.theta1
	}
	else
	{
		v.press = -(rho + x.press)*(x.theta1/x.theta0 - x.psi1/x.psi0);
		v.theta1 += -2.0*x.theta1/riso;
		v.psi1 += -2.0*x.psi1/riso;
		v.mass = 4.0*PI*SQR(riso)*rho*pow(x.psi0,6.)*sqrt(1.-2.*x.mass/riso/SQR(x.psi0));
		v.mbar = 4.0*PI*SQR(riso)*rhobar;
		v.qq += -2.*x.qq/riso;
		v.cc += -4.*x.cc/riso;
		v.mrest = 4.0*PI*SQR(riso)*rho0*pow(x.psi0,6.);
	}

	//   CCTK_VInfo(CCTK_THORNSTRING, "RHS out: %.15e %.15e %.15e %.15e %.15e %.15e\n", v.psi0, v.psi1, v.theta0, v.theta1, v.press, v.mass);
	return v;
}

void TOVStar::getRhoEpsOfPress(double press, double *rho0, double *eps, double *rho) const
{
	if(press == 0.)
	{
		*eps = 0.;
		*rho0 = 0.;
		*rho = 0.;
	}
	else
	{
		*rho0 = pow((press/tov_k), (1./tov_gamma));
		*eps = press/(*rho0)/(tov_gamma-1.);
		*rho = ((*rho0)*(1.+*eps));
	}
}

void TOVStar::addConformalRho(double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz) const
{
    DECLARE_CCTK_PARAMETERS;
    double distance = getDistanceFromCenter(x, y, z);
    rk4data_t interp = interpolateData(distance);

	double psi_this = interp.psi0;
    if(interp.press > 0.)
    {
		double lrho, leps, lrho0;
		getRhoEpsOfPress(interp.press, &lrho0, &leps, &lrho);
		double press = interp.press;

		//at some point: spin support?  this is not handled well.
		//

		double rhobar = (lrho+press)*pow(psi_this, 8.)/mbar;
		double sigmabar = rhobar*mbar/ibar;
		//ensure this is consistent with QCRhoPsi
		*dEbar = ((lrho+press)*SQR(w_linear)-press)*pow(psi_this, conformal_density_power);

		VectorAlgebra::Vector relative = VectorAlgebra::Vector(x, y, z) - position;
		VectorAlgebra::Vector Jbar = rhobar*momentum + sigmabar*crossProduct(angular_momentum, relative);
		//add some logic later for spin; this is nyi
		*dJbarx = Jbar.x;
		*dJbary = Jbar.y;
		*dJbarz = Jbar.z;
		*dlorentz = w_linear;
	}
	*dpsi += psi_this-1.;

	//CCTK_VInfo(CCTK_THORNSTRING, "addRho data at (%g, %g, %g) - radius %g: rho0 = %g, press = %g, eps = %g, psi = %g", x, y, z, distance, *drho0, *dpress, *deps, *dpsi_nonflat);
}

void TOVStar::addRho(double x, double y, double z, double *dpsi, double *drho0, double *deps, double *dpress, double *dv_x, double *dv_y, double *dv_z) const
{
    DECLARE_CCTK_PARAMETERS;
    double distance = getDistanceFromCenter(x, y, z);
    rk4data_t interp = interpolateData(distance);

    //this logic prevents ML from using combined psi as its initial guess for the conformal factor

    if(interp.press > 0.)
    {
		double lrho, leps, lrho0;
		getRhoEpsOfPress(interp.press, &lrho0, &leps, &lrho);
		*drho0 = lrho0;
		*deps = leps;
		*dpress = interp.press;

		double psi_this = interp.psi0;
		VectorAlgebra::Vector relative = VectorAlgebra::Vector(x,y,z) - position;
		VectorAlgebra::Vector linear_velocity = momentum/mbar/SQR(psi_this);
		VectorAlgebra::Vector velocity_from_angular = crossProduct(angular_momentum, relative)/ibar/SQR(psi_this);

		VectorAlgebra::Vector velocity = linear_velocity + velocity_from_angular;

		//at some point: spin support?  this is not handled well.
		*dv_x = velocity.x;
		*dv_y = velocity.y;
		*dv_z = velocity.z;
		*dpsi += interp.psi0 - 1.;
	}
	else if(!add_rhohats)
	{
		*dpsi += interp.psi0 - 1.;
	}

	//CCTK_VInfo(CCTK_THORNSTRING, "addRho data at (%g, %g, %g) - radius %g: rho0 = %g, press = %g, eps = %g, psi = %g", x, y, z, distance, *drho0, *dpress, *deps, *dpsi_nonflat);
}

void TOVStar::getQCNRhoPsi(double rr, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat) const
{
    rk4data_t interp = interpolateData(rr);
    double eps, rho0;
    getRhoEpsOfPress(interp.press, &rho0, &eps, rho);
    *qq = interp.qq;
    *cc = interp.cc;
	*nn = 4.*(*cc)*mbar/ibar;
    *psi_nonflat = interp.psi0 - 1.; 
    *rho = (*rho+interp.press)*SQR(w_linear)-interp.press;
}


bool TOVStar::isPointInsideStar(double x, double y, double z) const
{
    double distance = getDistanceFromCenter(x, y, z);
    return (distance <= radius);
}

double TOVStar::getDistanceFromCenter(double x, double y, double z) const
{
	VectorAlgebra::Vector relative = VectorAlgebra::Vector(x, y, z) - position;
    return relative.getMagnitude();
}

rk4data_t TOVStar::interpolateData(double distance) const
{
    rk4data_t interp;
    int last_point = inside_points - 1;
    int i_min, i_max;
    i_min = (int) (distance/dr);
    i_max = i_min + 1;
    if(i_max <= last_point)
    {
		// box desired value such that riso[i_min] <= riso < riso[i_min+1]
		//
		double driso, delta;

		driso = distance - tovdata[i_min].riso; // distance from desired value
		delta = tovdata[i_max].riso - tovdata[i_min].riso;  // step in riso 
		double minweight, maxweight;
		maxweight = (driso/delta);
		minweight = 1. - maxweight;
		interp = minweight*tovdata[i_min]+ maxweight*tovdata[i_max];
		if(std::isnan(interp.press))
		{
			CCTK_VInfo(CCTK_THORNSTRING, "i_min data at radius %g: riso = %g, press = %g, psi = %g, q = %g, c = %g", distance, tovdata[i_min].riso, tovdata[i_min].press, tovdata[i_min].psi0, tovdata[i_min].qq, tovdata[i_min].cc);
			CCTK_VInfo(CCTK_THORNSTRING, "i_max data at radius %g: riso = %g, press = %g, psi = %g, q = %g, c = %g", distance, tovdata[i_max].riso, tovdata[i_max].press, tovdata[i_max].psi0, tovdata[i_max].qq, tovdata[i_max].cc);
			CCTK_VInfo(CCTK_THORNSTRING, "linear interp data at radius %g: delta = %g; driso = %g", distance, delta, driso);
			CCTK_VInfo(CCTK_THORNSTRING, "point indices: last point is %d; i_min = %d; i_max = %d; radius/dr = %g; distance/dr = %g; rad-dist = %g", last_point, i_min, i_max, radius/dr, distance/dr, radius-distance);
		}
	}
    else
	{
		interp = tovdata[last_point];
		interp.press = 0.;
		interp.psi0 = 1. + (interp.psi0 - 1.)*(radius/distance);
		interp.theta0 = 1. - (1. - interp.theta0)*(radius/distance);
		interp.riso = distance;
		//qq, cc are automatically fine; psi1, theta1, mass0, mw2 should be regarded as unreliable
	}

//    CCTK_VInfo(CCTK_THORNSTRING, "Length of vector: %d. Inside points: %d", tovdata.size(), inside_points);
    return interp;
}

double TOVStar::getLapse(double x, double y, double z) const
{
    double distance = getDistanceFromCenter(x, y, z);
    rk4data_t interp = interpolateData(distance);
    return interp.theta0/interp.psi0;
} 
