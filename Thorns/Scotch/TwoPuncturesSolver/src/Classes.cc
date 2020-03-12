#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <assert.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"
#include "cctk.h"
#include "Classes.hh"
#include <new>
#include <ctime>
#include <cstdlib>

#define SQR(x) ((x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


Source::~Source() //required for destructors, even pure virtual ones
{
}

//PunctureSource definitions
PunctureSource::~PunctureSource()
{
}

PunctureSource::PunctureSource(int punc)
{
	DECLARE_CCTK_PARAMETERS;
	puncture = punc;
	r0[0] = (1.0-2.0*puncture)*par_b;
	r0[1] = 0.0;
	r0[2] = 0.0;
}

double SphericalHydroSource::getInitialGuess(double xx, double yy, double zz)
{
	return 0.;
}

SphericalSource::SphericalSource(int punc) : PunctureSource(punc)
{
	DECLARE_CCTK_PARAMETERS;
	pp[0] = par_Px[puncture];
	pp[1] = par_Py[puncture];
	pp[2] = par_Pz[puncture];
	ss[0] = par_Sx[puncture];
	ss[1] = par_Sy[puncture];
	ss[2] = par_Sz[puncture];
}


double PunctureSource::getRelativeDistance(double xx, double yy, double zz)
{
	double rr[3];
	getRelativeVector(xx, yy, zz, rr);
	return sqrt(scalarproduct(rr, rr, 3)); 
}

SphericalSource::~SphericalSource()
{
}

double SphericalSource::getNonFlatPsi(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx, yy, zz);
	return getNonFlatPsiOfR(rr);
}

void SphericalSource::addAij(double xx, double yy, double zz, double Aij[3][3])
{
	DECLARE_CCTK_PARAMETERS;
	double rr[3];
	getRelativeVector(xx, yy, zz, rr);
	double rmag = smoothR(sqrt(scalarproduct(rr, rr, 3)), epsilon);

	double ll[3];
	ll[0] = rr[0]/rmag;
	ll[1] = rr[1]/rmag;
	ll[2] = rr[2]/rmag;

	double junk, qq, cc, nn;
	getQCNRhoPsi(rmag, &qq, &cc, &nn, &junk, &junk);

	double l_cross_s[3];
	crossproduct(ll, ss, l_cross_s);

	double l_dot_p = scalarproduct(ll, pp, 3);
	double eta;
	
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(i == j)
			{
				eta = 1;
			}
			else
			{
				eta = 0.;
			}
			Aij[i][j] += (pp[i]*ll[j]+pp[j]*ll[i] + l_dot_p*(ll[i]*ll[j] - eta))*3.*qq/2./SQR(rmag);
			Aij[i][j] += (pp[i]*ll[j]+pp[j]*ll[i] + l_dot_p*(eta - 5*ll[i]*ll[j]))*3.*cc/QUAD(rmag);
			Aij[i][j] += -3.0*nn*(l_cross_s[i]*ll[j] + l_cross_s[j] * ll[i])/(rmag*rmag*rmag);
		}
	}
}

BlackHole::BlackHole(int punc) : SphericalSource(punc)
{
	DECLARE_CCTK_PARAMETERS;

	mass = par_m[puncture];
} 

/*
void BlackHole::addAij(double xx, double yy, double zz, double Aij[3][3])
{
	SphericalSource::addAij(xx, yy, zz, Aij);
	DECLARE_CCTK_PARAMETERS;
	double rmag = smoothR(getRelativeDistance(xx, yy, zz), epsilon);

	double ll[3];
	double rr[3];
	getRelativeVector(xx, yy, zz, rr);
	ll[0] = rr[0]/rmag;
	ll[1] = rr[1]/rmag;
	ll[2] = rr[2]/rmag;

	double l_cross_s[3];
	crossproduct(ll, ss, l_cross_s);

	//spin contribution
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			Aij[i][j] += -3.0*(l_cross_s[i]*ll[j] + l_cross_s[j] * ll[i])/(rmag*rmag*rmag);
		}
	}
}
*/

double BlackHole::addGaussian(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = smoothR(getRelativeDistance(xx, yy, zz), epsilon);
	double amplitude = 0.0;
        for ( int ng = 0; ng<cv_ngauss; ng++)
       	{
                 amplitude += cv_amp[ng]*exp(-0.5*SQR( (rr-cv_r0[ng])/cv_sigma[ng]) );
	}
	return amplitude;
}


SphericalHydroSource::SphericalHydroSource(int punc, void (*mlinfo)(int, double, double*, double*, double*, double*, double*), double (*mlmass)(int),  double (*mlradius)(int), int tov_offset) : SphericalSource(punc)
{
	DECLARE_CCTK_PARAMETERS;
	offset_number = puncture + tov_offset;
	getExternalQCNRhoPsi = mlinfo;
	mass = mlmass(offset_number);
	radius = mlradius(offset_number);
}
/*
double SphericalHydroSource::getNonFlatPsiOfR(double rr)
{
	double psi_nonflat, dummy;
	getQCRhoPsi(offset_number, rr, &dummy, &dummy, &dummy, &psi_nonflat);

	return psi_nonflat;
}
*/

double SphericalSource::getRho(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double rr = getRelativeDistance(xx,yy,zz);
	double rho, psi_nonflat, qq, cc, nn;
	getQCNRhoPsi(rr, &qq, &cc, &nn, &rho, &psi_nonflat);
	double psi = 1. + psi_nonflat;
//	return rho*pow(psi, 8.);
	return rho;
}

/*
double SphericalHydroSource::getQQ(double rr)
{
	double rho, psi_nonflat, qq, cc;
	getQCRhoPsi(offset_number, rr, &qq, &cc, &rho, &psi_nonflat);
	return qq;
}

double SphericalHydroSource::getCC(double rr)
{
	double rho, psi_nonflat, qq, cc;
	getQCRhoPsi(offset_number, rr, &qq, &cc, &rho, &psi_nonflat);
	return cc;
}
*/

double TPController::getAdmMassPlus(int ivar, int nvar, int n1, int n2, int n3, derivs v)
{
	return ps_plus->getAdmMass(ps_minus, ivar, nvar, n1, n2, n3, v);
}

double TPController::getAdmMassMinus(int ivar, int nvar, int n1, int n2, int n3, derivs v)
{
	return ps_minus->getAdmMass(ps_plus, ivar, nvar, n1, n2, n3, v);
}

double SphericalHydroSource::getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v)
{
	DECLARE_CCTK_PARAMETERS;
	return mass;
	/*
	int ppd = 100;
	double halfwidth = 2.*radius;
	double dh = 2.*halfwidth/ppd;
	double dA = dh*dh;
	double xmin = r0[0]-halfwidth;
	double ymin = r0[1]-halfwidth;
	double zmin = r0[2]-halfwidth;

	double xmax = xmin+2*halfwidth;
	double ymax = ymin+2*halfwidth;
	double zmax = zmin+2*halfwidth;

	double resultxy = 0.;
	double resultzx = 0.;
	double resultyz = 0.;
	for(int i = 0; i<ppd; i++)
	{
		for(int j = 0; j<ppd; j++)
		{
			double y1 = ymin+dh*i;
			double z1 = zmin+dh*j;

			double x2 = xmin+dh*i;
			double z2 = z1;

			double x3 = x2;
			double y3 = ymin+dh*j;

			resultxy += (-PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d3, x3, y3, zmin) + PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d3, x3, y3, zmax))*dA;
			resultyz += (-PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d1, xmin, y1, z1) + PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d1, xmax, y1, z1))*dA;
			resultxy += (-PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d2, x2, ymin, z2) + PunctIntPolForArbitVector(ivar, nvar, n1, n2, n3, v.d2, x2, ymax, z2))*dA;
		}
	}
	double result = -(resultxy + resultzx + resultyz)/(2.*Pi);
	return result;
	*/
	/*
	double u_here = ValueAtArbitPosition(ivar, nvar, n1, n2, n3, v, r0[0], r0[1], r0[2]);
	double psi_other_guess = other->getInitialGuess(r0[0], r0[1], r0[2]);
	double psi_here = 1 + psi_other_guess + u_here;
	double psi_isolated = 1 + getNonFlatPsi(r0[0], r0[1], r0[2]);
	double delta = psi_here - psi_isolated;
	double ratio = delta/psi_isolated;
	return mass*(1 - 3*ratio);
	*/
	/*
	srand(time(NULL));
	int points = 10000;
	double dV = 8*radius*radius*radius/points;
	double result = 0.;
	for(int i = 0; i < points; i++)
	{
		double xx = (2*radius*rand())/RAND_MAX + (r0[0]-radius);
		double yy = (2*radius*rand())/RAND_MAX + (r0[1]-radius);
		double zz = (2*radius*rand())/RAND_MAX + (r0[2]-radius);
		double u_here = ValueAtArbitPosition(ivar, nvar, n1, n2, n3, v, xx, yy, zz);
		double psi_other_guess = other->getInitialGuess(xx, yy, zz);
		double psi_here = 1 + psi_other_guess + u_here;
		double rhohat = getRhohat(xx, yy, zz);
		result += rhohat/psi_here/psi_here/psi_here*dV;
	}
	return result;
	*/
	/*
	int ppd = 50;
	double dh = 2*radius/ppd;
	double dV = dh*dh*dh;
	double staticresult = 0.;
	double traceresult = 0.;
	double fudge_factor = 1.005;

	for(int i = 0; i < ppd; i++)
	{
		for(int j = 0; j < ppd; j++)
		{
			for (int k = 0; k < ppd; k++)
			{
				double xx = r0[0]-radius+i*dh+dh/4.;
				double yy = r0[1]-radius+j*dh+dh/4.;
				double zz = r0[2]-radius+k*dh+dh/4.;
				double u_here = ValueAtArbitPosition(ivar, nvar, n1, n2, n3, v, xx, yy, zz);
				double psi_other_guess = other->getInitialGuess(xx, yy, zz);
				double psi_here = 1. + psi_other_guess + u_here;
				double rhohat = getRhohat(xx, yy, zz);
				staticresult += rhohat/psi_here/psi_here/psi_here*dV;
				double Aij[3][3];
				Aij[0][0] = 0.;
				Aij[0][1] = 0.;
				Aij[0][2] = 0.;
				Aij[1][0] = 0.;
				Aij[1][1] = 0.;
				Aij[1][2] = 0.;
				Aij[2][0] = 0.;
				Aij[2][1] = 0.;
				Aij[2][2] = 0.;
				addAij(xx, yy, zz, Aij);
				other->addAij(xx, yy, zz, Aij);
				double trace = 0. + 
				+	Aij[0][0]*Aij[0][0] 
				+	Aij[0][1]*Aij[0][1] 
				+	Aij[0][2]*Aij[0][2] 
				+	Aij[1][0]*Aij[1][0] 
				+	Aij[1][1]*Aij[1][1] 
				+	Aij[1][2]*Aij[1][2] 
				+	Aij[2][0]*Aij[2][0] 
				+	Aij[2][1]*Aij[2][1] 
				+	Aij[2][2]*Aij[2][2];
				if(trace > 1)
				{
					CCTK_VInfo(CCTK_THORNSTRING, "at (x,y,z) = (%g, %g, %g) [(i,j,k) = (%d, %d, %d)], trace = %g, psi = %g", xx, yy, zz,i, j, k, trace, psi_here);
				}
				traceresult += pow(psi_here, -7.)*trace/(16.*Pi)*dV;
			}
		}
	}
	staticresult *= fudge_factor;
	traceresult *= fudge_factor;
	double result = staticresult + traceresult;
	CCTK_VInfo(CCTK_THORNSTRING, "SHS mass info: static mass = %g; trace mass = %g", staticresult, traceresult);
	return result;
	*/
/*
	double error = .0001;
	double x0 = r0[0];
	double y0 = r0[1];
	double z0 = r0[2];

	//finding coordinates of boundary points
	double xp = x0 + radius;
	double xm = x0 - radius;
	double yp = y0 + radius;
	double ym = y0 - radius;
	double zp = z0 + radius;
	double zm = z0 - radius;

	double[3][3][3] values;
	for(int i = 0; i<5; i++)
	{
		for(int j = 0; j<5; j++)
		{
			for(int k = 0; k<5; k++)
			{
				values[i][j][k] = getAdmDensity(other, ivar, nvar, n1, n2, n3, v, x0+(i-2)*radius/2., y0+(j-2)*radius/2., z0+(k-2)*radius/2.);
			}
		}
	}

	simpson3d(values, x0, y0, z0, radius, other, params);	
	*/
}

/*
struct InterpParams
{
	int ivar, nvar, n1, n2, n3;
	derivs v;
} ;

double simpson3d(double[5][5][5] vals, double xmid, double ymid, double zmid, double dr, PunctureSource *other, InterpParmas params, double error)
{
	int coarseweights[5] = {1, 0, 4, 0, 1};
	int fineweights[5] = {1, 2, 4, 2, 1};
	double valcoarse = 0.;
	double valfine = 0.;

	for(int i = 0; i<5; i++)
	{
		for(int j = 0; j<5; j++)
		{
			for(int k = 0; k<5; k++)
			{
				valcoarse += coarseweights[i]*coarseweights[j]*coarseweights[k];
				valfine += fineweights[i]*fineweights[j]*fineweights[k];
			}
		}
	}
	if(abs(valcoarse - valfine) > error)
	{
		
	}
	else
	{
		return valfine;
	}
}

*/

double BlackHole::getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v)
{
	double u_here = ValueAtArbitPosition(ivar, nvar, n1, n2, n3, v, r0[0], r0[1], r0[2]);
	double psi_other_guess = other->getInitialGuess(r0[0], r0[1], r0[2]);
	double psi_nonsingular_here = 1 + psi_other_guess + u_here;
	return mass*(psi_nonsingular_here);
}

TPController::TPController()
{
	initializeSources();	
}

void TPController::deleteSources()
{
	delete ps_plus;
	delete ps_minus;
	std::vector<Source*>::iterator it;
	for(it = sources.begin(); it != sources.end(); ++it)
	{
		delete (*it);
	}
}

void TPController::reinitializeSources()
{
	deleteSources();
	initializeSources();
}

void TPController::initializeSources()
{
	ps_plus = PunctureSource::create(0);
	ps_minus = PunctureSource::create(1);
	//should be some logic here to instantiate other sources
}
	

TPController::~TPController()
{
	deleteSources();
}


TPController& TPController::create()
{
	static TPController controller;
	return controller;
}


void TPController::getAij(double xx, double yy, double zz, double Aij[3][3])
{
	Aij[0][0] = 0.;
	Aij[0][1] = 0.;
	Aij[0][2] = 0.;
	Aij[1][0] = 0.;
	Aij[1][1] = 0.;
	Aij[1][2] = 0.;
	Aij[2][0] = 0.;
	Aij[2][1] = 0.;
	Aij[2][2] = 0.;
	ps_plus->addAij(xx, yy, zz, Aij);
	ps_minus->addAij(xx, yy, zz, Aij);

	std::vector<Source*>::iterator it;
	for(it = sources.begin(); it != sources.end(); ++it)
	{
		(*it)->addAij(xx, yy, zz, Aij);
	}
}

double TPController::getNonFlatPsi(double xx, double yy, double zz)
{
	double psi_nf_plus = ps_plus->getNonFlatPsi(xx, yy, zz);
	double psi_nf_minus = ps_minus->getNonFlatPsi(xx, yy, zz);
	double psi_nonflat = psi_nf_plus + psi_nf_minus;
	std::vector<Source*>::iterator it;
	for(it = sources.begin(); it != sources.end(); ++it)
	{
		psi_nonflat += (*it)->getNonFlatPsi(xx, yy, zz);
	}
	return psi_nonflat;
}

double TPController::getRhohat(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	double power = conformal_density_power;
	double rhohat;
	if(add_rhohats)
	{
		double rhohat_plus = ps_plus->getRho(xx,yy,zz)*pow(1. + ps_plus->getNonFlatPsi(xx, yy, zz), power);
		double rhohat_minus = ps_minus->getRho(xx,yy,zz)*pow(1. + ps_minus->getNonFlatPsi(xx, yy, zz), power);

		std::vector<Source*>::iterator it;
		rhohat = rhohat_plus + rhohat_minus;
		for(it = sources.begin(); it != sources.end(); ++it)
		{
			rhohat += (*it)->getRho(xx, yy, zz)*pow(1.+(*it)->getNonFlatPsi(xx, yy, zz), power);
		}
	}
	else
	{
		double rho_plus = ps_plus->getRho(xx,yy,zz);
		double rho_minus = ps_minus->getRho(xx,yy,zz);
		
		double psi_nf_plus = ps_plus->getNonFlatPsi(xx, yy, zz);
		double psi_nf_minus = ps_minus->getNonFlatPsi(xx, yy, zz);

		std::vector<Source*>::iterator it;
		double rho = rho_plus + rho_minus;
		double psi = 1. + psi_nf_plus + psi_nf_minus;
		for(it = sources.begin(); it != sources.end(); ++it)
		{
			rho += (*it)->getRho(xx, yy, zz);
			psi += (*it)->getNonFlatPsi(xx, yy, zz);
		}
		rhohat = rho*pow(psi, power);
	}

	return rhohat;
}

//addGaussian is *not* defined for objects other than puncture sources
double TPController::addGaussian(double xx, double yy, double zz)
{
	return ps_plus->addGaussian(xx,yy,zz) + ps_minus->addGaussian(xx,yy,zz);
}

double TPController::getInitialGuess(double xx, double yy, double zz)
{
	double psi_nonflat = ps_plus->getInitialGuess(xx,yy,zz) + ps_minus->getInitialGuess(xx,yy,zz);
	std::vector<Source*>::iterator it;
	for(it = sources.begin(); it != sources.end(); ++it)
	{
		psi_nonflat += (*it)->getInitialGuess(xx, yy, zz);
	}
	return psi_nonflat;
}

double TPController::getAijAij(double xx, double yy, double zz)
{
	double Aij[3][3];

	getAij(xx, yy, zz, Aij);

	double trace = 0.0;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j<3; j++)
		{
			trace += Aij[i][j]*Aij[i][j];
		}
	}
	return trace;
}


double TPController::getR_plus(double xx, double yy, double zz)
{
	return ps_plus->getRelativeDistance(xx, yy, zz);
}

double TPController::getR_minus(double xx, double yy, double zz)
{
	return ps_minus->getRelativeDistance(xx, yy, zz);
}

void calcTestQCNRhoPsi(int puncture, double rr, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat)
{
	DECLARE_CCTK_PARAMETERS;
	double rshell = rr - cv_r0[puncture];
	double rshell2 = SQR(rshell);
	double sig2 = SQR(cv_sigma[puncture]);
	
	double chi2 = rshell2/sig2;
	double chi = sqrt(chi2);

	double gauschi = cv_amp[puncture]*exp(-0.5*chi2);
	double rh = gauschi*(1.0 + 2.0*(rshell/rr)-chi2 )/sig2/(2.0*Pi*pow(1.0+gauschi,5.0));
	psi_nonflat[0] = gauschi;
	rho[0] = rh;
	*qq = 0.;
	*cc = 0.;
	*nn = 0.;
}

double getTestMass(int puncture)
{
	return 0.0; //nyi
}

double getTestRadius(int puncture)
{
	DECLARE_CCTK_PARAMETERS;
	return cv_r0[puncture];
}

PunctureSource* PunctureSource::create(int punc)
{
	DECLARE_CCTK_PARAMETERS;

	if(CCTK_EQUALS(puncture_object[punc], "tov"))
	{
		int tov_offset=0;
		if(punc==1 && !CCTK_Equals(puncture_object[0], "tov"))
		{
			tov_offset=-1;
		}
		return new SphericalHydroSource(punc, calcTOVQCNRhoPsi, getTOVMass, getTOVRadius, tov_offset);
	}
	else if(CCTK_EQUALS(puncture_object[punc], "hydro test"))
	{
		return new SphericalHydroSource(punc, calcTestQCNRhoPsi, getTestMass, getTestRadius, 0);
	}

	else if(CCTK_EQUALS(puncture_object[punc], "none"))
	{
		return new EmptyPuncture(punc);
	}
	else
	{
		return new BlackHole(punc);
	}
}
