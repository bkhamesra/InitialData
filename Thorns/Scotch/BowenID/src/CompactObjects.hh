#ifndef BOWEN_COMPACT_CLASSES
#define BOWEN_COMPACT_CLASSES

#include <VectorAlgebra.hh>
#include "Utilities.hh"
#include <vector>
#include <cmath>
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))

struct tovdata_t
{
    double press, psi0, theta0, psi1, theta1, riso, mbar, qq, cc, mass, nn, mass0;
    //rho: T(n) \cdot n
    //psi0: zeroth derivative of conformal factor
    //theta0: zeroth derivative of psi*alpha
    //psi1, theta1: first derivatives
    //riso: r isotropic (to avoid confusion with r_schwarzchild)
    //mbar: normalization factor used in various places.  momentum = mbar*W^2*v, etc.
    //qq: a curvature function
    //cc: another curvature function
    //mass0: part of rest mass not dependent on w^2
    //mw2: part of mass proportional to w^(-2)

    tovdata_t();
    tovdata_t(const tovdata_t& a);
    tovdata_t& operator=(const tovdata_t& a);
    tovdata_t& operator*=(double a);
    tovdata_t operator*(double a) const;
    tovdata_t& operator/=(double a);
    tovdata_t operator/(double a) const;
    tovdata_t& operator+=(const tovdata_t &a);
    tovdata_t operator+(const tovdata_t &a) const;
};

tovdata_t operator*(double a, const tovdata_t& b);

class Source
{
	public:
		virtual double getRhohat(double xx, double yy, double zz) = 0;
		virtual double getConfDensity(double xx, double yy, double zz) = 0;
		virtual double getConfPressure(double xx, double yy, double zz) = 0;
		virtual double getConfK(double xx, double yy, double zz) = 0;
		virtual VectorAlgebra::Vector getJhat(double xx, double yy, double zz) = 0;
		virtual double getNonFlatPsi(double xx, double yy, double zz) = 0;
		virtual double getInitialPsiGuess(double xx, double yy, double zz) = 0;
		virtual double getLorentzFac(double xx, double yy, double zz) = 0;
		virtual double getMass() { return 0.; };
		virtual double getRestMass() { return getMass(); };


		//Aij is *added* to the existing matrix
		virtual VectorAlgebra::LinearOperator getAij(double xx, double yy, double zz) = 0;

		virtual ~Source() = 0;
};

class CompactSource : public Source
{
	protected:
		VectorAlgebra::Vector position;
		int index;
	public:
		VectorAlgebra::Vector getRelativeVector(double xx, double yy, double zz) { VectorAlgebra::Vector input(xx, yy, zz); return input-position; };
		double getRelativeDistance(double xx, double yy, double zz) { return smoothR(getRelativeVector(xx, yy, zz).getMagnitude()); };

		CompactSource(int index);
		virtual ~CompactSource() = 0;
		
};

class EmptySource : public CompactSource
{
	public:
		double getRhohat(double xx, double yy, double zz) { return 0.0;} ;
		double getConfDensity(double xx, double yy, double zz) { return 0.0;} ;
		double getConfPressure(double xx, double yy, double zz) { return 0.0;} ;
		double getConfK(double xx, double yy, double zz) {return 0.0;};
		VectorAlgebra::Vector getJhat(double xx, double yy, double zz) { return VectorAlgebra::zero; }; 
		double getNonFlatPsi(double xx, double yy, double zz) { return 0.0;} ;
		double getInitialPsiGuess(double xx, double yy, double zz) { return 0.0; };
		double getLorentzFac(double xx, double yy, double zz) { return 0.0; };

		double getAdmMass() { return 0.0; };

		//Aij is *added* to the existing matrix
		void addAij(double xx, double yy, double zz, double Aij[3][3]) {};
		~EmptySource() {};
		EmptySource(int index) : CompactSource(index) {};
};

class BowenSource : public CompactSource
{
	protected:
		VectorAlgebra::Vector momentum;
		VectorAlgebra::Vector spin;
		double mass;
		virtual double getMomentumDist(double rr) = 0;
		virtual double getSpinDist(double rr) = 0;

		double getMomentumWeightDeriv(double rr) { return 4.*M_PI*getMomentumDist(rr)*SQR(rr); };
		double getPSDDeriv(double rr) { return 2./3.*M_PI*getMomentumDist(rr)*QUAD(rr); };
		double getSpinWeightDeriv(double rr) { return 8./3.*M_PI*getSpinDist(rr)*QUAD(rr); };

		virtual double getMomentumWeight(double rr) = 0; //Q
		virtual double getPointSourceDeviation(double rr) = 0; //C
		virtual double getSpinWeight(double rr) = 0; //N
		virtual double getNonFlatPsiOfR(double rr) = 0;
		
		double simpsonIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmid, double rmax, double fmin, double fmid, double fmax, int depth);
		double radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax);
		double radius;

	public:
		virtual ~BowenSource() = 0;
		BowenSource(int idx);
		double getNonFlatPsi(double xx, double yy, double zz);
		VectorAlgebra::LinearOperator getAij(double xx, double yy, double zz);
		VectorAlgebra::Vector getJhat(double xx, double yy, double zz);
		double getMass() { return mass; };
};


class BlackHole : public BowenSource
{
	public:
		BlackHole(int idx);
		double getInitialPsiGuess(double xx, double yy, double zz) { return getNonFlatPsi(xx, yy, zz); };
		double getLorentzFac(double xx, double yy, double zz) { return 0.; };
		~BlackHole() {};
		double getAdmMass();
		double getRhohat(double xx, double yy, double zz) { return 0.; };
		double getConfDensity(double xx, double yy, double zz) { return 0.0;} ;
		double getConfPressure(double xx, double yy, double zz) { return 0.0;} ;
		double getConfK(double xx, double yy, double zz) { return 0.0;} ;
	protected:
		double getMomentumDist(double rr) { return 0.; }; //sigma_p = 0 outside the puncture
		double getSpinDist(double rr) { return 0.; }; //sigma_j = 0 outside the puncture
		double getMomentumWeight(double rr) { return 1.; }; //Q
		double getPointSourceDeviation(double rr) { return 0.; }; //C
		double getSpinWeight(double rr) { return 1.; };  //N
		double getNonFlatPsiOfR(double rr) { return mass/2./rr; };
};

class SphericalHydroSource : public BowenSource
{
	public:
		SphericalHydroSource(int idx) : BowenSource(idx) {};
		virtual ~SphericalHydroSource() = 0;
		double getRhohat(double xx, double yy, double zz); //defined in terms of getRho, getPress, getWW
		double getConfDensity(double xx, double yy, double zz)  ;
		double getConfPressure(double xx, double yy, double zz)  ;
		double getInitialPsiGuess(double xx, double yy, double zz) { return 0.;} ;
		double getLorentzFac(double xx, double yy, double zz) {double W2 = getWW2(xx,yy,zz);  return sqrt(W2);} ;
		double getMW() { return mass*sqrt(getWW2(position.x, position.y, position.z)); };

	protected:
		double mbar;
		double ibar;
		virtual double getRho(double rr) = 0; //this is *not* rest-mass density; this is rho0*(1+eps)
		virtual double getPress(double rr) = 0;
		double getWW2(double xx, double yy, double zz);
		double getMomentumDist(double rr); //sigma_p
		double getSpinDist(double rr) { return getMomentumDist(rr)*mbar/ibar; }; //sigma_j
};

class TOVStar : public SphericalHydroSource
{
	public:
		TOVStar(int index);
		~TOVStar() { rk4data.clear(); };
		
	protected:
		int inside_points;
		double eos_gamma;
		double eos_k;
		double dr;
		double mass0;
		std::vector<tovdata_t> rk4data;
		
		double getRho(double rr) { return getRhoFromEOS(interpolateData(rr).press); };
		double getPress(double rr) { return interpolateData(rr).press; };
		double getConfK(double xx, double yy, double zz) ;
		double getNonFlatPsiOfR(double rr) { return interpolateData(rr).psi0 - 1.; };
		double getPressFromEOS(double rho0) { return eos_k*pow(rho0, eos_gamma); };
		double getRho0FromEOS(double press) { if(press > 0) { return std::pow(press/eos_k, 1./eos_gamma); } else { return 0.; } };
		double getRhoFromEOS(double press) { if(press > 0) { return getRho0FromEOS(press)+press/(eos_gamma-1.); } else { return 0.; } };
		double getMomentumWeight(double rr);
		double getPointSourceDeviation(double rr);
		double getSpinWeight(double rr);
		tovdata_t interpolateData(double rr);
		tovdata_t getRHS(double rr, tovdata_t xx);
		tovdata_t tov_rkstep(double rr, double dr, tovdata_t in);

		double getRestMass() { return mass0; }
};

//a pressureless cloud
class GaussianCloud : public SphericalHydroSource
{
	public:
		GaussianCloud(int idx);
		~GaussianCloud() {};

	protected:
		double amplitude; //radius = sigma sqrt(3)
		double getPress(double rr) { return 0.; }
		double getRho(double rr);
		double getConfK(double xx, double yy, double zz) {return 0.;};
		double getNonFlatPsiOfR(double rr) { if(rr < radius) { return amplitude*exp(-1.5*SQR(rr/radius)); } else { return exp(-1.5)*amplitude*radius/rr; } }

		//this example (and these functions) unfinished
		double getMomentumWeight(double rr) { return 1.; } 
		double getSpinWeight(double rr) { return 1.; }
		double getPointSourceDeviation(double rr) { return 0.; }
};

#endif
