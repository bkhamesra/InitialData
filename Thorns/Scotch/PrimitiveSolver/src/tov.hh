#include <stdio.h>
#include <vector>
#include <VectorAlgebra.hh>

struct rk4data_t
{
    double press, psi0, theta0, psi1, theta1, riso, mbar, qq, cc, mass, mrest;
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

    rk4data_t() {press=psi0=theta0=psi1=theta1=mbar=mass=qq=cc=riso=mrest=0.0;}
    rk4data_t(const rk4data_t& a) {press = a.press;psi0 = a.psi0;theta0 = a.theta0;psi1 = a.psi1;theta1 = a.theta1; riso = a.riso; mbar = a.mbar;qq = a.qq; cc = a.cc; mass = a.mass; mrest=a.mrest; }

    rk4data_t& operator=(const rk4data_t& a) { press = a.press; psi0 = a.psi0; theta0 = a.theta0; psi1 = a.psi1; theta1 = a.theta1; riso = a.riso; mbar = a.mbar; qq = a.qq; cc = a.cc; mass = a.mass; mrest = a.mrest; return *this; }

    rk4data_t& operator*=(double a) {press *= a;psi0 *= a; psi1 *= a; theta0 *= a;theta1*=a; riso *= a; mbar *= a; qq *= a; cc *= a; mass *= a; mrest *= a; return *this;}

    rk4data_t operator*(double a) const {return rk4data_t(*this).operator*=(a);}

    rk4data_t& operator/=(double a) {press /= a;psi0 /= a;theta0 /= a;psi1 /= a; theta1 /= a; riso /= a; mbar /= a; qq /= a; cc /= a; mass /= a; mrest /= a;  return *this;}

    rk4data_t operator/(double a) const {return rk4data_t(*this).operator/=(a);}

    rk4data_t& operator+=(const rk4data_t &a) {press += a.press;psi0 += a.psi0;theta0 += a.theta0;psi1 += a.psi1; theta1 += a.theta1; riso += a.riso; mbar += a.mbar; qq += a.qq; cc += a.cc; mass += a.mass; mrest += a.mrest; return *this;}

    rk4data_t operator+(const rk4data_t &a) const {return rk4data_t(*this).operator+=(a);}
};

rk4data_t operator*(double a, const rk4data_t& b) {return rk4data_t(b).operator*=(a); }

class TOVStar
{
	public:
		TOVStar(int index, FILE *testsuitefh);
		~TOVStar() {};

		void addRho(double x, double y, double z, double *dpsi_nonflat, double *drho0, double *deps, double *dpress, double *dv_x, double *dv_y, double *dv_z) const;
		void addConformalRho(double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz) const;
		void getQCNRhoPsi(double rr, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat) const;

		double getLapse(double x, double y, double z) const;
		double getMass() { return mass*w_linear; };
		double getRadius() { return radius; };
		bool isPointInsideStar(double x, double y, double z) const;


	protected:
	        rk4data_t rkstep(double x, double h, rk4data_t y_in, rk4data_t (TOVStar::*f)(double, rk4data_t) const);

		double getDistanceFromCenter(double x, double y, double z) const;
		void getRhoEpsOfPress(double press, double *rho, double *eps, double *rho0) const;

		rk4data_t interpolateData(double distance) const;

		rk4data_t TOV_RHS(double riso, rk4data_t x) const;

		int inside_points;
		double tov_gamma;
		double tov_k;
		double rho_c;
		double rho_atmosphere;
		double dr;
		double radius;
		double mass;
		double mbar;
		double ibar;  //4*Cmax*Mbar
		double mrest;
		int star_index;
		int numpoints;
		std::vector<rk4data_t> tovdata;
//		double velocity[3];
		VectorAlgebra::Vector momentum;
		VectorAlgebra::Vector angular_momentum;
		VectorAlgebra::Vector position;
		double w_linear;
};

class TOVController
{
	
	public:
		static TOVController& create(int count);
		TOVStar getStar(int index);
		void getQCNRhoPsi(int index, double distance, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat);
		double getLapse(double x, double y, double z);
		void getRho(double x, double y, double z, double *dpsi, double *drho0, double *deps, double *dpress, double *dv_x, double *dv_y, double *dv_z);
		void getConformalRho(double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz);

		void updateTOVs();
		double getMass(int index);
		double getRadius(int index);

	private:
		int numstars;
		std::vector<TOVStar> stars;
		TOVController(int count);
		~TOVController();
		TOVController& operator=(TOVController const&){};
};
