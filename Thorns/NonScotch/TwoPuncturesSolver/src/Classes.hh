#ifndef TPO_CLASSES
#define TPO_CLASSES

#include <vector>
#define SQR(x) ((x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))
#define CUBE(x) ((x)*SQR(x))

class Source
{
	public:
		virtual double getRho(double xx, double yy, double zz) = 0;
		virtual double getNonFlatPsi(double xx, double yy, double zz) = 0;
		virtual double getInitialGuess(double xx, double yy, double zz) = 0;


		//Aij is *added* to the existing matrix
		virtual void addAij(double xx, double yy, double zz, double Aij[3][3]) = 0;

		virtual ~Source() = 0;
};

class PunctureSource : public Source
{
	protected:
		double r0[3];
		int puncture;
	public:
		static PunctureSource* create(int punc);
		void getRelativeVector(double xx, double yy, double zz, double rr[3]) { rr[0] = xx - r0[0]; rr[1] = yy - r0[1]; rr[2] = zz - r0[2]; };
		void getPosition(double *xx, double *yy, double *zz) { *xx = r0[0]; *yy = r0[1]; *zz = r0[2]; };
		double getRelativeDistance(double xx, double yy, double zz);
		virtual double addGaussian(double xx, double yy, double zz) = 0;
		virtual double getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v) = 0;

		PunctureSource(int punc);
		virtual ~PunctureSource() = 0;
		
};

class EmptyPuncture : public PunctureSource
{
	public:
		double getRho(double xx, double yy, double zz) { return 0.0;} ;
		double getNonFlatPsi(double xx, double yy, double zz) { return 0.0;} ;
		double getInitialGuess(double xx, double yy, double zz) { return 0.0; };

		double addGaussian(double xx, double yy, double zz) { return 0.0;};
		double getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v) { return 0.0; };

		//Aij is *added* to the existing matrix
		void addAij(double xx, double yy, double zz, double Aij[3][3]) {};
		~EmptyPuncture() {};
		EmptyPuncture(int punc) : PunctureSource(punc) {};
};

class SphericalSource : public PunctureSource
{
	protected:
		double pp[3];
		double ss[3];
		double mass;
		virtual void getQCNRhoPsi(double rr, double *qq, double *cc, double *nn, double *rho, double* psi) { *qq = 0.; *cc = 0.; *nn = 0.; *rho = 0.; *psi = 0.; };
		double getNonFlatPsiOfR(double rr) { double junk, psi; getQCNRhoPsi(rr, &junk, &junk, &junk, &junk, &psi); return psi; };
	public:
		double getRho(double xx, double yy, double zz);
		virtual ~SphericalSource() = 0;
		SphericalSource(int punc);
		double getNonFlatPsi(double xx, double yy, double zz);
		void addAij(double xx, double yy, double zz, double Aij[3][3]);
		virtual double getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v) = 0;
};


class BlackHole : public SphericalSource
{
	public:
		BlackHole(int punc);
		double getInitialGuess(double xx, double yy, double zz) { return getNonFlatPsi(xx, yy, zz); };
		double addGaussian(double xx, double yy, double zz);
		~BlackHole() {};
		double getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v);
	protected:
		void getQCNRhoPsi(double rr, double* qq, double* cc, double* nn, double* rho, double* psi) { *qq = 1.; *cc = 0.; *nn = 1.; *rho = 0.; *psi = mass/(2.*rr); };
};

class SphericalHydroSource : public SphericalSource
{
	public:
		SphericalHydroSource(int punc, void (*mlinfo)(int, double, double*, double*, double*, double*, double*), double (*mlmass)(int), double (*mlradius)(int), int tov_offset);
		double getInitialGuess(double xx, double yy, double zz);
		~SphericalHydroSource() {};
		double addGaussian(double xx, double yy, double zz) { return 0.0;};
		double getAdmMass(PunctureSource *other, int ivar, int nvar, int n1, int n2, int n3, derivs v);

	protected:
		int offset_number;
		void (*getExternalQCNRhoPsi)(int, double, double*, double*, double*, double*, double*);
		void getQCNRhoPsi(double rr, double* qq, double* cc, double* nn, double* rho, double* psi) { getExternalQCNRhoPsi(offset_number, rr, qq, cc, nn, rho, psi); };
		double radius;
};

class TPController
{
	private:
		PunctureSource *ps_plus, *ps_minus;
		std::vector<Source*> sources;
		TPController();
		~TPController();
		TPController& operator=(TPController const&){};
		void initializeSources();
		void deleteSources();
	public:
		void reinitializeSources();
		double getR_plus(double xx, double yy, double zz);
		double getR_minus(double xx, double yy, double zz);
		double getAdmMassPlus(int ivar, int nvar, int n1, int n2, int n3, derivs v);
		double getAdmMassMinus(int ivar, int nvar, int n1, int n2, int n3, derivs v);
		static TPController& create();
		void getAij(double xx, double yy, double zz, double Aij[3][3]);
		double getAijAij(double xx, double yy, double zz);
		double getInitialGuess(double xx, double yy, double zz);
		
		//trivial inline stuff
		double getNonFlatPsi(double xx, double yy, double zz);
		double getRhohat(double xx, double yy, double zz);

		double addGaussian(double xx, double yy, double zz);
		
};
#endif
