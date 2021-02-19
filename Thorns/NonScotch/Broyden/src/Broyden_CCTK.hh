#include <functional>
#include <string>
#include "Broyden.hh"
#include "cctk.h"


//D: StencilDescription
//V: Vector2D
//2: boolean
//S: StepState
struct ObjectDesc
{
	std::string object_type;
	CCTK_REAL xloc, yloc, zloc, target, integration_radius;
	CCTK_INT numpoints, polyhandle, object_index;

	//Required by getMassFunction
	void setParameter(const CCTK_REAL& setval) const;
	CCTK_REAL getMass(const CCTK_REAL& component) const;
};

using ObjectPair = std::pair<ObjectDesc, ObjectDesc>;

extern "C"
void Broyden_Run(CCTK_ARGUMENTS);

//Core functions for Broyden_Solve
bool getTestParameter();
void testOrSolve(const cGH* gh, const bool& runtest);
void testIntegrator(const cGH* gh);
void solveBinaryBroyden(const cGH* gh);

//Core functions for solveBinaryBroyden
Broyden::StencilDescription getStencilParams();
ObjectPair makeGhObjectPair(const cGH* gh);
Broyden::StepState getOutputFunction(const ObjectPair& opair, const Broyden::PoststepState& bowenpost);
Broyden::BroydenSolver::TestFunction checkStopConditions();

CCTK_REAL getMassValue(const CCTK_REAL& component, const ObjectDesc& ob);

ObjectPair makeObjectPair(const ObjectDesc& ob1, const ObjectDesc& ob2);
CCTK_INT getPolytropeHandle(const cGH* gh);
ObjectDesc makeObjectDescFromIndex(const CCTK_INT& idx, const CCTK_INT& polyhandle);

//Functions required by getOutputFunction
Broyden::StepState getMass(const ObjectPair& opair, const Broyden::PoststepState& bowenpost);
void setParameters(const ObjectPair& opair, const Broyden::Vector2D& xnow);
void recomputeSources();
Broyden::Vector2D makeVector(const CCTK_REAL& c1, const CCTK_REAL& c2);
Broyden::Vector2D getMassVector(const ObjectPair& opair, const Broyden::Vector2D& xnow);

//Utility, common functions
void printLn(const std::string& msg);
