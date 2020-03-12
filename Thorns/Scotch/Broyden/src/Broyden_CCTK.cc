#include <functional>
#include <iostream>
#include <sstream>
#include "Broyden_CCTK.hh"
#include "BasicIntegrator.hh"
#include "ParameterTools.hh"
#include "Broyden.hh"
#include "cctk.h"
#include "cctk_Parameters.h"

constexpr CCTK_REAL powi(CCTK_REAL base, CCTK_INT radix)
{
	return radix == 0 ? 1. : base*powi(base, radix-1);
}
using namespace Broyden;

extern "C"
void Broyden_Run(CCTK_ARGUMENTS)
{
	testOrSolve(cctkGH, getTestParameter());
}

bool getTestParameter()
{
	DECLARE_CCTK_PARAMETERS;
	return test_integrator;
}

void testOrSolve(const cGH* gh, const bool& runtest)
{
	return runtest ? testIntegrator(gh) : solveBinaryBroyden(gh);
}

void testIntegrator(const cGH* gh)
{
	//not expected to get the right result for a bh, but we're not worried about that
	recomputeSources();
	getMassValue(0., makeObjectDescFromIndex(0, getPolytropeHandle(gh)));
}

void solveBinaryBroyden(const cGH* gh)
{
	const auto sparams = getStencilParams();
	const auto opair = makeGhObjectPair(gh);
	const BroydenSolver::OutputFunction outputfunc = [opair](const PoststepState& ps) { return getMass(opair, ps); };
	const BroydenSolver::TestFunction testfunc = checkStopConditions();
	BroydenSolver::solveAuto(sparams, outputfunc, testfunc);
}

//Functions required only by testIntegrator
CCTK_REAL getMassValue(const CCTK_REAL& component, const ObjectDesc& ob)
{
	return ob.getMass(component);
}

//Functions required by solveBinaryBroyden
StencilDescription getStencilParams()
{
	DECLARE_CCTK_PARAMETERS;
	const StencilDescription desc = { Vector2D(initial_guess[0], initial_guess[1]), Vector2D(range[1], range[0]) };
	return desc;
}

ObjectPair makeGhObjectPair(const cGH* gh)
{
	const auto polyhandle = getPolytropeHandle(gh);
	return makeObjectPair(makeObjectDescFromIndex(0, polyhandle), makeObjectDescFromIndex(1, polyhandle));
}

StepState getMass(const ObjectPair& opair, const PoststepState& bowenpost)
{
	std::stringstream report;
	report << "iteration #" << bowenpost.iteration << " - xnow: " << to_string(bowenpost.xnow) << " - fpast: " << bowenpost.fpast << " - hpast: " << to_string(bowenpost.hpast) << std::endl;

	CCTK_INFO(report.str().c_str());
	setParameters(opair, bowenpost.xnow);
	recomputeSources();
	return makeStepState(bowenpost, getMassVector(opair, bowenpost.xnow));

}

BroydenSolver::TestFunction checkStopConditions()
{
	const std::function<bool(StepState)> test = [](StepState state) -> bool
	{
		DECLARE_CCTK_PARAMETERS;
		return state.iteration >= max_iterations || (std::abs(state.fnow.u) < tolerance[0] && std::abs(state.fnow.v) < tolerance[1]);
	};	
	return test;
}

ObjectPair makeObjectPair(const ObjectDesc& ob1, const ObjectDesc& ob2)
{
	return std::make_pair(ob1, ob2);
}

CCTK_INT getPolytropeHandle(const cGH* gh)
{
	return CCTK_IsThornActive("Whisky") ? *(CCTK_INT*)CCTK_VarDataPtr(gh, 0, "whisky::whisky_polytrope_handle") : -1;
}

ObjectDesc makeObjectDescFromIndex(const CCTK_INT& idx, const CCTK_INT& polyhandle)
{
		DECLARE_CCTK_PARAMETERS;
		if(idx >= 0)
		{
			const ObjectDesc object =  { compact_object[idx], object_rx[idx], object_ry[idx], object_rz[idx], target[idx], integration_radius[idx], integration_npoints[idx], polyhandle, idx };
			return object;
		}
		else
		{
			CCTK_WARN(CCTK_WARN_ABORT, "object index is invalid!");
		}
}

//Functions required by getMassFunction
void ObjectDesc::setParameter(const CCTK_REAL& setval) const
{

	//shoudl clean this up
	int paramset_ierr;
	if(CCTK_Equals(this->object_type.c_str(), "black hole"))
	{
		std::stringstream paramname;
		paramname << "bh_bare_mass[" << this->object_index << "]";
		std::stringstream vbout;
		vbout << "Setting bowenid::" << paramname.str() << " = " << setval;
		CCTK_INFO(vbout.str().c_str());
		paramset_ierr = ParameterTools::setRealParameter(paramname.str().c_str(), "bowenid", setval);
	}
	else if(CCTK_Equals(this->object_type.c_str(), "tov"))
	{
		std::stringstream paramname;
		paramname << "tov_rho_central[" << this->object_index << "]";
		std::stringstream vbout;
		vbout << "Setting bowenid::" << paramname.str() << " = " << setval;
		CCTK_VInfo(CCTK_THORNSTRING, vbout.str().c_str());
		paramset_ierr = ParameterTools::setRealParameter(paramname.str().c_str(), "bowenid",  setval);
	}
	else
	{
		std::stringstream err;
		err << "no object detected for object #" << this->object_index;
		CCTK_WARN(CCTK_WARN_ABORT, err.str().c_str());
	}
	if(paramset_ierr < 0)
		CCTK_WARN(CCTK_WARN_ABORT, "could not set parameter!");
}

CCTK_REAL ObjectDesc::getMass(const CCTK_REAL& component) const
{
	if(CCTK_Equals(this->object_type.c_str(), "black hole"))
	{
		std::stringstream vbout;
		const auto retval = ((1.+getU(this->xloc, this->yloc, this->zloc) + getExternalPsiGuess(this->object_index, this->xloc, this->yloc, this->zloc))*component);
		vbout << "Found mass of black hole #" << this->object_index << " = " << retval;
		CCTK_INFO(vbout.str().c_str());
		return (retval-this->target);
	}
	else if(CCTK_Equals(this->object_type.c_str(), "tov"))
	{
		if(this->polyhandle < 0)
		{
			return CCTK_WARN(CCTK_WARN_ABORT, "no valid polytrope handle for tov integration");
		}
		else
		{
			const auto handle = this->polyhandle;
			const std::function<CCTK_REAL(CCTK_REAL, CCTK_REAL, CCTK_REAL)> integrand = [handle](CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) -> CCTK_REAL
			{
				const auto cf = getConformalFactor(x, y, z);
				return getRho0W(handle, x, y, z, cf)*powi(cf, 6);
			};
			const auto csx = BasicIntegrator::CoordSpec(this->xloc-this->integration_radius, this->xloc+this->integration_radius, this->numpoints);
			const auto csy = BasicIntegrator::CoordSpec(this->yloc-this->integration_radius, this->yloc+this->integration_radius, this->numpoints);
			const auto csz = BasicIntegrator::CoordSpec(this->zloc-this->integration_radius, this->zloc+this->integration_radius, this->numpoints);
			const auto iparams = BasicIntegrator::IntegralParams(csx, csy, csz);
			std::stringstream pre_integral;
			pre_integral << "Performing basic integration for TOVStar #" << this->object_index;
			CCTK_INFO(pre_integral.str().c_str());
			std::stringstream post_integral;
			const auto iresult = (BasicIntegrator::integrate(iparams, integrand))*getSourceMass(this->object_index)/getRestMass(this->object_index);
			post_integral << "TOV Mass: " << iresult;
			CCTK_VInfo(CCTK_THORNSTRING, post_integral.str().c_str());
			return (iresult - this->target);
		}
	}
	else
	{
		std::stringstream err;
		err << "no object detected for object #" << this->object_index;
		CCTK_WARN(CCTK_WARN_ABORT, err.str().c_str());
	}
}

void setParameters(const ObjectPair& opair, const Vector2D& bpars)
{
	std::get<0>(opair).setParameter(bpars.u);
	std::get<1>(opair).setParameter(bpars.v);
}

void recomputeSources()
{
	remakeSources(1);
	recomputeTP(1);
}

Vector2D makeVector(const CCTK_REAL& c1, const CCTK_REAL& c2)
{
	return Vector2D(c1, c2);
}

Vector2D getMassVector(const ObjectPair& opair, const Vector2D& xnow)
{
	return makeVector(std::get<0>(opair).getMass(xnow.u), std::get<1>(opair).getMass(xnow.v));
}

StepState makeStepState(const PoststepState& bowenpost, const Vector2D& fnow)
{
	const StepState res = { bowenpost.xnow, bowenpost.xpast, fnow, bowenpost.fpast, bowenpost.hpast, bowenpost.iteration };
	return res;
}
