#ifndef BROYDEN_MATH_HDR
#define BROYDEN_MATH_HDR

#include <cmath>
#include <functional>
#include <string>
#include <sstream>
#include "Vector2D.hh"
#include "cctk.h"

namespace Broyden
{
	struct PoststepState
	{
		Vector2D xnow, xpast, fpast;
		Matrix2D hpast;
		CCTK_INT iteration;
	};

	struct StepState
	{
		Vector2D xnow, xpast, fnow, fpast;
		Matrix2D hpast;
		CCTK_INT iteration;
	};

	struct StencilVectors
	{
		Vector2D uplus, uminus, vplus, vminus;
	};

	struct StencilDescription
	{
		Vector2D center, percents;
	};

	PoststepState step(const StepState& data);
	StepState makeFromStencil(const StencilDescription& desc, const StencilVectors& svf);
	Vector2D makeXUPlus(const StencilDescription& desc);
	Vector2D makeXVPlus(const StencilDescription& desc);
	Vector2D makeXUMinus(const StencilDescription& desc);
	Vector2D makeXVMinus(const StencilDescription& desc);
	Vector2D projectXNow(const PoststepState& ps);
	Vector2D projectFNow(const StepState& ss);
	 PoststepState inflate(const Vector2D& xn);
	StencilVectors makeFStencil(const Vector2D& fuplus, const Vector2D& fuminus, const Vector2D& fvplus, const Vector2D& fvminus);
	StepState makeStepState(const PoststepState& ps, const Vector2D& fnow);

	struct BroydenSolver
	{
		using TestFunction = std::function<bool(StepState)>;
	   	TestFunction stop_condition;

		using OutputFunction = std::function<StepState(PoststepState)>;
	   	OutputFunction updater;

		BroydenSolver(const TestFunction& test, const OutputFunction& out) : stop_condition(test), updater(out) {}

		StepState takeFullStep(const StepState& data) const;

		static StepState solve(const BroydenSolver solver, const StepState& init);

		static Vector2D inflateAndDeflate(const Vector2D& invec, const OutputFunction& outfunc);

		static StepState makeState(const StencilDescription& desc, const OutputFunction& updater);

		static StepState solveAuto(const StencilDescription& desc, const OutputFunction& upd, const TestFunction& test);

		static StepState solveManual(const StepState& mss, const OutputFunction& upd, const TestFunction& test);
	};
}

#endif
