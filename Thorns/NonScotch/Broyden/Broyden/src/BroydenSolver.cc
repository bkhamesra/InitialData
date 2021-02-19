#include "Broyden.hh"

namespace Broyden
{
	StepState BroydenSolver::takeFullStep(const StepState& data) const
	{
		return updater(step(data));
	}

	StepState BroydenSolver::solve(const BroydenSolver solver, const StepState& sinit)
	{
		//stop_condition is therefore S -> 2
		if(solver.stop_condition(sinit))
		{
			return sinit;
		}
		else
		{
			return solve(solver, solver.updater(step(sinit)));
		}
	}
	
	Vector2D BroydenSolver::inflateAndDeflate(const Vector2D& invec, const OutputFunction& outfunc)
	{
		return projectFNow(outfunc(inflate(invec)));
	}
		
	StepState BroydenSolver::makeState(const StencilDescription& desc, const OutputFunction& updater)
	{
		const auto do_inflate = [desc, updater](Vector2D (*xfunc)(const StencilDescription&))
		{
			return inflateAndDeflate(xfunc(desc), updater);
		};
		const auto fuplus = do_inflate(makeXUPlus);
		const auto fvplus = do_inflate(makeXVPlus);
		const auto fuminus = do_inflate(makeXUMinus);
		const auto fvminus = do_inflate(makeXVMinus);

		return makeFromStencil(desc, makeFStencil(fuplus, fuminus, fvplus, fvminus));

	}
		
	StepState BroydenSolver::solveAuto(const StencilDescription& desc, const OutputFunction& upd, const TestFunction& test)
	{
		return solve(BroydenSolver(test, upd),makeState(desc, upd));
	}
		
	StepState BroydenSolver::solveManual(const StepState& mss, const OutputFunction& upd, const TestFunction& test)
	{
		return solve(BroydenSolver(test, upd), mss);
	}
}
