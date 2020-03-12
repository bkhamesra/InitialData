#include "Broyden.hh"
#include <iostream>
namespace Broyden
{
	PoststepState step(const StepState& data)
	{
		const Vector2D fdiff = data.fnow - data.fpast;
		const Vector2D funit = fdiff/Vector2D::norm(fdiff);	
		const Matrix2D hnow = data.hpast + Matrix2D::outer(((data.xnow - data.xpast - Matrix2D::applyMatrix(data.hpast, fdiff))/Vector2D::dotProduct(fdiff, fdiff)), fdiff);
		const Vector2D xnext = data.xnow - Matrix2D::applyMatrix(hnow, data.fnow);

		const PoststepState post = { xnext, data.xnow, data.fnow, hnow, data.iteration+1 };
//		std::cout << "step iteration: " << data.iteration << " - xnow: " << data.xnow << "- fnow: " << data.fnow << " - hpast: " << data.hpast << std::endl;
		return post;
	}

	StepState makeFromStencil(const StencilDescription& desc, const StencilVectors& svf)
	{
		const auto uplus = desc.center.u*(1.+desc.percents.u);
		const auto uminus = desc.center.u*(1.-desc.percents.u);
		const auto vplus = desc.center.v*(1.+desc.percents.v);
		const auto vminus = desc.center.v*(1.-desc.percents.v);
		const auto dfdu = (svf.uplus - svf.uminus)/(uplus - uminus);
		const auto dfdv = (svf.vplus - svf.vminus)/(vplus - vminus);

		const auto jpast = Matrix2D(dfdu.u, dfdu.v, dfdv.u, dfdv.v);
		const auto hpast = Matrix2D::invert(jpast);

/*		std::cout << "uplus: " << uplus
			<< " - uminus: " << uminus
			<< " - vplus: " << vplus
			<< " - vminus: " << vminus
			<< " - fuplus: " << svf.uplus
			<< " - fuminus: " << svf.uminus
			<< " - fvplus: " << svf.vplus
			<< " - fvminus: " << svf.vminus
			<< " - dfdu: " << dfdu
			<< " - dfdv: " << dfdv
			<< " - jpast: " << jpast
			<< " - hpast: " << hpast
			<< std::endl;
			*/
		const StepState state = { Vector2D(uplus, desc.center.v), Vector2D(desc.center.u, vplus), svf.uplus, svf.vplus, hpast, 0 };
		return state;
	}

	Vector2D makeXUPlus(const StencilDescription& desc)
	{
		return Vector2D(desc.center.u*(1.+desc.percents.u), desc.center.v);
	}

	Vector2D makeXVPlus(const StencilDescription& desc)
	{
		const auto xvplus = Vector2D(desc.center.u, desc.center.v*(1.+desc.percents.v));
//		std::cout << "in makeXVPlus: xvplus: " << xvplus << std::endl;
		return xvplus;
	}
	
	Vector2D makeXUMinus(const StencilDescription& desc)
	{
		return Vector2D(desc.center.u*(1.-desc.percents.u), desc.center.v);
	}

	Vector2D makeXVMinus(const StencilDescription& desc)
	{
		const auto xvminus = Vector2D(desc.center.u, desc.center.v*(1.-desc.percents.v));
//		std::cout << "in makeXVPlus: xvminus: " << xvminus << std::endl;
		return xvminus;
	}

	Vector2D projectXNow(const PoststepState& ps)
	{
		return ps.xnow;
	}

	StencilVectors makeFStencil(const Vector2D& fuplus, const Vector2D& fuminus, const Vector2D& fvplus, const Vector2D& fvminus)
	{
		const StencilVectors svf = { fuplus, fuminus, fvplus, fvminus };
		/*
		std::cout << "in makeFstencil:"
			<< " - fuplus: " << svf.uplus
			<< " - fuminus: " << svf.uminus
			<< " - fvplus: " << svf.vplus
			<< " - fvminus: " << svf.vminus
			<< std::endl;
			*/
		return svf;
	}

	StepState makeStepState(const PoststepState& ps, const Vector2D& fnow)
	{
		const StepState state = { ps.xnow, ps.xpast, fnow, ps.fpast, ps.hpast, ps.iteration };
		return state;
	}

	PoststepState inflate(const Vector2D& xnow)
	{
		const PoststepState ps = { xnow, Vector2D(), Vector2D(), Matrix2D(), 0 };
		return ps;
	}

	Vector2D projectFNow(const StepState& ss)
	{
		return ss.fnow;
	}

}


