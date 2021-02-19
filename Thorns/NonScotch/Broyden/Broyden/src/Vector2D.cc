#include "Vector2D.hh"
#include <cmath>

namespace Broyden
{
	std::string to_string(const Vector2D& vec)
	{
		std::stringstream out;
		out << "(" << vec.u << "," << vec.v << ")";
		return out.str();
	}
	
	std::string to_string(const Matrix2D& mat)
	{
		std::stringstream out;
		out << "[[" << mat.m11 << "," << mat.m12 << "], " << "[" << mat.m21 << "," << mat.m22 << "]]";
		return out.str();
	}

	/*
	std::string to_string(const StepState& step)
	{
		std::stringstream out;
		out << "hpast = " << to_string(step.hpast) << " - xnow = " << to_string(step.xnow) << " - " << "fnow = " << to_string(step.fnow);
		return out.str();
	}
	*/

	std::ostream& operator<<(std::ostream& out, const Vector2D& vec)
	{
		out << to_string(vec);
		return out;
	}

	std::ostream& operator<<(std::ostream& out, const Matrix2D& mat)
	{
		out << to_string(mat);
		return out;
	}
	
	CCTK_REAL Vector2D::norm(const Vector2D& vv)
	{
		return sqrt(dotProduct(vv, vv));
	}
}
