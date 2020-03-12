#ifndef VECTOR_2D_HDR
#define VECTOR_2D_HDR

#include "cctk.h"
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

namespace Broyden
{
	struct Vector2D
	{
		CCTK_REAL u, v;
		constexpr Vector2D(const CCTK_REAL& in1, const CCTK_REAL& in2): u(in1), v(in2) {}
		constexpr Vector2D() : u(0.), v(0.) {}
		constexpr static CCTK_REAL dotProduct(const Vector2D& vv, const Vector2D& ww);
		static CCTK_REAL norm(const Vector2D& vv);
	};

	struct Matrix2D
	{
		CCTK_REAL m11, m12, m21, m22;
		constexpr Matrix2D(const CCTK_REAL& in11, const CCTK_REAL& in12, const CCTK_REAL& in21, const CCTK_REAL& in22) : m11(in11), m12(in12), m21(in21), m22(in22) {}
		constexpr Matrix2D() : m11(0.), m12(0.), m21(0.), m22(0.) {}
		constexpr static Vector2D applyMatrix(const Matrix2D& aa, const Vector2D& vv);
		constexpr static Matrix2D outer(const Vector2D& vv, const Vector2D& ww);
		constexpr static Matrix2D invert(const Matrix2D& mat);
	};


	constexpr Matrix2D operator+(const Matrix2D& aa, const Matrix2D& bb);
	constexpr Matrix2D operator-(const Matrix2D& aa, const Matrix2D& bb);
	constexpr Matrix2D operator-(const Matrix2D& aa);
	constexpr Matrix2D operator*(const Matrix2D& aa, const CCTK_REAL alpha);
	constexpr Matrix2D operator*(const CCTK_REAL alpha, const Matrix2D& aa);
	constexpr Matrix2D operator/(const Matrix2D& aa, const CCTK_REAL alpha);

	constexpr Vector2D operator+(const Vector2D& vv, const Vector2D& ww);
	constexpr Vector2D operator-(const Vector2D& vv, const Vector2D& ww);
	constexpr Vector2D operator-(const Vector2D& aa);
	constexpr Vector2D operator*(const Vector2D& vv, const CCTK_REAL alpha);
	constexpr Vector2D operator*(const CCTK_REAL alpha, const Vector2D& vv);
	constexpr Vector2D operator/(const Vector2D& vv, const CCTK_REAL alpha);
	
	
	constexpr Matrix2D operator+(const Matrix2D& aa, const Matrix2D& bb)
	{
		return Matrix2D(
				aa.m11+bb.m11,
				aa.m12+bb.m12,
				aa.m21+bb.m21,
				aa.m22+bb.m22
				);
	}

   	constexpr Matrix2D operator-(const Matrix2D& aa)
	{
		return aa*(-1.);
	}

   	constexpr Matrix2D operator-(const Matrix2D& aa, const Matrix2D& bb)
	{
		return aa+(-bb);
	}

	constexpr Matrix2D operator*(const Matrix2D& aa, const CCTK_REAL alpha)
	{
		return Matrix2D( aa.m11*alpha, aa.m12*alpha, aa.m21*alpha, aa.m22*alpha);
	}

	constexpr Matrix2D operator*(const CCTK_REAL alpha, const Matrix2D& aa)
	{
		return aa*alpha;
	}

	constexpr Matrix2D operator/(const Matrix2D& aa, const CCTK_REAL alpha)
	{
		return aa*(1./alpha);
	}

	constexpr Vector2D operator+(const Vector2D& vv, const Vector2D& ww)
	{
		return Vector2D
			(
				vv.u+ww.u,
				vv.v+ww.v
			);
	}

	constexpr Vector2D operator-(const Vector2D& vv)
	{
		return vv*(-1.);
	}

	constexpr Vector2D operator-(const Vector2D& vv, const Vector2D& ww)
	{
		return vv+(-ww);
	}

	constexpr Vector2D operator*(const Vector2D& vv, const CCTK_REAL alpha)
	{
		return Vector2D
			(
				vv.u*alpha,
				vv.v*alpha 
			);
	}

	constexpr Vector2D operator*(const CCTK_REAL alpha, const Vector2D& vv)
	{
		return vv*alpha;
	}

	constexpr Vector2D operator/(const Vector2D& vv, const CCTK_REAL alpha)
	{
		return vv*(1./alpha);
	}

	constexpr Vector2D Matrix2D::applyMatrix(const Matrix2D& aa, const Vector2D& vv)
	{
		return Vector2D ( aa.m11*vv.u + aa.m12*vv.v, aa.m21*vv.u + aa.m22*vv.v);
	}

	constexpr Matrix2D Matrix2D::outer(const Vector2D& vv, const Vector2D& ww)
	{
		return Matrix2D
			(
			 vv.u*ww.u,
			 vv.u*ww.v,
			 vv.v*ww.u,
			 vv.v*ww.v 
			);
	}

	constexpr CCTK_REAL Vector2D::dotProduct(const Vector2D& vv, const Vector2D& ww)
	{
		return vv.u*ww.u+vv.v*ww.v;
	}

	constexpr Matrix2D Matrix2D::invert(const Matrix2D& mat)
	{
		return Matrix2D(mat.m11, -mat.m21, -mat.m12, mat.m22)/(mat.m11*mat.m22 - mat.m12*mat.m21);
	}
	std::string to_string(const Vector2D& vec);
	std::string to_string(const Matrix2D& mat);

	std::ostream& operator<<(std::ostream& out, const Vector2D& vec);
	std::ostream& operator<<(std::ostream& out, const Matrix2D& mat);
}

#endif
