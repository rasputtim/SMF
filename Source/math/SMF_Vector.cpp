/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the freeMem Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the freeMem Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#include "geometry/SMF_2DPoint.h"
#include "math/SMF_Vector.h"
#include "util/SMF_StringUtils.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_EulerAngles.h"
#include "math/SMF_Math.h"
#include "geometry/all.h"

namespace SMF {
namespace MATH{
const CVec2D CVec2D::origin( 0.0f, 0.0f );
const CVec2D CVec2D::one = CVec2D(1, 1);
const CVec2D CVec2D::unitX = CVec2D(1, 0);
const CVec2D CVec2D::unitY = CVec2D(0, 1);
const CVec2D CVec2D::nan = CVec2D(CMath::NAN_FLOAT, CMath::NAN_FLOAT);
const CVec2D CVec2D::inf = CVec2D(CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT);

const CVec4D CVec4D::origin( 0.0f, 0.0f, 0.0f, 0.0f );
const CVec4D CVec4D::zero = CVec4D(0, 0, 0, 0);
const CVec4D CVec4D::one = CVec4D(1, 1, 1, 1);
const CVec4D CVec4D::unitX = CVec4D(1, 0, 0, 0);
const CVec4D CVec4D::unitY = CVec4D(0, 1, 0, 0);
const CVec4D CVec4D::unitZ = CVec4D(0, 0, 1, 0);
const CVec4D CVec4D::unitW = CVec4D(0, 0, 0, 1);
const CVec4D CVec4D::nan = CVec4D(CMath::NAN_FLOAT, CMath::NAN_FLOAT, CMath::NAN_FLOAT, CMath::NAN_FLOAT);
const CVec4D CVec4D::inf = CVec4D(CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT);

CVec5D vec5_origin( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );
CVec6D vec6_origin( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );
CVec6D vec6_infinity( CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT );

const CVec3D CVec3D::zero(0.0f,0.0f, 0.0f);
const CVec3D CVec3D::origin(0.0f,0.0f, 0.0f);
const CVec3D CVec3D::one(1.0f, 1.0f, 1.0f);
const CVec3D CVec3D::unitX(1.0f, 0.0f, 0.0f);
const CVec3D CVec3D::unitY(0.0f, 1.0f, 0.0f);
const CVec3D CVec3D::unitZ(0.0f, 0.0f, 1.0f);
const CVec3D CVec3D::nan(CMath::NAN_FLOAT, CMath::NAN_FLOAT, CMath::NAN_FLOAT);
const CVec3D CVec3D::infinity(CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT, CMath::INFINITY_FLOAT);


//===============================================================
//
//	CVec2D
//
//===============================================================

CVec2D		CVec2D::operator=(const CPoint2D &a){
	x = a.x_;
	y = a.y_;
	return *this;
}

CVec2D::CVec2D(CPoint2D const& p){
x=p.x();
y=p.y(); 
}

CVec2D::CVec2D( CPoint2D const &vec1, CPoint2D const &vec2 ){
  *this = CVec2D(vec1)-CVec2D(vec2);

}
/*
=============
CVec2D::toString
=============
*/
const char *CVec2D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

//: Write "<CVec2D x,y> " to stream
ostream&  CVec2D::operator<<(ostream& s)
{
	 CVec2D const& p=*this;
  return s << "<CVec2D "<< p.x << ',' << p.y <<  "> ";
}



//: Read x y from stream

istream&  CVec2D::operator>>(istream& is)
{
	CVec2D & p=*this;
  return p.read(is);
}



//: Read from stream, possibly with formatting
//  Either just reads two blank-separated numbers,
//  or reads two comma-separated numbers,
//  or reads two numbers in parenthesized form "(123, 321)"

istream& CVec2D::read(istream& is)
{
  if (! is.good()) return is; // (TODO: should throw an exception)
  bool paren = false;
  float tx, ty;
  is >> std::ws; // jump over any leading whitespace
  if (is.eof()) return is; // nothing to be set because of EOF (TODO: should throw an exception)
  if (is.peek() == '(') { is.ignore(); paren=true; }
  is >> std::ws >> tx >> std::ws;
  if (is.eof()) return is;
  if (is.peek() == ',') is.ignore();
  is >> std::ws >> ty >> std::ws;
  if (paren) {
    if (is.eof()) return is;
    if (is.peek() == ')') is.ignore();
    else                  return is; // closing parenthesis is missing (TODO: throw an exception)
  }
  set(tx,ty);
  return is;
}
/*
=============
lerp

Linearly inperpolates one vector to another.
=============
*/
void CVec2D::lerp( const CVec2D &v1, const CVec2D &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + l * ( v2 - v1 );
	}
}

float CVec2D::distanceSq(const CVec2D &rhs) const
{
	float dx = x - rhs.x;
	float dy = y - rhs.y;
	return dx*dx + dy*dy;
}

float CVec2D::distance(const CVec2D &rhs) const
{
	return CMath::sqrt(distanceSq(rhs));
}

void CVec2D::orthogonalize(const CVec2D &a, CVec2D &b)
{
	SMF_ASSERT(!a.isZero());
	b -= (a*b) / a.getLenght() * a;
}


CVec2D CVec2D::orthogonalizePos()const{
	SMF_ASSERT(!isZero());
	CVec2D V;
	if (y==0 ) { //the vector is vertical
		V.set(y,x);
	}else {
		V.set(-y, x);
	}
	V.toNormal(); 
	return V;
}
CVec2D CVec2D::orthogonalizeNeg()const{
	SMF_ASSERT(!isZero());
	CVec2D V;
	if (x==0 ) { //the vector is vertical
		V.set(y,x);
	}else {
		V.set(y,-x);
	}
	V.toNormal(); 
	return V;
}

bool CVec2D::areOrthogonal(const CVec2D &a, const CVec2D &b, float epsilon)
{
	return a.isPerpendicular(b, epsilon);
}


void CVec2D::orthoNormalize(CVec2D &a, CVec2D &b)
{
	SMF_ASSERT(!a.isZero());
	a.toNormal();
	b -= (a*b) * a;
}
bool CVec2D::orientedCCW(const CVec2D &a, const CVec2D &b, const CVec2D &c)
{
	// Compute the determinant
	// | ax ay 1 |
	// | bx by 1 |
	// | cx cy 1 |
	// See Christer Ericson, Real-Time Collision Detection, p.32.
	return (a.x-c.x)*(b.y-c.y) - (a.y-c.y)*(b.x-c.x) >= 0.f;
}

//===============================================================
//
//	CVec3D
//
//===============================================================

/*
=============
CVec3D::toYaw
=============
*/
float CVec3D::toYaw() const {
	float yaw;
	
	if ( ( y == 0.0f ) && ( x == 0.0f ) ) {
		yaw = 0.0f;
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}
	}

	return yaw;
}

/*
=============
CVec3D::toPitch
=============
*/
float CVec3D::toPitch() const {
	float	forward;
	float	pitch;
	
	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		forward = ( float )CMath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}

	return pitch;
}

/*
=============
CVec3D::toAngles
=============
*/
CEulerAngles CVec3D::toAngles() const {
	float forward;
	float yaw;
	float pitch;
	
	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		yaw = 0.0f;
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}

		forward = ( float )CMath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}

	return CEulerAngles( -pitch, yaw, 0.0f );
}

/*
=============
CVec3D::toPolar
=============
*/
CPolar3D CVec3D::toPolar() const {
	float forward;
	float yaw;
	float pitch;
	
	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		yaw = 0.0f;
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}

		forward = ( float )CMath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}
	return CPolar3D( CMath::sqrt( x * x + y * y + z * z ), yaw, -pitch );
}

/*
=============
CVec3D::toMat3
=============
*/
CMat3D CVec3D::toMat3() const {
	CMat3D	mat;
	float	d;

	mat[0] = *this;
	d = x * x + y * y;
	if ( !d ) {
		mat[1][0] = 1.0f;
		mat[1][1] = 0.0f;
		mat[1][2] = 0.0f;
	} else {
		d = CMath::invSqrt( d );
		mat[1][0] = -y * d;
		mat[1][1] = x * d;
		mat[1][2] = 0.0f;
	}
	mat[2] = cross( mat[1] );

	return mat;
}

/*
=============
CVec3D::toString
=============
*/
const char *CVec3D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
=============
lerp

Linearly inperpolates one vector to another.
=============
*/
void CVec3D::lerp( const CVec3D &v1, const CVec3D &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + l * ( v2 - v1 );
	}
}

/*
=============
sLerp

Spherical linear interpolation from v1 to v2.
Vectors are expected to be normalized.
=============
*/
#define LERP_DELTA 1e-6

void CVec3D::sLerp( const CVec3D &v1, const CVec3D &v2, const float t ) {
	float omega, cosom, sinom, scale0, scale1;

	if ( t <= 0.0f ) {
		(*this) = v1;
		return;
	} else if ( t >= 1.0f ) {
		(*this) = v2;
		return;
	}

	cosom = v1 * v2;
	if ( ( 1.0f - cosom ) > LERP_DELTA ) {
		omega = CMath::acos( cosom );
		sinom = CMath::sin( omega );
		scale0 = CMath::sin( ( 1.0f - t ) * omega ) / sinom;
		scale1 = CMath::sin( t * omega ) / sinom;
	} else {
		scale0 = 1.0f - t;
		scale1 = t;
	}

	(*this) = ( v1 * scale0 + v2 * scale1 );
}

/*
=============
projectSelfOntoSphere

Projects the z component onto a sphere.
=============
*/
void CVec3D::projectSelfOntoSphere( const float radius ) {
	float rsqr = radius * radius;
	float len = getLenght();
	if ( len  < rsqr * 0.5f ) {
		z = sqrt( rsqr - len );
	} else {
		z = rsqr / ( 2.0f * sqrt( len ) );
	}
}

bool CVec3D::isFinite() const
{
	return MATH::isFinite(x) && MATH::isFinite(y) && MATH::isFinite(z);
}
bool CVec3D::isNormalized(float epsilonSq) const
{
	return fabs(getLengthSqr()-1.f) <= epsilonSq;
}

bool CVec3D::isZero(float epsilonSq) const
{
	return fabs(getLengthSqr()) <= epsilonSq;
}

//==============distance=========================



float CVec3D::distance(const CVec3D &este) const
{
	return CMath::sqrt(distanceSq(este));
}

float CVec3D::distance(const CLine &este) const { return este.distance(*this); }
float CVec3D::distance(const CRay &este) const { return este.distance(*this); }
float CVec3D::distance(const CLineSegment &este) const { return este.distance(*this); }
//float CVec3D::distance(const CPlane &este) const { return este.distance(*this); }
float CVec3D::distance(const CTriangle &este) const { return este.distance(*this); }
float CVec3D::distance(const CAABBox &este) const { return este.distance(*this); }
float CVec3D::distance(const COBBox &este) const { return este.distance(*this); }
float CVec3D::distance(const CSphere &este) const { return este.distance(*this); }

float CVec3D::distanceSq(const CVec3D &rhs) const
{
	float dx = x - rhs.x;
	float dy = y - rhs.y;
	float dz = z - rhs.z;
	return dx*dx + dy*dy + dz*dz;
}

CVec3D CVec3D::perpendicular(const CVec3D &hint, const CVec3D &hint2) const
{
	//s SMF_ASSERT(!this->isZero());
	//s SMF_ASSERT(hint.isNormalized());
	//s SMF_ASSERT(hint2.isNormalized());
	CVec3D v = this->cross(hint);
	float len = v.getLenght();
	v.toNormal();
	if (len == 0)
		return hint2;
	else
		return v;
}

CVec3D CVec3D::anotherPerpendicular(const CVec3D &hint, const CVec3D &hint2) const
{
	//s SMF_ASSERT(!this->isZero());
	//s SMF_ASSERT(hint.isNormalized());
	//s SMF_ASSERT(hint2.isNormalized());
	CVec3D firstPerpendicular = perpendicular(hint, hint2);
	CVec3D v = this->cross(firstPerpendicular);
	return v.normalized();
}

float CVec3D::scaleToLength(float newLength)
{
	float length = getLengthSqr();
	if (length < 1e-6f)
	{
		set(newLength, 0, 0); // Will always produce a vector of the requested length.
		return 0.f;
	}

	length = CMath::sqrt(length);
	float scalar = newLength / length;
	x *= scalar;
	y *= scalar;
	z *= scalar;
	return length;
}

CVec3D CVec3D::scaledToLength(float newLength) const
{
	SMF_ASSERT(!isZero());

	CVec3D v = *this;
	v.scaleToLength(newLength);
	return v;
}

CVec3D CVec3D::projectTo(const CVec3D &direction) const
{
	SMF_ASSERT(!direction.isZero());
	return direction * (*this*direction) / direction.getLengthSqr();
}

CVec3D CVec3D::projectToNorm(const CVec3D &direction) const
{
	SMF_ASSERT(direction.isNormalized());
	return direction * (*this*direction);
}





void CVec3D::orthoNormalize(CVec3D &a, CVec3D &b)
{
	SMF_ASSERT(!a.isZero());
	SMF_ASSERT(!b.isZero());
	a.toNormal();
	b -= b.projectToNorm(a);
	b.toNormal();
}

void CVec3D::orthoNormalize(CVec3D &a, CVec3D &b, CVec3D &c)
{
	SMF_ASSERT(!a.isZero());
	a.toNormal();
	b -= b.projectToNorm(a);
	SMF_ASSERT(!b.isZero());
	b.toNormal();
	c -= c.projectToNorm(a);
	c -= c.projectToNorm(b);
	SMF_ASSERT(!c.isZero());
	c.toNormal();
}

bool MUST_USE_RESULT CVec3D::areOrthonormal(const CVec3D &a, const CVec3D &b, float epsilon)
{
	return a.isPerpendicular(b, epsilon) && a.isNormalized(epsilon*epsilon) && b.isNormalized(epsilon*epsilon);
}

bool MUST_USE_RESULT CVec3D::areOrthonormal(const CVec3D &a, const CVec3D &b, const CVec3D &c, float epsilon)
{
	return a.isPerpendicular(b, epsilon) &&
	       a.isPerpendicular(c, epsilon) &&
	       b.isPerpendicular(c, epsilon) &&
	       a.isNormalized(epsilon*epsilon) &&
	       b.isNormalized(epsilon*epsilon) &&
	       c.isNormalized(epsilon*epsilon);
}

//===============================================================
//
//	CVec4D
//
//===============================================================
#if 0
CVec4D operator*( const float scalar, const CVec4D b ) {
	CVec4D ret;
	if(CMath::MATH_AUTOMATIC_SSE){
	CSIMD::getProcessor()->vector4D_ScaleOf(&ret, &b,scalar);
	}else{
	CSIMD::getGenProcessor()->vector4D_ScaleOf(&ret, &b, scalar);
	//return CVec4D( b.x * scalar, b.y * scalar, b.z * scalar, b.w * scalar );
	}
	return ret;
}
#endif
CVec2D CVec4D::xy() const
{
	return CVec2D(x, y);
}

CVec3D CVec4D::xyz() const
{
	return CVec3D(x, y, z);
}

float CVec4D::getLengthSqr3() const
{
	return x*x + y*y + z*z;

}

float CVec4D::getLength3() const
{
	return CMath::sqrt(x*x + y*y + z*z);

}

const char *CVec4D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
=============
lerp

Linearly inperpolates one vector to another.
=============
*/
void CVec4D::lerp( const CVec4D &v1, const CVec4D &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + (( v2 - v1 )*l);
	}
}

bool CVec4D::isWZeroOrOne(float epsilon) const
{
	return CMath::equalsAbs(w, 0.f, epsilon) || CMath::equalsAbs(w, 1.f, epsilon);
}



//===============================================================
//
//	CVec5D
//
//===============================================================

/*
=============
CVec5D::toString
=============
*/
const char *CVec5D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
=============
CVec5D::lerp
=============
*/
void CVec5D::lerp( const CVec5D &v1, const CVec5D &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		x = v1.x + l * ( v2.x - v1.x );
		y = v1.y + l * ( v2.y - v1.y );
		z = v1.z + l * ( v2.z - v1.z );
		s = v1.s + l * ( v2.s - v1.s );
		t = v1.t + l * ( v2.t - v1.t );
	}
}


//===============================================================
//
//	CVec6D
//
//===============================================================

/*
=============
CVec6D::toString
=============
*/
const char *CVec6D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}


//===============================================================
//
//	CVecXD
//
//===============================================================

float	CVecXD::temp[VECX_MAX_TEMP+4];
float *	CVecXD::tempPtr = (float *) ( ( (uintptr_t) CVecXD::temp + 15 ) & ~15 );
int		CVecXD::tempIndex = 0;

/*
=============
CVecXD::toString
=============
*/
const char *CVecXD::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

} //end MATH
} //end SMF
