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


#ifndef _SMF__SMF_MATH_VECTOR_H__
#define _SMF__SMF_MATH_VECTOR_H__

#include "../SMF_Config.h"
#include "SMF_MathDefs.h"
#include "SMF_Math.h"
#include "../util/SMF_Heap.h"
#include "../util/SMF_Random.h"
#include "../exceptions/all.h"
#include "SMF_Simd.h"

namespace SMF{
namespace GEO{
class CLine;
class CRay;
class CLineSegment;
class CTriangle;
class CAABBox;
class COBBox;
class CSphere;
class CPoint2D;
}
namespace MATH{
class CSIMDProcessor;

/*
===============================================================================

  Vector classes

===============================================================================
*/

#define VECTOR_EPSILON		0.001f

class CEulerAngles;
class CPolar3D;
class CMat3D;
class CPolar3D;

/*
A vector space V is a set that is closed under finite vector addition and scalar multiplication. The basic example is n-dimensional Euclidean space R^n, where every element is represented by a list of n real numbers, scalars are real numbers, addition is componentwise, and scalar multiplication is multiplication on each term separately.

For a general vector space, the scalars are members of a field F, in which case V is called a vector space over F.

Euclidean n-space R^n is called a real vector space, and C^n is called a complex vector space.

In order for V to be a vector space, the following conditions must hold for all elements X,Y,Z in V and any scalars r,s in F:

http://www.euclideanspace.com/maths/algebra/vectors/index.htm
http://mathworld.wolfram.com/VectorSpace.html
http://mathworld.wolfram.com/EuclideanSpace.html
*/


/**
 * \class CVec2D
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de 2  dimensões
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 * \note
 A vector space V is a set that is closed under finite vector addition and scalar multiplication. The basic example is n-dimensional Euclidean space R^n, where every element is represented by a list of n real numbers, scalars are real numbers, addition is componentwise, and scalar multiplication is multiplication on each term separately.

For a general vector space, the scalars are members of a field F, in which case V is called a vector space over F.

Euclidean n-space R^n is called a real vector space, and C^n is called a complex vector space.

In order for V to be a vector space, the following conditions must hold for all elements X,Y,Z in V and any scalars r,s in F:

http://www.euclideanspace.com/maths/algebra/vectors/index.htm
http://mathworld.wolfram.com/VectorSpace.html
http://mathworld.wolfram.com/EuclideanSpace.html

  */
class SMF_API CVec2D {
public:
	float			x;
	float			y;
					/**
					\if pt_br
					\brief inicializa o Vetor n origem (0,0)
					\elseif us_en
					\brief vec(0,0)
					\endif
					**/
					CVec2D();
					explicit CVec2D( const float x, const float y );
					//contructs a vector from other two
					//
					explicit CVec2D( CVec2D const &vec1, CVec2D const &vec2 );
					CVec2D( CPoint2D const &vec1, CPoint2D const &vec2 );

					//: Copy constructor fro CVec2D
					explicit CVec2D(CPoint2D const& p);


		/**
	\brief seta os pontos do vetor
	\param x ponto 1
	\param y ponto 2
	**/
	void 			set( const float x, const float y );
	/**
	\brief zera o vetor
	**/
	void			toZero();
	//: Write "<CVec2D x,y> " to stream

	ostream& operator<<(ostream& s);

	//: Read from stream, possibly with formatting
	//  Either just reads two blank-separated numbers,
	//  or reads two comma-separated numbers,
	//  or reads two numbers in parenthesized form "(123, 321)"

	istream& operator>>(istream& s);

	float			operator[]( int index ) const;
	float &			operator[]( int index );

	/**
	\brief constrói um vetor igual ao passado como parâmetro
	\param a vetor a copiar
	\return vetor igual
	**/
	CVec2D		operator=(const CVec2D &a);
	/**
	\brief constrói um vetor Posicional igual ao ponto passado como parâmetro
	\param a ponto a ser transformado em vetor posicional
	\return vetor igual
	**/
	CVec2D		operator=(const CPoint2D &a);

	CVec2D		operator-() const;
	/*
	\brief produto escalar (ou produto interno) de dois vetores, denotado a . b
	\note  a · b = |a| × |b| × cos(θ)  ou a · b = ax * bx + ay * by
	\note http://www.mathsisfun.com/algebra/vectors-dot-product.html
	\note 0 = the vectors are orthogonal
	\return retorna um valor que é o produto escalar do vetor passado pelo vetor interno
	\param a vetor a ser multiplicado por este
	*/
	float			operator*( const CVec2D &a ) const;
	float			scalar(CVec2D &b)const { return (*this)*b;}
	float			inner(CVec2D &b)const { return (*this)*b;}
	float			dot(CVec2D &b)const { return (*this)*b;}
	/*
	\brief produto vetor por uma constante
	\return retorna um vetor que é o produto deste vetor pela constante passada
	\param a constante a ser multiplicada por este vetor
	*/
	CVec2D		operator*( const float a ) const;
	/*
	\brief divisão do vetor por uma constante
	\return retorna um vetor que é a divisão deste vetor pela constante passada
	\param a constante a ser dividida por este vetor
	*/
	CVec2D		operator/( const float a ) const;
	/**
	\brief divisão de dois vetores
	\return retorna  o vetor resultante da divisão do vetor passado com o vetor deste objeto
	\param a vetor a ser dividido deste
	**/

	CVec2D &		operator/( const CVec2D &a );
    /**
	\brief soma de dois vetores. a + b
	\return retorna um novo vetor que é a soma do vetor passado com o vetor deste objeto
	\param b vetor a ser somado com este
	**/
	CVec2D		operator+( const CVec2D &b ) const;
	/**
	\brief subtração de dois vetores. a - b
	\return retorna um novo vetor que é a subtração do vetor passado com o vetor deste objeto
	\param b vetor a ser subtraído deste este
	**/
	CVec2D		operator-( const CVec2D &b ) const;
	/**
	\brief soma de dois vetores
	\return retorna uma referência para este vetor, que se torna a soma do vetor passado com o vetor deste objeto
	\param a vetor a ser somado com este
	**/
	CVec2D &		operator+=( const CVec2D &a );
	/**
	\brief subtração de dois vetores
	\return retorna uma referência para este vetor, que se torna a subtração do vetor passado com o vetor deste objeto
	\param a vetor a ser subtraído deste
	**/
	CVec2D &		operator-=( const CVec2D &a );

	/**
	\brief divisão de dois vetores
	\return retorna uma referência para este vetor, que se torna a divisão do vetor passado com o vetor deste objeto
	\param a vetor a ser dividido deste
	**/
	CVec2D &		operator/=( const CVec2D &a );
		/*
	\brief divisão deste vetor por uma constante
	\return retorna este vetor, que se torne divisão dele pela constante passada
	\param a constante a ser dividida por este vetor
	*/
	CVec2D &		operator/=( const float a );
	/*
	\brief produto deste vetor por uma constante
	\return retorna este vetor, que se torne o produto dele pela constante passada
	\param a constante a ser multiplicada por este vetor
	*/
	CVec2D &		operator*=( const float a );
	/**
	\brief produto vetorial por constante
	\return retorna um novo vetor, que é o produto do vetor passado pela constante passada
	\param a constante a multiplicar pelo vetorpassado
	\param b vetor a ser multiplicado pelo vetor passado
	**/
	friend CVec2D	operator*( const float a, const CVec2D b );

	bool			compare( const CVec2D &a ) const;							// exact compare, no epsilon
	bool			compare( const CVec2D &a, const float epsilon ) const;		// compare with epsilon
	bool			operator==(	const CVec2D &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const CVec2D &a ) const;						// exact compare, no epsilon

	/// Multiplies this vector by a vector, element-wise.
	/// \note Mathematically, the multiplication of two vectors is not defined in linear space structures,
	///	 but this function is provided here for syntactical convenience.
	/// \return (x*v.x, y*v.y).
	CVec2D mul(const CVec2D &v) const;
	/**
	\brief retorna a norma, ou tamanho do vetor
	\note representada por |v|
	\note http://mathworld.wolfram.com/Norm.html
	**/
	float			getLenght() const;
	/**
	\brief retorna a norma, ou tamanho do vetor
	\note representada por |v|
	\note http://mathworld.wolfram.com/Norm.html
	**/
	float			getLengthFast() const;
	/**
	\brief retorna x^2 + y^2
	**/
	float			getLengthSqr() const;
	/**
	\brief normaliza este vetor. transformando o no vetor unitário  ou versor
	\return retorna o comprimento (length) do vetor
	\note  V^ = V / |v|
	**/
	void			toNormal();
	/**
	\brief normaliza o vetor. calcula o vetor unitário ou versor
	\return retorna o comprimento (length)
	\note  V^ = V / |v|
	**/
	float			normalizeFast();
	/**
	\brief retorna o angulo, em graus entre o vetor passado e este vetor
	\return Ângulo entre os dois vetores, em graus
	\param a vetor ao qual se deseja saber o ângulo
	**/
	float	   getAngle(const CVec2D & a)const;
	/// Returns the angle between this vector and the specified vector, in radians.
	/** \note This function takes into account that this vector or the other vector can be unnormalized, and normalizes the computations.
			If you are computing the angle between two normalized vectors, it is better to use angleBetweenNorm().
		\see angleBetweenNorm(). */
	float angleBetween(const CVec2D &other) const;

	/// Returns the angle between this vector and the specified normalized vector, in radians.
	/** \param normalizedVector The direction vector to compute the angle against. This vector must be normalized.
		\note This vector must be normalized to call this function.
		\see angleBetween(). */
	float angleBetweenNorm(const CVec2D &normalizedVector) const;

	//: Return a normalised version of this vector.
    //  If a is zero length, return (0,0).
	CVec2D	normalized();
	/// \return Min(x, y).
	/** \see minElementIndex(). */
	float minElement() const;
	/// Returns the index that has the smallest value in this vector.
	/** \see minElement(). */
	int minElementIndex() const;
	/// \return Max(x, y).
	/** \see maxElementIndex(). */
	float maxElement() const;
	/// Returns the index that has the smallest value in this vector.
	/** \see maxElement(). */
	int maxElementIndex() const;

	/**
	\brief verifica se dois vetores são perpendiculares
	\return true se os vetores são perpendiculares e falso se não são perpendiculares
	\param a vetor a verificar se é perpendicular a este
	**/
	/// Tests if two vectors are perpendicular to each other.
	/** \see (), isZero(), isPerpendicular(), compare(). */
	bool isPerpendicular(const CVec2D &other, float epsilon =CMath::EPSILON_SuperLow) const;


		/// Tests if the length of this vector is one, up to the given epsilon.
	/** \see isZero(), isFinite(), isPerpendicular(). */
	bool isNormalized(float epsilonSq = 1e-6f) const;

	/// Tests if this is the null vector, up to the given epsilon.
	/** \see (), isFinite(), isPerpendicular(). */
	bool isZero(float epsilonSq = 1e-6f) const;

	/// Tests if this vector contains valid finite elements.
	/** \see (), isZero(), isPerpendicular(). */
	bool isFinite() const;

	/// are two vectors parallel, i.e., is one a scalar multiple of the other?
	// If the third argument is specified, it is taken as the "tolerance", i.e.
	// in that case this function returns true if the vectors are almost parallel.
	/*
	magnitude of u: um = u.Lenght()
    magnitude of v: vm = v.Lenght()

	If two vectors are parallel, then the dot product is equal to the product of the magnitudes:
	ux*vx + uy*vy = um * uv
**/
   bool isParallel(CVec2D const& b, double eps=0.0)const;



	CVec2D &		truncate( float length );	// cap length
	void			clamp( const CVec2D &min, const CVec2D &max );
	void			snap();				// snap to closest integer value
	void			snapInt();			// snap towards integer (floor)
	/**
	\brief retorna a dimensão do vetor
	**/
	int				getDimension() const;
	/**
	\brief retorna um ponteiro para o valor x deste vetor
	**/
	const float *	toFloatPtr() const;
	/**
	\brief retorna um ponteiro para o valor x deste vetor
	**/
	float *			toFloatPtr();
	/**
	\brief retorna uma string C que descreve o vetor
	**/
	const char *	toString( int precision = 2 ) const;
	//Linearly inperpolates one vector to another.
	void			lerp( const CVec2D &v1, const CVec2D &v2, const float l );
	/**
	\brief cross product of two 2D vectors
	\note The magnitude of the cross product can be interpreted as the positive area of enclosed parallellogram.
	\note \f$ \vec{V} x \vec{U} = \| \vec{V} \|. \| \vec{U} \| .sin (\theta) =   V_x.U_y - V_y.U_x   \f$
	\note Properties:

	\f$ If \quad  \vec{V} \f$  and  \f$ \vec{U} \f$  are collinear, \f$ \vec{V} × \vec{U}\f$  is a null vector.

	\f$ If \quad \vec{V} \f$ is null, \f$ \vec{V} × \vec{U} \f$  is a null vector.

	\f$ If \quad \vec{U} \f$ is null,  \f$ \vec{V} × \vec{U} \f$  is a null vector.

	\f$\vec{V} × \vec{U} \f$ is normal to the plane containing the vectors \f$\vec{V} \f$  and  \f$\vec{U}\f$ .

	\f$\vec{V} × \vec{U} =−(\vec{U} × \vec{V} )=(−\vec{U} )×\vec{V}  \f$

	\see http://en.wikipedia.org/wiki/Cross_product#Geometric_meaning
	**/
	float cross(CVec2D const& b) { return ((x *b.y) - (y*b.x)); }


	/// Makes the given vectors linearly independent.
	/** This function directly follows the Gram-Schmidt procedure on the input vectors.
		The vector a is kept unmodified, and vector b is modified to be perpendicular to a.
		\note If any of the input vectors is zero, then the resulting set of vectors cannot be made orthogonal.
		\see areOrthogonal(), orthoNormalize(), areOrthonormal(). */
	static void orthogonalize(const CVec2D &a, CVec2D &b);
	/**
	\return a vector normal to this vector on the anti-clockwise direction
	\note this vector is considered a position Vector, ie origin (0,0)
	Firstly, one must know that: if two vectors are perpendicular, their dot product equals zero.
	The normal vector (x',y') is perpendicular to the line connecting (x1,y1) and (x2,y2). This line has direction (x2-x1,y2-y1), or (dx,dy).
	So,

	(x',y').(dx,dy) = 0
	x'.dx + y'.dy = 0

	The are plenty of pairs (x',y') that satisfy the above equation. But the best pair that ALWAYS satisfies is either (dy,-dx) or (-dy,dx)

	**/
	CVec2D orthogonalizePos()const;
	/**
	\return a vector normal to this vector on the clockwise direction
	**/
	CVec2D orthogonalizeNeg()const;

	/// Returns true if the given vectors are orthogonal to each other.
	// are two vectors orthogonal, i.e., is their dot product zero?
    // If the third argument is specified, it is taken as the "tolerance", i.e.
    // in that case this function returns true if the vectors are almost orthogonal.

	/** \see orthogonalize(), orthoNormalize(), areOrthonormal(). */
	static bool areOrthogonal(const CVec2D &a, const CVec2D &b, float epsilon =CMath::EPSILON_SuperLow);

	/// Makes the given vectors linearly independent and normalized in length.
	/** This function directly follows the Gram-Schmidt procedure on the input vectors.
		The vector a is first normalized, and vector b is modified to be perpendicular to a, and also normalized.
		\note If either of the input vectors is zero, then the resulting set of vectors cannot be made orthonormal.
		\see orthogonalize(), areOrthogonal(), areOrthonormal(). */
	static void orthoNormalize(CVec2D &a, CVec2D &b);

	/// Tests if the triangle a->b->c is oriented counter-clockwise.
	/** Returns true if the triangle a->b->c is oriented counter-clockwise, when viewed in the XY-plane
		where x spans to the right and y spans up.
		Another way to think of this is that this function returns true, if the point C lies to the left
		of the directed line AB. */
	static bool orientedCCW(const CVec2D &a, const CVec2D &b, const CVec2D &c);
		/// Specifies a compile-time constant CVec2D with value (0, 0).
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D origin;
	/// Specifies a compile-time constant CVec2D with value (1, 1). [similarOverload: zero]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D one;
	/// Specifies a compile-time constant CVec2D with value (1, 0).
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D unitX;
	/// Specifies a compile-time constant CVec2D with value (0, 1). [similarOverload: unitX]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D unitY;
	/// A compile-time constant CVec2D with value (NaN, NaN).
	/** For this constant, each element has the value of quiet NaN, or Not-A-Number.
		\note Never compare a CVec2D to this value! Due to how IEEE floats work, for each float x, both expressions "x == nan" and "x != nan" return false!
			  That is, nothing is equal to NaN, not even NaN itself!
		\note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D nan;
	/// A compile-time constant CVec2D with value (+infinity, +infinity). [similarOverload: nan]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec2D inf;


		/// Computes the distance between this and the given CVec2D.
	/** \see distanceSq(), getLenght(), getLengthSqr(). */
	float distance(const CVec2D &point) const;

	/// Computes the squared distance between this and the given point.
	/** Calling this function is faster than calling distance(), since this function avoids computing a square root.
		If you only need to compare distances to each other, but are not interested in the actual distance values,
		you can compare by using distanceSq(), instead of distance(), since sqrt() is an order-preserving
		(monotonous and non-decreasing) function.
		\see distance(), getLenght(), getLengthSqr(). */
	float distanceSq(const CVec2D &point) const;

	 //: Read from stream, possibly with formatting
  //  Either just reads two blank-separated numbers,
  //  or reads two comma-separated numbers,
  //  or reads two numbers in parenthesized form "(123, 321)"
  istream& read(istream& is);
	//  +-+-+ CVec2D simple I/O +-+-+

	static void test_vector_2d();

public:


};




/**
http://mathworld.wolfram.com/ZeroVector.html
\note A zero vector, denoted 0, is a vector of length 0, and thus has all components equal to zero.
      It is the additive identity of the additive group of vectors.
**/

SMF_INLINE_FORCED CVec2D::CVec2D() {
	x=0;
	y=0;
}



SMF_INLINE_FORCED CVec2D::CVec2D( const float x, const float y ) {
	this->x = x;
	this->y = y;
}

SMF_INLINE_FORCED CVec2D::CVec2D( CVec2D const &vec1, CVec2D const &vec2 ){
  *this = vec1-vec2;

}


SMF_INLINE_FORCED void CVec2D::set( const float x, const float y ) {
	this->x = x;
	this->y = y;
}

SMF_INLINE_FORCED void CVec2D::toZero() {
	x = y = 0.0f;
}

SMF_INLINE_FORCED bool CVec2D::compare( const CVec2D &a ) const {
	return ( ( x == a.x ) && ( y == a.y ) );
}

SMF_INLINE_FORCED bool CVec2D::compare( const CVec2D &a, const float epsilon ) const {
	if ( CMath::fabs( x - a.x ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( y - a.y ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED CVec2D		CVec2D::operator=(const CVec2D &a){
	x = a.x;
	y = a.y;
	return *this;
}


SMF_INLINE_FORCED bool CVec2D::operator==( const CVec2D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CVec2D::operator!=( const CVec2D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED float CVec2D::operator[]( int index ) const {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float& CVec2D::operator[]( int index ) {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED CVec2D CVec2D::mul(const CVec2D &rhs) const
{
	return CVec2D(x * rhs.x, y * rhs.y);
}

SMF_INLINE_FORCED float CVec2D::getLenght() const {
	float result  = ( float )CMath::sqrt( x * x + y * y );
	return result;
}

SMF_INLINE_FORCED float CVec2D::getLengthFast() const {
	float sqrLength;

	sqrLength = x * x + y * y;
	// reciprocal square root, returns huge number when x == 0.0
	return sqrLength * CMath::rSqrt( sqrLength );
}

SMF_INLINE_FORCED float CVec2D::getLengthSqr() const {
	return ( x * x + y * y );
}

SMF_INLINE_FORCED void CVec2D::toNormal() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y;
	// inverse square root with 32 bits precision, returns huge number when x == 0.0
	invLength = CMath::invSqrt( sqrLength );
	x *= invLength;
	y *= invLength;

}
SMF_INLINE_FORCED CVec2D CVec2D::normalized() {
	CVec2D vec=*this;
	 vec.toNormal();
	 return vec;
}
SMF_INLINE_FORCED float CVec2D::normalizeFast() {
	float lengthSqr, invLength;

	lengthSqr = x * x + y * y;
	invLength = CMath::rSqrt( lengthSqr );
	x *= invLength;
	y *= invLength;
	return invLength * lengthSqr;
}

SMF_INLINE_FORCED float CVec2D::minElement() const
{
	return MIN(x, y);
}

SMF_INLINE_FORCED int CVec2D::minElementIndex() const
{
	return (x <= y) ? 0 : 1;
}

SMF_INLINE_FORCED float CVec2D::maxElement() const
{
	return  MAX(x, y);
}

SMF_INLINE_FORCED int CVec2D::maxElementIndex() const
{
	return (x > y) ? 0 : 1;
}
SMF_INLINE_FORCED float		CVec2D::getAngle(const CVec2D & b)const{

CVec2D const& a=*this;
#if 0
float ang1 = CMath::acos(x/a.getLenght());
float ang2 = CMath::acos(b.x/b.getLenght());
float ang3 = CMath::fabs(ang1-ang2);
//return ang3;
#endif
  float ab = (a*b);
  float a_b = (  CMath::sqrt( CMath::fabs(a.getLengthSqr() * b.getLengthSqr()) ));

  float param=  ab / a_b;
  float ang4= CMath::acos(param);
  return ang4;

#if 0

	CVec2D este = *this;
	float dot = este * a;
	int lent = CMath::abs(getLengthFast() * a.getLengthFast());
	float param = dot / lent;
	float result = CMath::acos(param) * 180.0 / CMath::PI;
	return result;
#endif
}
SMF_INLINE_FORCED float CVec2D::angleBetween(const CVec2D &other) const
{
	return CMath::acos((*this)*(other) / CMath::sqrt(getLengthSqr() * other.getLengthSqr()));
}

SMF_INLINE_FORCED float CVec2D::angleBetweenNorm(const CVec2D &other) const
{
	SMF_ASSERT(this->isNormalized());
	SMF_ASSERT(other.isNormalized());
	return CMath::acos((*this)*(other));
}
SMF_INLINE_FORCED bool CVec2D::isPerpendicular(const CVec2D &other, float epsilon) const
{
	return fabs((*this)*(other)) <= epsilon;
}
SMF_INLINE_FORCED bool CVec2D::isParallel(CVec2D const& b, double eps)const
{
  float cross = (x *b.y) - (y*b.x);//this->cross(b); // should be zero
  if (eps <= 0 || cross == float(0)) return cross == float(0);
  eps *= eps * getLenght() * b.getLenght();
  return cross*cross < eps;

#if 0
	CVec2D cross = (*this).mul(b);
  float soma = cross.x+cross.y;
  float soma2 =(this->getLenght() * b.getLenght());
  //if (eps <= 0 || cross == CVec2D(0,0)) return cross == CVec2D(0,0);
  return (soma - soma2)<=eps;
#endif
}
SMF_INLINE_FORCED bool	CVec2D::isNormalized(float epsilonSq) const
{
	return fabs(getLengthSqr()-1.f) <= epsilonSq;
}

SMF_INLINE_FORCED bool	CVec2D::isZero(float epsilonSq) const
{
	return fabs(getLengthSqr()) <= epsilonSq;
}

SMF_INLINE_FORCED bool	CVec2D::isFinite() const
{
	return MATH::isFinite(x) && MATH::isFinite(y);
}
SMF_INLINE_FORCED CVec2D &CVec2D::truncate( float length ) {
	float length2;
	float ilength;

	if ( !length ) {
		toZero();
	}
	else {
		length2 = getLengthSqr();
		if ( length2 > length * length ) {
			ilength = length * CMath::invSqrt( length2 );
			x *= ilength;
			y *= ilength;
		}
	}

	return *this;
}

SMF_INLINE_FORCED void CVec2D::clamp( const CVec2D &min, const CVec2D &max ) {
	if ( x < min.x ) {
		x = min.x;
	} else if ( x > max.x ) {
		x = max.x;
	}
	if ( y < min.y ) {
		y = min.y;
	} else if ( y > max.y ) {
		y = max.y;
	}
}

SMF_INLINE_FORCED void CVec2D::snap() {
	//Rounds x downward, returning the largest integral value that is not greater than x.
	x = CMath::floor( x + 0.5f );
	y = CMath::floor( y + 0.5f );
}

SMF_INLINE_FORCED void CVec2D::snapInt() {
	x = float( int( x ) );
	y = float( int( y ) );
}

SMF_INLINE_FORCED CVec2D CVec2D::operator-() const {
	return CVec2D( -x, -y );
}

SMF_INLINE_FORCED CVec2D CVec2D::operator-( const CVec2D &a ) const {
	float xx=x - a.x;
	float yy=y - a.y;
	return CVec2D(xx , yy );
}

SMF_INLINE_FORCED float CVec2D::operator*( const CVec2D &a ) const {
	return x * a.x + y * a.y;
}

SMF_INLINE_FORCED CVec2D CVec2D::operator*( const float a ) const {
	return CVec2D( x * a, y * a );
}

SMF_INLINE_FORCED CVec2D CVec2D::operator/( const float a ) const {
	float inva = 1.0f / a;
	return CVec2D( x * inva, y * inva );
}

SMF_INLINE_FORCED CVec2D operator*( const float a, const CVec2D b ) {
	return CVec2D( b.x * a, b.y * a );
}

SMF_INLINE_FORCED CVec2D CVec2D::operator+( const CVec2D &a ) const {
	return CVec2D( x + a.x, y + a.y );
}

SMF_INLINE_FORCED CVec2D &CVec2D::operator+=( const CVec2D &a ) {
	x += a.x;
	y += a.y;

	return *this;
}

SMF_INLINE_FORCED CVec2D &CVec2D::operator/=( const CVec2D &a ) {
	x /= a.x;
	y /= a.y;

	return *this;
}
SMF_INLINE_FORCED CVec2D &CVec2D::operator/( const CVec2D &a ) {
	CVec2D vec(x / a.x,
	y / a.y);

	return vec;
}

SMF_INLINE_FORCED CVec2D &CVec2D::operator/=( const float a ) {
	float inva = 1.0f / a;
	x *= inva;
	y *= inva;

	return *this;
}

SMF_INLINE_FORCED CVec2D &CVec2D::operator-=( const CVec2D &a ) {
	x -= a.x;
	y -= a.y;

	return *this;
}

SMF_INLINE_FORCED CVec2D &CVec2D::operator*=( const float a ) {
	x *= a;
	y *= a;

	return *this;
}

SMF_INLINE_FORCED int CVec2D::getDimension() const {
	return 2;
}

SMF_INLINE_FORCED const float *CVec2D::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CVec2D::toFloatPtr() {
	return &x;
}

class CVec4D;
/**
 * \class CVec3D
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de 3  dimensões
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CVec3D {
public:
	float			x;
	float			y;
	float			z;

					CVec3D();
					explicit CVec3D( const float xyz ) { set( xyz, xyz, xyz ); }
					explicit CVec3D( const float x, const float y, const float z );

	void 			set( const float x, const float y, const float z );

	void			toZero();
	/// \return Min(x, y, z).
	/** \see minElementIndex(). */
	float minElement() const;
	/// Returns the index that has the smallest value in this vector.
	/** \see minElement(). */
	int minElementIndex() const;
	/// \return MAX(x, y, z).
	/** \see maxElementIndex(). */
	float maxElement() const;
	/// Returns the index that has the smallest value in this vector.
	/** \see maxElement(). */
	int maxElementIndex() const;
	/// Takes the element-wise absolute value of this vector.
	/** \return CVec3D(|x|, |y|, |z|).
		\see Neg(). */
	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CVec3D		operator-() const;
	CVec3D &		operator=( const CVec3D &a );		// required because of a msvc 6 & 7 bug
	/**
	\brief dot
	/// Computes the dot product of this and the given vector \f$( \vec{U} )\f$.
	    The dot product has a geometric interpretation of measuring how close two direction vectors are to pointing
		in the same direction, computing angles between vectors, or the length of a projection of one vector to another.
		\return \f$  _x.U_x + _y.U_y + _z.U_z  \f$
		\note DOT: \f$ \quad \vec{V} \cdot \vec{U} = V_x.U_x + V_y.U_y + V_z.U_z = \| \vec{V} \|. \| \vec{U} \| .cos (\theta) \f$
		\p where θ is the angle between \f$ \vec{V} \f$  and \f$ \vec{U} \f$ .

Properties:

	 \f$ If \quad \vec{V} \f$  and \f$ \vec{U} \f$  are perpendicular, \f$\vec{U} \cdot \vec{V} \f$  is equal to zero.

	 \f$ If \quad \vec{V} \f$  is null, \f$\vec{U} \cdot \vec{V} \f$  is equal to zero.

	 \f$ If \quad \vec{U} \f$  is null, \f$\vec{U} \cdot \vec{V} \f$  is equal to zero.

	 \f$ \vec{V} \cdot \vec{V} = \| \vec{V} \| ^2 \f$

	 \f$ \vec{V} \cdot \vec{U} = \vec{U} \cdot \vec{V} \f$

	 \f$ \vec{V} \cdot (-\vec{U}) = (-\vec{V}) \cdot \vec{U} = - ( \vec{V} \cdot \vec{U} ) \f$

	\see angleBetween(), projectTo(), projectToNorm(), cross(), OuterProduct(), ScalarTripleProduct().

	\see http://plaza.obu.edu/corneliusk/mp/vrdp.pdf
	**/
	float		operator*( const CVec3D &v ) const;
	float		scalar(const CVec3D &v ) const { return (*this)*v;};
	float		inner(const CVec3D &v ) const  { return (*this)*v;};

	CVec3D		operator*( const float a ) const;
	CVec3D		operator/( const float a ) const;
	CVec3D		operator+( const CVec3D &a ) const;
	CVec3D		operator-( const CVec3D &a ) const;
	CVec3D &		operator+=( const CVec3D &a );
	CVec3D &		operator-=( const CVec3D &a );
	CVec3D &		operator/=( const CVec3D &a );
	CVec3D &		operator/=( const float a );
	CVec3D &		operator*=( const float a );

	friend CVec3D	operator*( const float a, const CVec3D b );
	/// Takes the element-wise absolute value of this vector.
	/** \return CVec3D(|x|, |y|, |z|).
		\see Neg(). */
	CVec3D abs() const;
	bool			compare( const CVec3D &a ) const;							// exact compare, no epsilon
	bool			compare( const CVec3D &a, const float epsilon ) const;		// compare with epsilon
	bool			operator==(	const CVec3D &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const CVec3D &a ) const;						// exact compare, no epsilon
	// http://en.wikipedia.org/wiki/Degenerate_form
	/**
	\brief corrige casos de axiais normais degenerados
	\return verdadeiro se foi feita alguma correção no vetor
	\note casos:
				x=0, y=0, z!=1   (z>0)
				x=0, y=0, z!=-1  (z<0)
				x=0, y!=1, z=0   (y>0)
				x=0, y!=-1, z=0  (y<0)
				x!=1, y=0, z=0   (x>0)
				x!=-1, y=0, z=0  (x<0)

				ou

				Abx(x)=1 e (y != 0 || z !=0)
				abs(y)=1 e (x != 0 || z !=0)
				abs(z)=1 e (x != 0 || y !=0)
	**/
	bool			fixDegenerateNormal();	// fix degenerate axial cases
	/**
	\brief corrige números muito pequenos para zero
	\return verdadeiro se alguma normal foi corrigidas
	**/
	bool			fixDenormals();			// change tiny numbers to zero
	/**
	\if pt_br
	\brief produto vetorial  este cross a
	\f[
	\vec{V} \times \vec{U} = \left (
	\begin{matrix}
	  V_y.U_z - V_z.U_y \\
	  V_z.U_x - V_x.U_z \\
	  V_x.U_y - V_y.U_x
	 \end{matrix}
	\right )
	\f]
	\return um novo vetor que é a multiplicação do vetor passado por este
	\elseif us_en
	\note
	\vec{V} \times \vec{U} = \left (
	\begin{matrix}
	  V_y.U_z - V_z.U_y \\
	  V_z.U_x - V_x.U_z \\
	  V_x.U_y - V_y.U_x
	 \end{matrix}
	\right )
	\note Properties
	If V⃗  and U⃗  are collinear, V⃗ ×U⃗  is a null vector.
	If V⃗  is null, V⃗ ×U⃗  is a null vector.
	If U⃗  is null, V⃗ ×U⃗  is a null vector.
	V⃗ ×U⃗  is normal to the plane containing the vectors V⃗  and U⃗ .
	V⃗ ×U⃗ =−(U⃗ ×V⃗ )=(−U⃗ )×V⃗
	\endif
	**/
	CVec3D		cross( const CVec3D &a ) const;
	/**
	\brief produto vetorial a x b
	\param a primeiro vetor a ser multiplicado
	\param b segundo vetor a ser multiplicado
	\return este vetor, o qual  se torna a multiplicação de a x b
	**/
	CVec3D &	cross( const CVec3D &a, const CVec3D &b );
	/**
	\brief retorna a norma, ou tamanho do vetor
	\note representada por |v|
	\note http://mathworld.wolfram.com/Norm.html
	**/
	float			getLenght() const;
	/**
	\brief retorna a norma, ou tamanho do vetor
	\note representada por |v|
	\note http://mathworld.wolfram.com/Norm.html
	**/
	float			getLengthFast() const;
	/**
	\brief retorna x^2 + y^2 + z^2
	**/
	float			getLengthSqr() const;
	/**
	\brief normaliza o vetor. calculo o vetor unitáriocalcula o vetor unitário ou versor
	\return retorna o comprimento (length)
	\note  V^ = V / |v|
	**/
	void			toNormal();
	SMF_INLINE void normalize(){ toNormal();};
	/**
	\brief normaliza o vetor. calcula o vetor unitário ou versor
	\return retorna o comprimento (length)
	\note  V^ = V / |v|
	**/
	float			normalizeFast();

	/// Returns a normalized copy of this vector.
	/** \note If the vector is zero and cannot be normalized, the vector (1, 0, 0) is returned, and an error message is printed.
			If you do not want to generate an error message on failure, but want to handle the failure yourself, use the
			toNormal() function instead.
		\see toNormal(). */
	CVec3D normalized() const;
	/// Makes the given vectors linearly independent and normalized in length.
	/** This function directly follows the Gram-Schmidt procedure on the input vectors.
		The vector a is first normalized, and vector b is modified to be perpendicular to a, and also normalized.
		Finally, if specified, the vector c is adjusted to be perpendicular to a and b, and normalized.
		\note If any of the input vectors is zero, then the resulting set of vectors cannot be made orthonormal.
		\see orthogonalize(), areOrthogonal(), areOrthonormal(). */
	static void orthoNormalize(CVec3D &a, CVec3D &b);
	static void orthoNormalize(CVec3D &a, CVec3D &b, CVec3D &c);

	/// Returns true if the given vectors are orthogonal to each other and all of length 1.
	/** \see orthogonalize(), areOrthogonal(), orthoNormalize(), areCollinear(). */
	static MUST_USE_RESULT bool areOrthonormal(const CVec3D &a, const CVec3D &b, float epsilon =CMath::EPSILON_SuperLow);
	static MUST_USE_RESULT bool areOrthonormal(const CVec3D &a, const CVec3D &b, const CVec3D &c, float epsilon =CMath::EPSILON_SuperLow);

	/**
	\brief retorna o angulo, em graus entre o vetor passado e este vetor
	\return Ângulo entre os dois vetores, em graus
	\param a vetor ao qual se deseja saber o ângulo
	**/
	float		    getAngle(CVec3D & a);
	// cap length
	CVec3D &		truncate( float length );
	void			clamp( const CVec3D &min, const CVec3D &max );
	CVec3D			clamp( const CVec3D &min, const CVec3D &max ) const;
	/// snap to closest integer value
	void			snap();
	/// snap towards integer (floor)
	void			snapInt();

	int				getDimension() const;

	float			toYaw() const;
	float			toPitch() const;
	CEulerAngles	toAngles() const;
	CPolar3D		toPolar() const;
	/// This vector became column 0 of the generated matriz (CMat3D is column major)
	CMat3D			toMat3() const;		// vector should be normalized
	const CVec2D &	toVec2() const;

	CVec2D &		toVec2();
	CVec2D &		project_xy()const;
	CVec2D &		project_zy()const;
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;
	/// Generates a new CVec3D by filling its entries by the given scalar.
	/** \see CVec3D::CVec3D(float scalar), setFromScalar(). */
	static MUST_USE_RESULT CVec3D fromScalar(float scalar);

	/// Fills each entry of this CVec3D by the given scalar.
	/** \see CVec3D::CVec3D(float scalar), fromScalar(). */
	void setFromScalar(float scalar);

	void			normalVectors( CVec3D &left, CVec3D &down ) const;	// vector should be normalized
	void			orthogonalBasis( CVec3D &left, CVec3D &up ) const;

	void			projectOntoPlane( const CVec3D &normal, const float overBounce = 1.0f );
	bool			projectAlongPlane( const CVec3D &normal, const float epsilon, const float overBounce = 1.0f );
	void			projectSelfOntoSphere( const float radius );

	//Linearly inperpolates one vector to another.
	void			lerp( const CVec3D &v1, const CVec3D &v2, const float l );
	void			sLerp( const CVec3D &v1, const CVec3D &v2, const float l );
	/**
	\brief verifica se dois vetores são perpendiculares
	\return true se os vetores são perpendiculares e falso se não são perpendiculares
	\param a vetor a verificar se é perpendicular a este
	**/
	/// Tests if two vectors are perpendicular to each other.
	/** \see (), isZero(), isPerpendicular(), compare(). */
	bool	isPerpendicular(const CVec3D &other, float epsilon =CMath::EPSILON_SuperLow) const;

	/// Tests if the points p1, p2 and p3 lie on a straight line, up to the given epsilon.
	/** \see areOrthogonal(), areOrthonormal(), Line::areCollinear(). */
	static MUST_USE_RESULT bool areCollinear(const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, float epsilon =CMath::EPSILON_SuperLow);

	/// Computes a new normalized direction vector that is perpendicular to this vector and the specified hint vector.
	/** If this vector points toward the hint vector, the vector hint2 is returned instead.
		\see anotherPerpendicular(), cross(). */
	CVec3D perpendicular(const CVec3D &hint = CVec3D(0,1,0), const CVec3D &hint2 = CVec3D(0,0,1)) const;

	/// Returns another vector that is perpendicular to this vector and the vector returned by perpendicular().
	/** The set (this, perpendicular(), anotherPerpendicular()) forms a right-handed normalized 3D basis.
		\see perpendicular(), cross(). */
	CVec3D anotherPerpendicular(const CVec3D &hint = CVec3D(0,1,0), const CVec3D &hint2 = CVec3D(0,0,1)) const;

	/// Tests if this vector contains valid finite elements.
	/**
	\brief testa se este vetor comtém elementos finitos e válidos
	\see (), isZero(), isPerpendicular(). */
	bool			isFinite() const;

	/// Tests if the length of this vector is one, up to the given epsilon.
	/**
	\brief Testa se o length do vetor é um, até o epsilon dado
	\param epsilonSq (http://pt.wikipedia.org/wiki/%C3%89psilon_de_m%C3%A1quina)
	\see isZero(), isFinite(), isPerpendicular().
	*/
	bool			isNormalized(float epsilonSq = 1e-6f) const;

	/// Tests if this is the null vector, up to the given epsilon.
	/**
	\brief Testa se este vetor é nulo, até o epsilon dado
	\param epsilonSq (http://pt.wikipedia.org/wiki/%C3%89psilon_de_m%C3%A1quina)
	\see (), isFinite(), isPerpendicular(). *
	*/
	bool			isZero(float epsilonSq = 1e-6f) const;

//static members

	/// Specifies a compile-time constant CVec3 with value (0.0f, 0.0f, 0.0f).
	/**
	\brief Especifica uma constante Vector3 com valor (0.0f, 0.0f, 0.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D zero;
	static const CVec3D origin;
	/// Specifies a compile-time constant CVec3D with value (1.0f, 1.0f, 1.0f). [similarOverload: zero]
	/**
	\brief Especifica uma constante Vector3 com valor (1.0f, 1.0f, 1.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D one;
	/// Specifies a compile-time constant CVec3D with value (1.0f, 0.0f, 0.0f).
	/**
	\brief Especifica uma constante Vector3 com valor (1.0f, 0.0f, 0.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D unitX;
	/// Specifies a compile-time constant CVec3D with value (0.0f, 1.0f, 0.0f). [similarOverload: unitX]
	/**
	\brief Especifica uma constante Vector3 com valor (0.0f, 1.0f, 0.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D unitY;
	/// Specifies a compile-time constant CVec3D with value (0.0f, 0.0f, 1.0f). [similarOverload: unitX]
	/**
	/**
	\brief Especifica uma constante Vector3 com valor (0.0f, 0.0f, 1.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D unitZ;
	/// A compile-time constant CVec3D with value (NaN, NaN, NaN).
	/**
	    \brief uma constante CVec3D com valor (NaN, NaN, NaN).
	    For this constant, each element has the value of quiet NaN, or Not-A-Number.
		\warning Nunca compare um CVec3D com este valor! "x == nan" and "x != nan" sempre retornará falso.
		\warning Never compare a CVec3D to this value! Due to how IEEE floats work, for each float x, both expressions "x == nan" and "x != nan" return false!
			  That is, nothing is equal to NaN, not even NaN itself!
		\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
		\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec3D nan;
	/// A compile-time constant CVec3D with value (+infinity, +infinity, +infinity). [similarOverload: nan]
	/**
	\brief Especifica uma constante Vector3 com valor (+infinity, +infinity, +infinity).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units!
	**/
	static const CVec3D infinity;
	//=============MIN e MAX===============

	/// Returns an element-wise minimum of this and the vector (ceil, ceil, ceil).
	/** Each element that is larger than ceil is replaced by ceil. */
	CVec3D Min(float ceil) const;
	/// Returns an element-wise minimum of this and the given vector.
	/** Each element that is larger than ceil is replaced by ceil.
		\see MAX(), clamp(). */
	CVec3D Min(const CVec3D &ceil) const;
	/// Returns an element-wise maximum of this and the vector (floor, floor, floor).
	/** Each element that is smaller than floor is replaced by floor. */
	CVec3D Max(float floor) const;
	/// Returns an element-wise maximum of this and the given vector.
	/** Each element that is smaller than floor is replaced by floor.
		\see Min(), clamp(). */
	CVec3D Max(const CVec3D &floor) const;


	//========distance =============================

	/// Computes the distance between this point and the given object.
	/** This function finds the nearest point to this point on the given object, and computes its distance
		to this point.
		If this point lies inside the given object, a distance of 0 is returned.
		\todo add CVec3D::distance(Polygon/Circle/Disc/Frustum/Polyhedron).
		\see distanceSq(), getLenght(), getLengthSqr(). */
	float distance(const CVec3D &point) const;
	float distance(const CLine &line) const;
	float distance(const CRay &ray) const;
	float distance(const CLineSegment &lineSegment) const;
	//float distance(const CPlane &plane) const;
	float distance(const CTriangle &triangle) const;
	float distance(const CAABBox &aabb) const;
	float distance(const COBBox &obb) const;
	float distance(const CSphere &sphere) const;
		/// Computes the squared distance between this and the given point.
	/** Calling this function is faster than calling distance(), since this function avoids computing a square root.
		If you only need to compare distances to each other, but are not interested in the actual distance values,
		you can compare by using distanceSq(), instead of distance(), since sqrt() is an order-preserving
		(monotonous and non-decreasing) function.
		\see distance(), getLenght(), getLengthSqr(). */
	float distanceSq(const CVec3D &point) const;

	/// Scales this vector so that its new length is as given.
	/** Calling this function is effectively the same as normalizing the vector first and then multiplying by newLength.
		In the case of failure, this vector is set to (newLength, 0, 0), so calling this function will never result in an
		unnormalized vector.
		\note This function operates in-place.
		\return The old length of this vector. If this function returns 0, the scaling failed, and this vector is arbitrarily
			reset to (newLength, 0, 0). In case of failure, no error message is generated. You are expected to handle the failure
			yourself.
		\see scaledToLength(). */
	float scaleToLength(float newLength);

	/// Returns a scaled copy of this vector which has its new length as given.
	/** This function assumes the length of this vector is not zero. In the case of failure, an error message is printed,
		and the vector (newLength, 0, 0) is returned.
		\see scaleToLength(). */
	CVec3D scaledToLength(float newLength) const;

		/// Projects this vector onto the given unnormalized direction vector.
	/** \param direction The direction vector to project this vector onto. This function will normalize this
			vector, so you can pass in an unnormalized vector.
		\see projectToNorm(). */
	CVec3D projectTo(const CVec3D &direction) const;

	/// Projects this vector onto the given normalized direction vector.
	/** \param direction The vector to project onto. This vector must be normalized.
		\see projectTo(). */
	CVec3D projectToNorm(const CVec3D &direction) const;


};




SMF_INLINE_FORCED CVec3D::CVec3D() {
}

SMF_INLINE_FORCED CVec3D::CVec3D( const float x, const float y, const float z ) {
	this->x = x;
	this->y = y;
	this->z = z;
}

SMF_INLINE_FORCED float CVec3D::operator[]( const int index ) const {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float &CVec3D::operator[]( const int index ) {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED void CVec3D::set( const float x, const float y, const float z ) {
	this->x = x;
	this->y = y;
	this->z = z;
}

SMF_INLINE_FORCED void CVec3D::toZero() {
	x = y = z = 0.0f;
}
SMF_INLINE_FORCED CVec3D CVec3D::abs() const
{
	return CVec3D(fabs(x), fabs(y), fabs(z));
}

SMF_INLINE_FORCED float CVec3D::minElement() const
{
	return MIN(MIN(x, y), z);
}

SMF_INLINE_FORCED int CVec3D::minElementIndex() const
{
	if (x <= y && x <= z)
		return 0;
	else
		return (y <= z) ? 1 : 2;
}

SMF_INLINE_FORCED float CVec3D::maxElement() const
{
	return MAX(MAX(x, y), z);
}

SMF_INLINE_FORCED int CVec3D::maxElementIndex() const
{
	if (x >= y && x >= z)
		return 0;
	else
		return (y >= z) ? 1 : 2;
}


SMF_INLINE_FORCED CVec3D CVec3D::operator-() const {
	return CVec3D( -x, -y, -z );
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator=( const CVec3D &a ) {
	x = a.x;
	y = a.y;
	z = a.z;
	return *this;
}

SMF_INLINE_FORCED CVec3D CVec3D::operator-( const CVec3D &a ) const {
	return CVec3D( x - a.x, y - a.y, z - a.z );
}

SMF_INLINE_FORCED float CVec3D::operator*( const CVec3D &a ) const {
	return x * a.x + y * a.y + z * a.z;
}

SMF_INLINE_FORCED CVec3D CVec3D::operator*( const float a ) const {
	return CVec3D( x * a, y * a, z * a );
}

SMF_INLINE_FORCED CVec3D CVec3D::operator/( const float a ) const {
	float inva = 1.0f / a;
	return CVec3D( x * inva, y * inva, z * inva );
}

SMF_INLINE_FORCED CVec3D operator*( const float a, const CVec3D b ) {
	return CVec3D( b.x * a, b.y * a, b.z * a );
}

SMF_INLINE_FORCED CVec3D CVec3D::operator+( const CVec3D &a ) const {
	return CVec3D( x + a.x, y + a.y, z + a.z );
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator+=( const CVec3D &a ) {
	x += a.x;
	y += a.y;
	z += a.z;

	return *this;
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator/=( const CVec3D &a ) {
	x /= a.x;
	y /= a.y;
	z /= a.z;

	return *this;
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator/=( const float a ) {
	float inva = 1.0f / a;
	x *= inva;
	y *= inva;
	z *= inva;

	return *this;
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator-=( const CVec3D &a ) {
	x -= a.x;
	y -= a.y;
	z -= a.z;

	return *this;
}

SMF_INLINE_FORCED CVec3D &CVec3D::operator*=( const float a ) {
	x *= a;
	y *= a;
	z *= a;

	return *this;
}

SMF_INLINE_FORCED bool CVec3D::compare( const CVec3D &a ) const {
	return ( ( x == a.x ) && ( y == a.y ) && ( z == a.z ) );
}

SMF_INLINE_FORCED bool CVec3D::compare( const CVec3D &a, const float epsilon ) const {
	if ( CMath::fabs( x - a.x ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( y - a.y ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( z - a.z ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED bool CVec3D::operator==( const CVec3D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CVec3D::operator!=( const CVec3D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED float CVec3D::normalizeFast() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z;
	invLength = CMath::rSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;
	return invLength * sqrLength;
}

SMF_INLINE_FORCED bool CVec3D::fixDegenerateNormal() {
	if ( x == 0.0f ) {
		if ( y == 0.0f ) {
			if ( z > 0.0f ) {
				if ( z != 1.0f ) {
					z = 1.0f;
					return true;
				}
			} else {
				if ( z != -1.0f ) {
					z = -1.0f;
					return true;
				}
			}
			return false;
		} else if ( z == 0.0f ) {
			if ( y > 0.0f ) {
				if ( y != 1.0f ) {
					y = 1.0f;
					return true;
				}
			} else {
				if ( y != -1.0f ) {
					y = -1.0f;
					return true;
				}
			}
			return false;
		}
	} else if ( y == 0.0f ) {
		if ( z == 0.0f ) {
			if ( x > 0.0f ) {
				if ( x != 1.0f ) {
					x = 1.0f;
					return true;
				}
			} else {
				if ( x != -1.0f ) {
					x = -1.0f;
					return true;
				}
			}
			return false;
		}
	}
	if ( CMath::fabs( x ) == 1.0f ) {
		if ( y != 0.0f || z != 0.0f ) {
			y = z = 0.0f;
			return true;
		}
		return false;
	} else if ( CMath::fabs( y ) == 1.0f ) {
		if ( x != 0.0f || z != 0.0f ) {
			x = z = 0.0f;
			return true;
		}
		return false;
	} else if ( CMath::fabs( z ) == 1.0f ) {
		if ( x != 0.0f || y != 0.0f ) {
			x = y = 0.0f;
			return true;
		}
		return false;
	}
	return false;
}

SMF_INLINE_FORCED bool CVec3D::fixDenormals() {
	bool denormal = false;
	if ( fabs( x ) < 1e-30f ) {
		x = 0.0f;
		denormal = true;
	}
	if ( fabs( y ) < 1e-30f ) {
		y = 0.0f;
		denormal = true;
	}
	if ( fabs( z ) < 1e-30f ) {
		z = 0.0f;
		denormal = true;
	}
	return denormal;
}

SMF_INLINE_FORCED CVec3D CVec3D::cross( const CVec3D &a ) const {
	return CVec3D( y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x );
}

SMF_INLINE_FORCED CVec3D &CVec3D::cross( const CVec3D &a, const CVec3D &b ) {
	x = a.y * b.z - a.z * b.y;
	y = a.z * b.x - a.x * b.z;
	z = a.x * b.y - a.y * b.x;

	return *this;
}

SMF_INLINE_FORCED float CVec3D::getLenght() const {
	return ( float )CMath::sqrt( x * x + y * y + z * z );
}

SMF_INLINE_FORCED float CVec3D::getLengthSqr() const {
	return ( x * x + y * y + z * z );
}

SMF_INLINE_FORCED float CVec3D::getLengthFast() const {
	float sqrLength;

	sqrLength = x * x + y * y + z * z;
	return sqrLength * CMath::rSqrt( sqrLength );
}

SMF_INLINE_FORCED void CVec3D::toNormal() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z;
	invLength = CMath::invSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;

}

SMF_INLINE_FORCED CVec3D CVec3D::normalized() const
{
	CVec3D copy = *this;
	float oldLength = copy.getLenght();
	copy.toNormal();
	//s SMF_ASSERT(oldLength > 0.f && "CVec3D::normalized() failed!");
	//s MARK_UNUSED(oldLength);
	return copy;
}
SMF_INLINE_FORCED CVec3D &CVec3D::truncate( float length ) {
	float length2;
	float ilength;

	if ( !length ) {
		toZero();
	}
	else {
		length2 = getLengthSqr();
		if ( length2 > length * length ) {
			ilength = length * CMath::invSqrt( length2 );
			x *= ilength;
			y *= ilength;
			z *= ilength;
		}
	}

	return *this;
}

SMF_INLINE_FORCED float		CVec3D::getAngle(CVec3D & b){
  CVec3D & a=*this;
 float ab = (a*b);
  float a_b = (  CMath::sqrt( CMath::fabs(a.getLengthSqr() * b.getLengthSqr()) ));
  float param=  ab / a_b;
  float ang4= CMath::acos(param);
  return ang4;
}



SMF_INLINE_FORCED bool CVec3D::isPerpendicular(const CVec3D &other, float epsilon) const
{
	return fabs((*this)*other) <= epsilon * getLenght() * other.getLenght();
}

SMF_INLINE_FORCED bool MUST_USE_RESULT CVec3D::areCollinear(const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, float epsilon)
{
	return (p2-p1).cross(p3-p1).getLengthSqr() <= epsilon;
}
SMF_INLINE_FORCED void CVec3D::clamp( const CVec3D &min, const CVec3D &max ) {
	if ( x < min.x ) {
		x = min.x;
	} else if ( x > max.x ) {
		x = max.x;
	}
	if ( y < min.y ) {
		y = min.y;
	} else if ( y > max.y ) {
		y = max.y;
	}
	if ( z < min.z ) {
		z = min.z;
	} else if ( z > max.z ) {
		z = max.z;
	}
}


SMF_INLINE_FORCED CVec3D CVec3D::Min(float ceil) const
{
	return CVec3D(MIN(x, ceil), MIN(y, ceil), MIN(z, ceil));
}

SMF_INLINE_FORCED CVec3D CVec3D::Min(const CVec3D &ceil) const
{
	return CVec3D(MIN(x, ceil.x), MIN(y, ceil.y), MIN(z, ceil.z));
}

SMF_INLINE_FORCED CVec3D CVec3D::Max(float floor) const
{
	return CVec3D(MAX(x, floor), MAX(y, floor), MAX(z, floor));
}

SMF_INLINE_FORCED CVec3D CVec3D::Max(const CVec3D &floor) const
{
	return CVec3D(MAX(x, floor.x), MAX(y, floor.y), MAX(z, floor.z));
}
SMF_INLINE_FORCED CVec3D CVec3D::clamp( const CVec3D &min, const CVec3D &max ) const {

	return Min(max).Max(min);
}
SMF_INLINE_FORCED void CVec3D::snap() {
	x = CMath::floor( x + 0.5f );
	y = CMath::floor( y + 0.5f );
	z = CMath::floor( z + 0.5f );
}

SMF_INLINE_FORCED void CVec3D::snapInt() {
	x = float( int( x ) );
	y = float( int( y ) );
	z = float( int( z ) );
}

SMF_INLINE_FORCED int CVec3D::getDimension() const {
	return 3;
}

SMF_INLINE_FORCED const CVec2D &CVec3D::toVec2() const {
	return *reinterpret_cast<const CVec2D *>(this);
}

SMF_INLINE_FORCED CVec2D &CVec3D::project_xy() const
{
      CVec2D v;

      v.x = x;
      v.y = y;

      return v;
}


SMF_INLINE_FORCED CVec2D &CVec3D::project_zy() const
    {
      CVec2D v;

      v.x = z;
      v.y = y;

      return v;
    }

SMF_INLINE_FORCED CVec2D &CVec3D::toVec2() {
	return *reinterpret_cast<CVec2D *>(this);
}

SMF_INLINE_FORCED const float *CVec3D::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CVec3D::toFloatPtr() {
	return &x;
}

SMF_INLINE_FORCED CVec3D MUST_USE_RESULT CVec3D::fromScalar(float scalar)
{
	return CVec3D(scalar, scalar, scalar);
}

SMF_INLINE_FORCED void CVec3D::setFromScalar(float scalar)
{
	x = scalar;
	y = scalar;
	z = scalar;
}
SMF_INLINE_FORCED void CVec3D::normalVectors( CVec3D &left, CVec3D &down ) const {
	float d;

	d = x * x + y * y;
	if ( !d ) {
		left[0] = 1;
		left[1] = 0;
		left[2] = 0;
	} else {
		d = CMath::invSqrt( d );
		left[0] = -y * d;
		left[1] = x * d;
		left[2] = 0;
	}
	down = left.cross( *this );
}

SMF_INLINE_FORCED void CVec3D::orthogonalBasis( CVec3D &left, CVec3D &up ) const {
	float l, s;

	if ( CMath::fabs( z ) > 0.7f ) {
		l = y * y + z * z;
		s = CMath::invSqrt( l );
		up[0] = 0;
		up[1] = z * s;
		up[2] = -y * s;
		left[0] = l * s;
		left[1] = -x * up[2];
		left[2] = x * up[1];
	}
	else {
		l = x * x + y * y;
		s = CMath::invSqrt( l );
		left[0] = -y * s;
		left[1] = x * s;
		left[2] = 0;
		up[0] = -z * left[1];
		up[1] = z * left[0];
		up[2] = l * s;
	}
}

SMF_INLINE_FORCED void CVec3D::projectOntoPlane( const CVec3D &normal, const float overBounce ) {
	float backoff;

	backoff = *this * normal;

	if ( overBounce != 1.0 ) {
		if ( backoff < 0 ) {
			backoff *= overBounce;
		} else {
			backoff /= overBounce;
		}
	}

	*this -= backoff * normal;
}

SMF_INLINE_FORCED bool CVec3D::projectAlongPlane( const CVec3D &normal, const float epsilon, const float overBounce ) {
	CVec3D cross;
	float len;

	cross = this->cross( normal ).cross( (*this) );
	// normalize so a fixed epsilon can be used
	cross.toNormal();
	len = normal * cross;
	if ( CMath::fabs( len ) < epsilon ) {
		return false;
	}
	cross *= overBounce * ( normal * (*this) ) / len;
	(*this) -= cross;
	return true;
}


/**
 * \class CVec4D
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de 4  dimensões
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

  */
class SMF_API ALIGNTO16 CVec4D {
public:

		//NAMELESS_UNION_BEGIN // Allow nonstandard nameless struct in union extension on MSC.

	union
	{
		struct
		{
			/// The x component.
			/** A float4 is 16 bytes in size. This element lies in the memory offsets 0-3 of this class. */
			float x;
			/// The y component. [similarOverload: x]
			/** This element is packed to the memory offsets 4-7 of this class. */
			float y;
			/// The z component. [similarOverload: x]
			/** This element is packed to the memory offsets 8-11 of this class. */
			float z;
			/// The w component. [similarOverload: x]
			/** This element is packed to the memory offsets 12-15 of this class. */
			float w;
		};
		sf_m128 v;
	};
//	NAMELESS_UNION_END
					CVec4D();
					explicit CVec4D( const float x, const float y, const float z, const float w );
					/// Constructs a new CVec3D with the value (xyz.x, xyz.y, xyz.z, w).
					/** \see x, y, z, w. */
					CVec4D(const CVec3D &xyz, float w);

	void 			set( const float x, const float y, const float z, const float w );
	void			toZero();
	/// Returns the (x, y) part of this vector.
	CVec2D xy() const;

	/// Returns the (x, y, z) part of this vector.
	CVec3D xyz() const;
	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CVec4D			operator-() const;
	float			operator*( const CVec4D &a ) const;
	CVec4D			operator*( const float a ) const;
	CVec4D			operator/( const float a ) const;
	CVec4D			operator+( const CVec4D &a ) const;
	CVec4D			operator-( const CVec4D &a ) const;
	CVec4D &		operator+=( const CVec4D &a );
	CVec4D &		operator-=( const CVec4D &a );
	CVec4D &		operator/=( const CVec4D &a );
	CVec4D &		operator/=( const float a );
	CVec4D &		operator*=( const float a );

		/// Computes the dot product of the (x, y, z) parts of this and the given CVec4D.
	/** \note This function ignores the w component of this vector (assumes w=0).
		\see Dot4(), Cross3(). */
	float Dot3(const CVec3D &rhs) const;
	float Dot3(const CVec4D &rhs) const;

	/**
	\brief Multiplies the x, y, z and w components of the vector by the given scalar. Note that if w != 0,
	this does NOT scale the length of the homogeneous 3D vector.
	**/
	//friend CVec4D  operator*(float a, CVec4D bb );

	/// Returns true if this vector is equal to the given vector, up to given per-element epsilon.
	bool			compare( const CVec4D &a ) const;							// exact compare, no epsilon
	bool			compare( const CVec4D &a, const float epsilon ) const;		// compare with epsilon
	bool			compare(float x, float y, float z, float w, float epsilon = 1e-3f) const;

	bool			operator==(	const CVec4D &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const CVec4D &a ) const;						// exact compare, no epsilon
	/** Computes the length of this vector.
	\return Sqrt(x*x + y*y + z*z + w*w).
	\see lengthSqr3(), length3(), getLengthSqr(). 
	*/
	float			getLenght() const;
	/** Computes the squared length of this vector.
	 Calling this function is faster than calling Length4(), since this function avoids computing a square root.
		If you only need to compare lengths to each other, but are not interested in the actual length values,
		you can compare by using LengthSq4(), instead of Length4(), since Sqrt() is an order-preserving
		(monotonous and non-decreasing) function.
	\return x*x + y*y + z*z + w*w.
	\see length3(), lengthSqr3(), getLength().
	*/
	float			getLengthSqr() const;
	/** Computes the squared length of the (x, y, z) part of this vector.
	 Calling this function is faster than calling getLength3(), since this function avoids computing a square root.
		If you only need to compare lengths to each other, but are not interested in the actual length values,
		you can compare by using getLengthSq3(), instead of getLength3(), since Sqrt() is an order-preserving
		(monotonous and non-decreasing) function.
		\note This function ignores the w component of this vector.
		\return x*x + y*y + z*z.
		\see getLength3(), getLengthSqr(), getLength(). 
	*/
	float getLengthSqr3() const;

	/** Computes the length of the (x, y, z) part of this vector.
	\note This function ignores the w component of this vector.
	\return Sqrt(x*x + y*y + z*z).
	\see gatLengthSqr3(), gatLengthSqr(), getLength().
	*/
	float getLength3() const;



	void			toNormal();
	float			normalizeFast();		// returns length

	int				getDimension() const;

	const CVec2D &	toVec2() const;
	CVec2D &		toVec2();
	const CVec3D &	toVec3() const;
	CVec3D &		toVec3();
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const sf_m128*   toM128Ptr()const;
	sf_m128*   toM128Ptr();
	const char *	toString( int precision = 2 ) const;

	void			lerp( const CVec4D &v1, const CVec4D &v2, const float l );

	/** Tests if the length of the (x, y, z) part of this vector is one, up to the given epsilon.
	*/
	bool isNormalized3(float epsilonSq = 1e-6f) const;

	/** Returns true if the length of this vector is 1, up to the given epsilon.
	 This function takes into account all the four components of this vector when calculating the norm.
	*/
	bool isNormalized4(float epsilonSq = 1e-6f) const;



	CVec4D(sf_m128 vec):v(vec) {}
//================Static======================
		/// Specifies a compile-time constant CVec4D with value (0, 0, 0, 0).
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D zero;
	static const CVec4D origin;

	/// Specifies a compile-time constant CVec4D with value (1, 1, 1, 1). [similarOverload: zero]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D one;

	/// Specifies a compile-time constant CVec4D with value (1, 0, 0, 0).
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D unitX;

	/// Specifies a compile-time constant CVec4D with value (0, 1, 0, 0). [similarOverload: unitX]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D unitY;

	/// Specifies a compile-time constant CVec4D with value (0, 0, 1, 0). [similarOverload: unitX]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D unitZ;

	/// Specifies a compile-time constant CVec4D with value (0, 0, 0, 1). [similarOverload: unitX]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D unitW;

	/// A compile-time constant CVec4D with value (NaN, NaN, NaN, NaN).
	/** For this constant, each element has the value of quiet NaN, or Not-A-Number.
		\note Never compare a CVec4D to this value! Due to how IEEE floats work, for each float x, both expressions "x == nan" and "x != nan" return false!
			  That is, nothing is equal to NaN, not even NaN itself!
		\note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D nan;

	/// A compile-time constant CVec4D with value (+infinity, +infinity, +infinity, +infinity). [similarOverload: nan]
	/** \note Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! */
	static const CVec4D inf;
	/// Tests if the (x, y, z) parts of two vectors are perpendicular to each other.
	bool isPerpendicular3(const CVec4D &other, float epsilon = 1e-6f) const;

		/// Returns true if the w component of this float4 is either 0 or 1.
	/** This is a required condition for several functions to work correctly.
		*/
	bool isWZeroOrOne(float epsilon =CMath::EPSILON_SuperLow) const;

};


SMF_INLINE_FORCED CVec4D::CVec4D() {
}

SMF_INLINE_FORCED CVec4D::CVec4D( const float x, const float y, const float z, const float w ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}
SMF_INLINE_FORCED CVec4D::CVec4D(const CVec3D &xyz, float w_)
:x(xyz.x), y(xyz.y), z(xyz.z), w(w_)
{
}

SMF_INLINE_FORCED void CVec4D::set( const float x, const float y, const float z, const float w ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

SMF_INLINE_FORCED void CVec4D::toZero() {
	x = y = z = w = 0.0f;
}

SMF_INLINE_FORCED float CVec4D::operator[]( int index ) const {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float& CVec4D::operator[]( int index ) {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED CVec4D CVec4D::operator-() const {
	return CVec4D( -x, -y, -z, -w );
}

SMF_INLINE_FORCED CVec4D CVec4D::operator-( const CVec4D &a ) const {
	return CVec4D( x - a.x, y - a.y, z - a.z, w - a.w );
}

SMF_INLINE_FORCED float CVec4D::operator*( const CVec4D &a ) const {
	return x * a.x + y * a.y + z * a.z + w * a.w;
}

SMF_INLINE_FORCED CVec4D CVec4D::operator*( const float a ) const {
	return CVec4D( x * a, y * a, z * a, w * a );
}

SMF_INLINE_FORCED CVec4D CVec4D::operator/( const float a ) const {
	float inva = 1.0f / a;
	return CVec4D( x * inva, y * inva, z * inva, w * inva );
}

SMF_INLINE_FORCED CVec4D CVec4D::operator+( const CVec4D &a ) const {
	return CVec4D( x + a.x, y + a.y, z + a.z, w + a.w );
}

SMF_INLINE_FORCED CVec4D &CVec4D::operator+=( const CVec4D &a ) {
	x += a.x;
	y += a.y;
	z += a.z;
	w += a.w;

	return *this;
}

SMF_INLINE_FORCED CVec4D &CVec4D::operator/=( const CVec4D &a ) {
	x /= a.x;
	y /= a.y;
	z /= a.z;
	w /= a.w;

	return *this;
}

SMF_INLINE_FORCED CVec4D &CVec4D::operator/=( const float a ) {
	float inva = 1.0f / a;
	x *= inva;
	y *= inva;
	z *= inva;
	w *= inva;

	return *this;
}

SMF_INLINE_FORCED CVec4D &CVec4D::operator-=( const CVec4D &a ) {
	x -= a.x;
	y -= a.y;
	z -= a.z;
	w -= a.w;

	return *this;
}

SMF_INLINE_FORCED CVec4D &CVec4D::operator*=( const float a ) {
	x *= a;
	y *= a;
	z *= a;
	w *= a;

	return *this;
}
SMF_INLINE_FORCED bool CVec4D::isNormalized4(float epsilonSq) const
{
	return CMath::fabs(getLengthSqr()-1.f) <= epsilonSq;
}

SMF_INLINE_FORCED bool CVec4D::isNormalized3(float epsilonSq) const
{
	return CMath::fabs(getLengthSqr3()-1.f) <= epsilonSq;
}
SMF_INLINE_FORCED bool CVec4D::compare( const CVec4D &a ) const {
	return ( ( x == a.x ) && ( y == a.y ) && ( z == a.z ) && w == a.w );
}

SMF_INLINE_FORCED bool CVec4D::compare( const CVec4D &a, const float epsilon ) const {
	if ( CMath::fabs( x - a.x ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( y - a.y ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( z - a.z ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( w - a.w ) > epsilon ) {
		return false;
	}

	return true;
}
SMF_INLINE_FORCED bool CVec4D::compare(float x_, float y_, float z_, float w_, float epsilon) const
{
	return fabs(x - x_) < epsilon &&
		   fabs(y - y_) < epsilon &&
		   fabs(z - z_) < epsilon &&
		   fabs(w - w_) < epsilon;
}

SMF_INLINE_FORCED bool CVec4D::operator==( const CVec4D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CVec4D::operator!=( const CVec4D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED float CVec4D::getLenght() const {
	return ( float )CMath::sqrt( x * x + y * y + z * z + w * w );
}

SMF_INLINE_FORCED float CVec4D::getLengthSqr() const {
	return ( x * x + y * y + z * z + w * w );
}

SMF_INLINE_FORCED void CVec4D::toNormal() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z + w * w;
	invLength = CMath::invSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;
	w *= invLength;

}

SMF_INLINE_FORCED float CVec4D::normalizeFast() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z + w * w;
	invLength = CMath::rSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;
	w *= invLength;
	return invLength * sqrLength;
}

SMF_INLINE_FORCED int CVec4D::getDimension() const {
	return 4;
}

SMF_INLINE_FORCED const CVec2D &CVec4D::toVec2() const {
	return *reinterpret_cast<const CVec2D *>(this);
}

SMF_INLINE_FORCED CVec2D &CVec4D::toVec2() {
	return *reinterpret_cast<CVec2D *>(this);
}

SMF_INLINE_FORCED const CVec3D &CVec4D::toVec3() const {
	return *reinterpret_cast<const CVec3D *>(this);
}

SMF_INLINE_FORCED CVec3D &CVec4D::toVec3() {
	return *reinterpret_cast<CVec3D *>(this);
}

SMF_INLINE_FORCED const float *CVec4D::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CVec4D::toFloatPtr() {
	return &x;
}
SMF_INLINE_FORCED const sf_m128*   CVec4D::toM128Ptr()const{
	return &v;
}
SMF_INLINE_FORCED 	sf_m128*   CVec4D::toM128Ptr(){
	return &v;
}

SMF_INLINE_FORCED bool CVec4D::isPerpendicular3(const CVec4D &other, float epsilon) const
{
	return fabs(this->Dot3(other)) < epsilon;
}

SMF_INLINE_FORCED float CVec4D::Dot3(const CVec3D &rhs) const
{
if(CMath::MATH_AUTOMATIC_SSE){
//s	return dot3_float(v, CVec4D(rhs, 0.f));
}else{
	return x * rhs.x + y * rhs.y + z * rhs.z;
}
}

SMF_INLINE_FORCED float CVec4D::Dot3(const CVec4D &rhs) const
{
	if(CMath::MATH_AUTOMATIC_SSE){
//s		return dot3_float(v, rhs.v);
	}else{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}
}


/**
 * \class CVec5D
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de 2  dimensões
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

  */
class SMF_API CVec5D {
public:
	float			x;
	float			y;
	float			z;
	float			s;
	float			t;

					CVec5D();
					explicit CVec5D( const CVec3D &xyz, const CVec2D &st );
					explicit CVec5D( const float x, const float y, const float z, const float s, const float t );

	float			operator[]( int index ) const;
	float &			operator[]( int index );
	CVec5D &		operator=( const CVec3D &a );
	void			toNormal();
	float			normalizeFast();		// returns length

	int				getDimension() const;

	const CVec3D &	toVec3() const;
	CVec3D &		toVec3();
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	void			lerp( const CVec5D &v1, const CVec5D &v2, const float l );
};

extern CVec5D vec5_origin;
#define vec5_zero vec5_origin

SMF_INLINE_FORCED CVec5D::CVec5D() {
}

SMF_INLINE_FORCED CVec5D::CVec5D( const CVec3D &xyz, const CVec2D &st ) {
	x = xyz.x;
	y = xyz.y;
	z = xyz.z;
	s = st[0];
	t = st[1];
}

SMF_INLINE_FORCED CVec5D::CVec5D( const float x, const float y, const float z, const float s, const float t ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->s = s;
	this->t = t;
}

SMF_INLINE_FORCED float CVec5D::operator[]( int index ) const {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float& CVec5D::operator[]( int index ) {
	return ( &x )[ index ];
}

SMF_INLINE_FORCED CVec5D &CVec5D::operator=( const CVec3D &a ) {
	x = a.x;
	y = a.y;
	z = a.z;
	s = t = 0;
	return *this;
}

SMF_INLINE_FORCED int CVec5D::getDimension() const {
	return 5;
}

SMF_INLINE_FORCED const CVec3D &CVec5D::toVec3() const {
	return *reinterpret_cast<const CVec3D *>(this);
}

SMF_INLINE_FORCED CVec3D &CVec5D::toVec3() {
	return *reinterpret_cast<CVec3D *>(this);
}

SMF_INLINE_FORCED const float *CVec5D::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CVec5D::toFloatPtr() {
	return &x;
}
SMF_INLINE_FORCED void CVec5D::toNormal() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z + s * s + t * t;
	invLength = CMath::invSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;
	s *= invLength;
	t *= invLength;


}

SMF_INLINE_FORCED float CVec5D::normalizeFast() {
	float sqrLength, invLength;

	sqrLength = x * x + y * y + z * z + s * s + t * t;
	invLength = CMath::rSqrt( sqrLength );
	x *= invLength;
	y *= invLength;
	z *= invLength;
	s *= invLength;
	t *= invLength;
	return invLength * sqrLength;
}


/**
 * \class CVec6D
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de 6  dimensões
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

  */
class SMF_API CVec6D {
public:
					CVec6D();
					explicit CVec6D( const float *a );
					explicit CVec6D( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );

	void 			set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );
	void			toZero();

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CVec6D			operator-() const;
	CVec6D			operator*( const float a ) const;
	CVec6D			operator/( const float a ) const;
	float			operator*( const CVec6D &a ) const;
	CVec6D			operator-( const CVec6D &a ) const;
	CVec6D			operator+( const CVec6D &a ) const;
	CVec6D &		operator*=( const float a );
	CVec6D &		operator/=( const float a );
	CVec6D &		operator+=( const CVec6D &a );
	CVec6D &		operator-=( const CVec6D &a );

	friend CVec6D	operator*( const float a, const CVec6D b );

	bool			compare( const CVec6D &a ) const;							// exact compare, no epsilon
	bool			compare( const CVec6D &a, const float epsilon ) const;		// compare with epsilon
	bool			operator==(	const CVec6D &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const CVec6D &a ) const;						// exact compare, no epsilon

	float			getLenght() const;
	float			getLengthSqr() const;
	void			toNormal();
	float			normalizeFast();		// returns length

	int				getDimension() const;

	const CVec3D &	subVec3( int index ) const;
	CVec3D &		subVec3( int index );
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

private:
	float			p[6];
};

extern CVec6D vec6_origin;
#define vec6_zero vec6_origin
extern CVec6D vec6_infinity;

SMF_INLINE_FORCED CVec6D::CVec6D() {
}

SMF_INLINE_FORCED CVec6D::CVec6D( const float *a ) {
	memcpy( p, a, 6 * sizeof( float ) );
}

SMF_INLINE_FORCED CVec6D::CVec6D( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
	p[0] = a1;
	p[1] = a2;
	p[2] = a3;
	p[3] = a4;
	p[4] = a5;
	p[5] = a6;
}

SMF_INLINE_FORCED CVec6D CVec6D::operator-() const {
	return CVec6D( -p[0], -p[1], -p[2], -p[3], -p[4], -p[5] );
}

SMF_INLINE_FORCED float CVec6D::operator[]( const int index ) const {
	return p[index];
}

SMF_INLINE_FORCED float &CVec6D::operator[]( const int index ) {
	return p[index];
}

SMF_INLINE_FORCED CVec6D CVec6D::operator*( const float a ) const {
	return CVec6D( p[0]*a, p[1]*a, p[2]*a, p[3]*a, p[4]*a, p[5]*a );
}

SMF_INLINE_FORCED float CVec6D::operator*( const CVec6D &a ) const {
	return p[0] * a[0] + p[1] * a[1] + p[2] * a[2] + p[3] * a[3] + p[4] * a[4] + p[5] * a[5];
}

SMF_INLINE_FORCED CVec6D CVec6D::operator/( const float a ) const {
	float inva;

	SMF_ASSERT( a != 0.0f );
	inva = 1.0f / a;
	return CVec6D( p[0]*inva, p[1]*inva, p[2]*inva, p[3]*inva, p[4]*inva, p[5]*inva );
}

SMF_INLINE_FORCED CVec6D CVec6D::operator+( const CVec6D &a ) const {
	return CVec6D( p[0] + a[0], p[1] + a[1], p[2] + a[2], p[3] + a[3], p[4] + a[4], p[5] + a[5] );
}

SMF_INLINE_FORCED CVec6D CVec6D::operator-( const CVec6D &a ) const {
	return CVec6D( p[0] - a[0], p[1] - a[1], p[2] - a[2], p[3] - a[3], p[4] - a[4], p[5] - a[5] );
}

SMF_INLINE_FORCED CVec6D &CVec6D::operator*=( const float a ) {
	p[0] *= a;
	p[1] *= a;
	p[2] *= a;
	p[3] *= a;
	p[4] *= a;
	p[5] *= a;
	return *this;
}

SMF_INLINE_FORCED CVec6D &CVec6D::operator/=( const float a ) {
	float inva;

	SMF_ASSERT( a != 0.0f );
	inva = 1.0f / a;
	p[0] *= inva;
	p[1] *= inva;
	p[2] *= inva;
	p[3] *= inva;
	p[4] *= inva;
	p[5] *= inva;
	return *this;
}

SMF_INLINE_FORCED CVec6D &CVec6D::operator+=( const CVec6D &a ) {
	p[0] += a[0];
	p[1] += a[1];
	p[2] += a[2];
	p[3] += a[3];
	p[4] += a[4];
	p[5] += a[5];
	return *this;
}

SMF_INLINE_FORCED CVec6D &CVec6D::operator-=( const CVec6D &a ) {
	p[0] -= a[0];
	p[1] -= a[1];
	p[2] -= a[2];
	p[3] -= a[3];
	p[4] -= a[4];
	p[5] -= a[5];
	return *this;
}

SMF_INLINE_FORCED CVec6D operator*( const float a, const CVec6D b ) {
	return b * a;
}

SMF_INLINE_FORCED bool CVec6D::compare( const CVec6D &a ) const {
	return ( ( p[0] == a[0] ) && ( p[1] == a[1] ) && ( p[2] == a[2] ) &&
			( p[3] == a[3] ) && ( p[4] == a[4] ) && ( p[5] == a[5] ) );
}

SMF_INLINE_FORCED bool CVec6D::compare( const CVec6D &a, const float epsilon ) const {
	if ( CMath::fabs( p[0] - a[0] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[1] - a[1] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[2] - a[2] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[3] - a[3] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[4] - a[4] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[5] - a[5] ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED bool CVec6D::operator==( const CVec6D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CVec6D::operator!=( const CVec6D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CVec6D::set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
	p[0] = a1;
	p[1] = a2;
	p[2] = a3;
	p[3] = a4;
	p[4] = a5;
	p[5] = a6;
}

SMF_INLINE_FORCED void CVec6D::toZero() {
	p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = 0.0f;
}

SMF_INLINE_FORCED float CVec6D::getLenght() const {
	return ( float )CMath::sqrt( p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5] );
}

SMF_INLINE_FORCED float CVec6D::getLengthSqr() const {
	return ( p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5] );
}

SMF_INLINE_FORCED void CVec6D::toNormal() {
	float sqrLength, invLength;

	sqrLength = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5];
	invLength = CMath::invSqrt( sqrLength );
	p[0] *= invLength;
	p[1] *= invLength;
	p[2] *= invLength;
	p[3] *= invLength;
	p[4] *= invLength;
	p[5] *= invLength;

}

SMF_INLINE_FORCED float CVec6D::normalizeFast() {
	float sqrLength, invLength;

	sqrLength = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5];
	invLength = CMath::rSqrt( sqrLength );
	p[0] *= invLength;
	p[1] *= invLength;
	p[2] *= invLength;
	p[3] *= invLength;
	p[4] *= invLength;
	p[5] *= invLength;
	return invLength * sqrLength;
}

SMF_INLINE_FORCED int CVec6D::getDimension() const {
	return 6;
}

SMF_INLINE_FORCED const CVec3D &CVec6D::subVec3( int index ) const {
	return *reinterpret_cast<const CVec3D *>(p + index * 3);
}

SMF_INLINE_FORCED CVec3D &CVec6D::subVec3( int index ) {
	return *reinterpret_cast<CVec3D *>(p + index * 3);
}

SMF_INLINE_FORCED const float *CVec6D::toFloatPtr() const {
	return p;
}

SMF_INLINE_FORCED float *CVec6D::toFloatPtr() {
	return p;
}



/**
 * \class CVecXD
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um Vetor de dimensões arbitrárias
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * \note  A classeCVecXD não pode ser utilizada em multithread devido ao pool de memória temporário
 * \note  O Vetor é armazenado num peáço de memódia alinhado e preenchido de 16bits
 *  \note http://0xc0de.wordpress.com/2008/10/31/alinhamento-alignment-e-preenchimento-padding/
 * \note  The vector lives on 16 sf_u8 aligned and 16 sf_u8 padded memory.
 * \note  due to the temporary memory pool CVecXD cannot be used by multiple threads
 *
 * Contact: Rasputtim@hotmail.com
 *

  **/


class SMF_API CVecXD {
	friend class CMatXD;

public:
					CVecXD();
					explicit CVecXD( int length );
					explicit CVecXD( int length, float *data );
					~CVecXD();

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CVecXD		operator-() const;
	CVecXD &		operator=( const CVecXD &a );
	CVecXD		operator*( const float a ) const;
	CVecXD		operator/( const float a ) const;
	float			operator*( const CVecXD &a ) const;
	CVecXD		operator-( const CVecXD &a ) const;
	CVecXD		operator+( const CVecXD &a ) const;
	CVecXD &		operator*=( const float a );
	CVecXD &		operator/=( const float a );
	CVecXD &		operator+=( const CVecXD &a );
	CVecXD &		operator-=( const CVecXD &a );

	friend CVecXD	operator*( const float a, const CVecXD b );
    /// exact compare, no epsilon
	bool			compare( const CVecXD &a ) const;
	/// compare with epsilon
	bool			compare( const CVecXD &a, const float epsilon ) const;
	/// exact compare, no epsilon
	bool			operator==(	const CVecXD &a ) const;
	/// exact compare, no epsilon
	bool			operator!=(	const CVecXD &a ) const;

	void			setSize( int size );
	void			changeSize( int size, bool makeZero = false );
	int				getSize() const { return size; }
	void			setData( int length, float *data );
	void			toZero();
	void			zero( int length );
	void			random( int seed, float l = 0.0f, float u = 1.0f );
	void			random( int length, int seed, float l = 0.0f, float u = 1.0f );
	void			negate();
	void			clamp( float min, float max );
	CVecXD &		swapElements( int e1, int e2 );

	float			getLenght() const;
	float			getLengthSqr() const;
	CVecXD			toNormal() const;
	float			normalizeSelf();

	int				getDimension() const;

	const CVec3D &	subVec3( int index ) const;
	CVec3D &		subVec3( int index );
	const CVec6D &	subVec6( int index ) const;
	CVec6D &		subVec6( int index );
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;
    ///Shoud be private but matriz use it
	void			setTempSize( int size );
private:
    /// size of the vector
	int				size;
    /// if -1 p points to data set with setData
	int				alloced;
	/// memory the vector is stored
	float *			p;
	/// used to store intermediate results
	static float	temp[VECX_MAX_TEMP+4];
	/// pointer to 16 sf_u8 aligned temporary memory
	static float *	tempPtr;
	/// index into memory pool, wraps around
	static int		tempIndex;



};


SMF_INLINE_FORCED CVecXD::CVecXD() {
	size = alloced = 0;
	p = NULL;
}

SMF_INLINE_FORCED CVecXD::CVecXD( int length ) {
	size = alloced = 0;
	p = NULL;
	setSize( length );
}

SMF_INLINE_FORCED CVecXD::CVecXD( int length, float *data ) {
	size = alloced = 0;
	p = NULL;
	setData( length, data );
}

SMF_INLINE_FORCED CVecXD::~CVecXD() {
	// if not temp memory
	if ( p && ( p < CVecXD::tempPtr || p >= CVecXD::tempPtr + VECX_MAX_TEMP ) && alloced != -1 ) {
		mem_Free16( p );
	}
}

SMF_INLINE_FORCED float CVecXD::operator[]( const int index ) const {
	SMF_ASSERT( index >= 0 && index < size );
	return p[index];
}

SMF_INLINE_FORCED float &CVecXD::operator[]( const int index ) {
	SMF_ASSERT( index >= 0 && index < size );
	return p[index];
}

SMF_INLINE_FORCED CVecXD CVecXD::operator-() const {
	int i;
	CVecXD m;

	m.setTempSize( size );
	for ( i = 0; i < size; i++ ) {
		m.p[i] = -p[i];
	}
	return m;
}

SMF_INLINE_FORCED CVecXD &CVecXD::operator=( const CVecXD &a ) {
	setSize( a.size );
#ifdef VECX_SIMD
	SIMDProcessor->copy16( p, a.p, a.size );
#else
	memcpy( p, a.p, a.size * sizeof( float ) );
#endif
	CVecXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CVecXD CVecXD::operator+( const CVecXD &a ) const {
	CVecXD m;

	SMF_ASSERT( size == a.size );
	m.setTempSize( size );
#ifdef VECX_SIMD
	SIMDProcessor->add16( m.p, p, a.p, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		m.p[i] = p[i] + a.p[i];
	}
#endif
	return m;
}

SMF_INLINE_FORCED CVecXD CVecXD::operator-( const CVecXD &a ) const {
	CVecXD m;

	SMF_ASSERT( size == a.size );
	m.setTempSize( size );
#ifdef VECX_SIMD
	SIMDProcessor->sub16( m.p, p, a.p, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		m.p[i] = p[i] - a.p[i];
	}
#endif
	return m;
}

SMF_INLINE_FORCED CVecXD &CVecXD::operator+=( const CVecXD &a ) {
	SMF_ASSERT( size == a.size );
#ifdef VECX_SIMD
	SIMDProcessor->addAssign16( p, a.p, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		p[i] += a.p[i];
	}
#endif
	CVecXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CVecXD &CVecXD::operator-=( const CVecXD &a ) {
	SMF_ASSERT( size == a.size );
#ifdef VECX_SIMD
	SIMDProcessor->subAssign16( p, a.p, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		p[i] -= a.p[i];
	}
#endif
	CVecXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CVecXD CVecXD::operator*( const float a ) const {
	CVecXD m;

	m.setTempSize( size );
#ifdef VECX_SIMD
	SIMDProcessor->mul16( m.p, p, a, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		m.p[i] = p[i] * a;
	}
#endif
	return m;
}

SMF_INLINE_FORCED CVecXD &CVecXD::operator*=( const float a ) {
#ifdef VECX_SIMD
	SIMDProcessor->mulAssign16( p, a, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		p[i] *= a;
	}
#endif
	return *this;
}

SMF_INLINE_FORCED CVecXD CVecXD::operator/( const float a ) const {
	SMF_ASSERT( a != 0.0f );
	return (*this) * ( 1.0f / a );
}

SMF_INLINE_FORCED CVecXD &CVecXD::operator/=( const float a ) {
	SMF_ASSERT( a != 0.0f );
	(*this) *= ( 1.0f / a );
	return *this;
}

SMF_INLINE_FORCED CVecXD operator*( const float a, const CVecXD b ) {
	return b * a;
}

SMF_INLINE_FORCED float CVecXD::operator*( const CVecXD &a ) const {
	int i;
	float sum = 0.0f;

	SMF_ASSERT( size == a.size );
	for ( i = 0; i < size; i++ ) {
		sum += p[i] * a.p[i];
	}
	return sum;
}

SMF_INLINE_FORCED bool CVecXD::compare( const CVecXD &a ) const {
	int i;

	SMF_ASSERT( size == a.size );
	for ( i = 0; i < size; i++ ) {
		if ( p[i] != a.p[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CVecXD::compare( const CVecXD &a, const float epsilon ) const {
	int i;

	SMF_ASSERT( size == a.size );
	for ( i = 0; i < size; i++ ) {
		if ( CMath::fabs( p[i] - a.p[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CVecXD::operator==( const CVecXD &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CVecXD::operator!=( const CVecXD &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CVecXD::setSize( int newSize ) {
	int alloc = ( newSize + 3 ) & ~3;
	if ( alloc > alloced && alloced != -1 ) {
		if ( p ) {
			mem_Free16( p );
		}
		p = (float *) mem_Alloc16( alloc * sizeof( float ) );
		alloced = alloc;
	}
	size = newSize;
	VECX_CLEAREND();
}

SMF_INLINE_FORCED void CVecXD::changeSize( int newSize, bool makeZero ) {
	int alloc = ( newSize + 3 ) & ~3;
	if ( alloc > alloced && alloced != -1 ) {
		float *oldVec = p;
		p = (float *) mem_Alloc16( alloc * sizeof( float ) );
		alloced = alloc;
		if ( oldVec ) {
			for ( int i = 0; i < size; i++ ) {
				p[i] = oldVec[i];
			}
			mem_Free16( oldVec );
		}
		if ( makeZero ) {
			// zero any new elements
			for ( int i = size; i < newSize; i++ ) {
				p[i] = 0.0f;
			}
		}
	}
	size = newSize;
	VECX_CLEAREND();
}

SMF_INLINE_FORCED void CVecXD::setTempSize( int newSize ) {

	size = newSize;
	alloced = ( newSize + 3 ) & ~3;
	SMF_ASSERT( alloced < VECX_MAX_TEMP );
	if ( CVecXD::tempIndex + alloced > VECX_MAX_TEMP ) {
		CVecXD::tempIndex = 0;
	}
	p = CVecXD::tempPtr + CVecXD::tempIndex;
	CVecXD::tempIndex += alloced;
	VECX_CLEAREND();
}

SMF_INLINE_FORCED void CVecXD::setData( int length, float *data ) {
	if ( p && ( p < CVecXD::tempPtr || p >= CVecXD::tempPtr + VECX_MAX_TEMP ) && alloced != -1 ) {
		mem_Free16( p );
	}

#ifdef _WIN32
	SMF_ASSERT( ( ( (int) data ) & 15 ) == 0 ); // data must be 16 sf_u8 aligned
#else
	SMF_ASSERT( ( ( (uintptr_t) data ) & 15 ) == 0 ); // data must be 16 byte aligned
#endif


	p = data;
	size = length;
	alloced = -1;
	VECX_CLEAREND();
}

SMF_INLINE_FORCED void CVecXD::toZero() {
#ifdef VECX_SIMD
	SIMDProcessor->zero16( p, size );
#else
	memset( p, 0, size * sizeof( float ) );
#endif
}

SMF_INLINE_FORCED void CVecXD::zero( int length ) {
	setSize( length );
#ifdef VECX_SIMD
	SIMDProcessor->zero16( p, length );
#else
	memset( p, 0, size * sizeof( float ) );
#endif
}

SMF_INLINE_FORCED void CVecXD::random( int seed, float l, float u ) {
	int i;
	float c;
	CRandom rnd( seed );

	c = u - l;
	for ( i = 0; i < size; i++ ) {
		p[i] = l + rnd.randomFloat() * c;
	}
}

SMF_INLINE_FORCED void CVecXD::random( int length, int seed, float l, float u ) {
	int i;
	float c;
	CRandom rnd( seed );

	setSize( length );
	c = u - l;
	for ( i = 0; i < size; i++ ) {
		p[i] = l + rnd.randomFloat() * c;
	}
}

SMF_INLINE_FORCED void CVecXD::negate() {
#ifdef VECX_SIMD
	SIMDProcessor->negate16( p, size );
#else
	int i;
	for ( i = 0; i < size; i++ ) {
		p[i] = -p[i];
	}
#endif
}

SMF_INLINE_FORCED void CVecXD::clamp( float min, float max ) {
	int i;
	for ( i = 0; i < size; i++ ) {
		if ( p[i] < min ) {
			p[i] = min;
		} else if ( p[i] > max ) {
			p[i] = max;
		}
	}
}

SMF_INLINE_FORCED CVecXD &CVecXD::swapElements( int e1, int e2 ) {
	float tmp;
	tmp = p[e1];
	p[e1] = p[e2];
	p[e2] = tmp;
	return *this;
}

SMF_INLINE_FORCED float CVecXD::getLenght() const {
	int i;
	float sum = 0.0f;

	for ( i = 0; i < size; i++ ) {
		sum += p[i] * p[i];
	}
	return CMath::sqrt( sum );
}

SMF_INLINE_FORCED float CVecXD::getLengthSqr() const {
	int i;
	float sum = 0.0f;

	for ( i = 0; i < size; i++ ) {
		sum += p[i] * p[i];
	}
	return sum;
}

SMF_INLINE_FORCED CVecXD CVecXD::toNormal() const {
	int i;
	CVecXD m;
	float invSqrt, sum = 0.0f;

	m.setTempSize( size );
	for ( i = 0; i < size; i++ ) {
		sum += p[i] * p[i];
	}
	invSqrt = CMath::invSqrt( sum );
	for ( i = 0; i < size; i++ ) {
		m.p[i] = p[i] * invSqrt;
	}
	return m;
}

SMF_INLINE_FORCED float CVecXD::normalizeSelf() {
	float invSqrt, sum = 0.0f;
	int i;
	for ( i = 0; i < size; i++ ) {
		sum += p[i] * p[i];
	}
	invSqrt = CMath::invSqrt( sum );
	for ( i = 0; i < size; i++ ) {
		p[i] *= invSqrt;
	}
	return invSqrt * sum;
}

SMF_INLINE_FORCED int CVecXD::getDimension() const {
	return size;
}

SMF_INLINE_FORCED CVec3D &CVecXD::subVec3( int index ) {
	SMF_ASSERT( index >= 0 && index * 3 + 3 <= size );
	return *reinterpret_cast<CVec3D *>(p + index * 3);
}

SMF_INLINE_FORCED const CVec3D &CVecXD::subVec3( int index ) const {
	SMF_ASSERT( index >= 0 && index * 3 + 3 <= size );
	return *reinterpret_cast<const CVec3D *>(p + index * 3);
}

SMF_INLINE_FORCED CVec6D &CVecXD::subVec6( int index ) {
	SMF_ASSERT( index >= 0 && index * 6 + 6 <= size );
	return *reinterpret_cast<CVec6D *>(p + index * 6);
}

SMF_INLINE_FORCED const CVec6D &CVecXD::subVec6( int index ) const {
	SMF_ASSERT( index >= 0 && index * 6 + 6 <= size );
	return *reinterpret_cast<const CVec6D *>(p + index * 6);
}

SMF_INLINE_FORCED const float *CVecXD::toFloatPtr() const {
	return p;
}

SMF_INLINE_FORCED float *CVecXD::toFloatPtr() {
	return p;
}



//===============================================================
/**
 * \class 	CPolar3D

 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Coordenadas Polares 3D.
 *
 * \elseif us_en
 * \brief 	3D Polar Coordinates.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPolar3D {
public:
	float			radius, theta, phi;

					CPolar3D();
					explicit CPolar3D( const float radius, const float theta, const float phi );

	void 			set( const float radius, const float theta, const float phi );

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CPolar3D		operator-() const;
	CPolar3D &		operator=( const CPolar3D &a );

	CVec3D			toVec3() const;
};

SMF_INLINE_FORCED CPolar3D::CPolar3D() {
}

SMF_INLINE_FORCED CPolar3D::CPolar3D( const float radius, const float theta, const float phi ) {
	SMF_ASSERT( radius > 0 );
	this->radius = radius;
	this->theta = theta;
	this->phi = phi;
}

SMF_INLINE_FORCED void CPolar3D::set( const float radius, const float theta, const float phi ) {
	SMF_ASSERT( radius > 0 );
	this->radius = radius;
	this->theta = theta;
	this->phi = phi;
}

SMF_INLINE_FORCED float CPolar3D::operator[]( const int index ) const {
	return ( &radius )[ index ];
}

SMF_INLINE_FORCED float &CPolar3D::operator[]( const int index ) {
	return ( &radius )[ index ];
}

SMF_INLINE_FORCED CPolar3D CPolar3D::operator-() const {
	return CPolar3D( radius, -theta, -phi );
}

SMF_INLINE_FORCED CPolar3D &CPolar3D::operator=( const CPolar3D &a ) {
	radius = a.radius;
	theta = a.theta;
	phi = a.phi;
	return *this;
}

SMF_INLINE_FORCED CVec3D CPolar3D::toVec3() const {
	float sp, cp, st, ct;
	CMath::sincos( phi, sp, cp );
	CMath::sincos( theta, st, ct );
 	return CVec3D( cp * radius * ct, cp * radius * st, radius * sp );
}


/*
===============================================================================

	Old 3D vector macros, should no longer be used.

===============================================================================
*/

#define DotProduct( a, b)			((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define VectorSubtract( a, b, c )	((c)[0]=(a)[0]-(b)[0],(c)[1]=(a)[1]-(b)[1],(c)[2]=(a)[2]-(b)[2])
#define VectorAdd( a, b, c )		((c)[0]=(a)[0]+(b)[0],(c)[1]=(a)[1]+(b)[1],(c)[2]=(a)[2]+(b)[2])
#define	VectorScale( v, s, o )		((o)[0]=(v)[0]*(s),(o)[1]=(v)[1]*(s),(o)[2]=(v)[2]*(s))
#define	VectorMA( v, s, b, o )		((o)[0]=(v)[0]+(b)[0]*(s),(o)[1]=(v)[1]+(b)[1]*(s),(o)[2]=(v)[2]+(b)[2]*(s))
#define VectorCopy( a, b )			((b)[0]=(a)[0],(b)[1]=(a)[1],(b)[2]=(a)[2])
} //end MATH
} //end SMF
#endif /* !__SMF_MATH_VECTOR_H__ */
