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

#ifndef _SMF__MATH_PLANE_H__
#define _SMF__MATH_PLANE_H__
#include "../SMF_Config.h"
#include "SMF_GeoDefs.h"
#include "../math/SMF_Vector.h"
#include "../math/SMF_Matriz.h"

namespace SMF {
namespace MATH{
class CVec3D;
class CMat3D;

}
namespace GEO{

using namespace MATH;

class CPolygon;
class CCircle;
class CPolyhedron;
/**
 * \class CPlane
 *
 * \ingroup SMF_Geometric
 * \image html pics/plane.png
 * \brief Implementa um Plano 3D com a equação geral: (a * x) + (b * y) + (c * z) + d = 0
 *  sendo P=(x0, y0, z0) um ponto qualquer no plano e
 *  sendo n->=(a, b, c) um vetor ortogonal ao plano (vetor normal), onde a sua raiz é o ponto (x,y,z) e a cabeça é o ponto (a,b,c) e
 *  sendo d = -ax0-by0-cz0
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 */
class SMF_API CPlane {
public:
					CPlane();
					/**
					\brief inicializa um plano 3D-> a * x + b * y + c * z + d = 0
					\param a parâmetro a da equação do plano
					\param b parâmetro b da equação do plano
					\param c parâmetro c da equação do plano
					\param d parâmetro d da equação do plano
					\note a,b,c pontos do vetor ortogonal ao plano
					**/
					CPlane( float a, float b, float c, float d );
					/**
					\brief inicializa um plano 3D-> a * x + b * y + c * z + d = 0
					\param normal Vetor normal que contem os parâmetros a, b e c da equação do plano
					\param d parâmetro d da equação do plano
					\note d = -dist
					**/
					CPlane( const CVec3D &normal, const float distPoint );
					/// Constructs a plane by specifying a single point on the plane, and the surface normal.
					/** \param normal The direction the plane is facing. This vector must have been normalized in advance.
					\see set(). */
					 CPlane(const CVec3D &normal, const CVec3D &distPoint );
					/// Sets this plane by specifying three points on the plane.
					/** The normal of the plane will point to the halfspace from which the points are observed to be oriented in
					counter-clockwise order.
					\note The points v1, v2 and v3 must not all lie on the same line. */

	/// Constructs a plane by specifying three points on the plane.
	/** The normal of the plane will point to
		the halfspace from which the points are observed to be oriented in counter-clockwise order.
		\note The points v1, v2 and v3 must not all lie on the same line.
		\see set(). */
	CPlane(const CVec3D &v1, const CVec3D &v2, const CVec3D &v3);

					/// Constructs a plane by specifying a line that lies on the plane, and the plane normal.
	/** \param line The line object that is to be contained in the newly constructed plane.
		\param normal The direction the plane if facing. This vector must have been normalized in advance. The normal
			of the line must not be collinear with the direction of this normal. If a line segment is specified, the line
			segment must not be degenerate. */
	CPlane(const CRay &line, const CVec3D &normal);
	CPlane(const CLine &line, const CVec3D &normal);
	CPlane(const CLineSegment &line, const CVec3D &normal);

					 void set(const CVec3D &v1, const CVec3D &v2, const CVec3D &v3);
	/// Sets this plane by specifying a single point on the plane, and the surface normal.
	/** \param normal The direction the plane is facing. This vector must have been normalized in advance. */
	void set(const CVec3D &point, const CVec3D &normal);


	float			operator[]( int index ) const;
	float &			operator[]( int index );
	/**
	\brief inverte o plano
	**/
	CPlane			operator-() const;						// flips plane
	CPlane &		operator=( const CVec3D &v );		// sets normal and sets CPlane::d to zero
	/**
	\brief soma as quações dos planos
	**/
	CPlane			operator+( const CPlane &p ) const;	// add plane equations
	/**
	\brief subtrai as quações dos planos
	**/
	CPlane			operator-( const CPlane &p ) const;	// subtract plane equations

	CPlane &		operator*=( const CMat3D &m );			// getNormal() *= m

	bool			compare( const CPlane &p ) const;						// exact compare, no epsilon
	bool			compare( const CPlane &p, const float epsilon ) const;	// compare with epsilon
	bool			compare( const CPlane &p, const float normalEps, const float distEps ) const;	// compare with epsilon
	bool			operator==(	const CPlane &p ) const;					// exact compare, no epsilon
	bool			operator!=(	const CPlane &p ) const;					// exact compare, no epsilon
	/**
	\brief zera a equação do plano
	**/
	void			toZero();							// zero plane
	/**
	\brief seta o vetor normal do plano (a,b,c)
	**/
	void			setNormal( const CVec3D &normal );		// sets the normal
	/**
	\brief retorna uma referência para o vetor normal do plano (a,b,c)
	**/
	const CVec3D &	getNormal() const;					// reference to const normal
	/**
	\brief retorna uma referência para o vetor normal do plano (a,b,c)
	**/
	CVec3D &		getNormal();							// reference to normal
	/**
	\brief normaliza o plano: |a, b, c| == 1. não ajusta d
	\param fixDegenerate corrige, ou não normal degenerada
	\return o comprimento
	**/
	float			normalize( bool fixDegenerate = true );	// only normalizes the plane normal, does not adjust d
	/**
	\brief corrige o vetor normal, caso ele esteja degenerado
	\return verdadeiro, se foi feita alguma correção
	**/
	bool			fixDegenerateNormal();			// fix degenerate normal
	/**
	\return verdadeiro se alguma degeneração foi corrigida
	**/
	bool			fixDegeneracies( float distEpsilon );	// fix degenerate normal and dist
	/**
	\brief retorna o parâmetro d do plano ((a * x) + (b * y) + (c * z) + d = 0)
	**/
	float			getDist() const;						// returns: -d
	/**
	\brief seta o parâmetro d : d=-dist
	\param dist d= - dist
	**/
	void			setDist( const float dist );			// sets: d = -dist
	/**
	\brief retorna o tipo de plano
	\return retorna o tipo do plano conforme abaixo:
	// plane types
	* PLANETYPE_X					0
	* PLANETYPE_Y					1
	* PLANETYPE_Z					2
	* PLANETYPE_NEGX				3
	* PLANETYPE_NEGY				4
	* PLANETYPE_NEGZ				5
	* PLANETYPE_TRUEAXIAL			6	// all types < 6 are true axial planes
	* PLANETYPE_ZEROX				6
	* PLANETYPE_ZEROY				7
	* PLANETYPE_ZEROZ				8
	* PLANETYPE_NONAXIAL			9
	**/
	int				getType() const;						// returns plane type
	/**
	\brief cria um plano através de três pontos, p1, p2, p3, sendo p2 o ponto central
	\param p1 ponto 1 qualquer que pertence ao plano
	\param p2 ponto 2 qualquer que pertence ao plano
	\param p3 ponto 3 qualquer que pertence ao plano
	\param fixDegenerate se deve ou não corrigir a matriz caso haja degeneração
	\return retorna falso se não conseguir calcular d
	**/
	bool			fromPoints( const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, bool fixDegenerate = true );
	/**
	\brief cria um plano através de dois vetores pertecentes ao plano e um ponto p
	\param dir1 vetor 1 qualquer que pertence ao plano
	\param dir2 vetor 2 qualquer que pertence ao plano
	\param p ponto  qualquer que pertence ao plano
	\param fixDegenerate se deve ou não corrigir a matriz caso haja degeneração
	\return retorna falso se não conseguir calcular d
	**/

	bool			fromVecs( const CVec3D &dir1, const CVec3D &dir2, const CVec3D &p, bool fixDegenerate = true );
	/**
	\brief cria um plano através do vetor Normal e de um ponto p
	\param dir1 vetor normal ao plano
	\param p ponto  qualquer que pertence ao plano
	\param fixDegenerate se deve ou não corrigir a matriz caso haja degeneração
	\return retorna falso se não conseguir calcular d
	**/

	bool			fromPointNormal( const CVec3D &dir1,  const CVec3D &p, bool fixDegenerate = true );

	void			fitThroughPoint( const CVec3D &p );	// assumes normal is valid
	bool			heightFit( const CVec3D *points, const int numPoints );
	CPlane			translate( const CVec3D &translation ) const;
	CPlane &		translateSelf( const CVec3D &translation );
	CPlane			rotate( const CVec3D &origin, const CMat3D &axis ) const;
	CPlane &		rotateSelf( const CVec3D &origin, const CMat3D &axis );
	/**
	Computes the dot product of a plane and a vector.
	Given a plane (a, b, c, d) and a 3-D vector (x, y, z), this method's return value is a*x + b*y + c*z + d*1.
	The distance method is useful for determining the plane's relationship with a coordinate in 3-D space.
	**/
	float			getDistance( const CVec3D &v ) const;
	/**
	Computes the dot product of a plane and a vector.
	Given a plane (a, b, c, d) and a 4-D vector (x, y, z, w), this method's return value is a*x + b*y + c*z + d*w. The dot method is useful for determining the plane's relationship with a homogeneous coordinate.
	For example, it can determine whether a coordinate is on a particular plane, or on which side of a particular plane a coordinate lies.**/
	float			getDistance( const CVec4D &v ) const;
	/**
	Given a plane (a, b, c, d) and a 3-D vector (x, y, z), this method's return value is a*x + b*y + c*z + d*0.
	The DotNormal method is useful for calculating the angle between the normal of the plane and another normal.
	\return value that represents the dot product of the plane and the 3-D vector.
	\note quando d = 0, o plano passa pela origem 0 = (0,0,0).
	**/
	float			getDistanceNormal( const CVec3D &v ) const;

	int				side( const CVec3D &v, const float epsilon = 0.0f ) const;
	/**
	\brief   verifica se a reta delimitada pelos pontos start e end intercepta o plano.
	\param start ponto inicial da reta
	\param end ponto final da reta
	\return   verdadeiro se a reta intercepta o plano e falso caso contrário.
	**/
	bool			lineIntersection( const CVec3D &start, const CVec3D &end ) const;
					// intersection point is start + dir * scale
	bool			rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale ) const;
	/**
	\brief verifica se dois planos se interceptam
	\param plane plano que se dejeja verificar a intersecção
	\param statPoint
	\param dir
	\return true se os planos se interceptam e false caso contrário
	**/
	bool			planeIntersection( const CPlane &plane, CVec3D &startPoint, CVec3D &dir ) const;
	/**
	\brief retorna a dimenção da classe
	\return dimenção da classe
	**/
	int				getDimension() const;

	const CVec4D &	ToVec4() const;
	CVec4D &		ToVec4();
	/**
	\brief retorna um ponteiro para os dados da classe (a,b,c,d)
	\return ponteiro para os dados da classe
	**/
	const float *	toFloatPtr() const;
	/**
	\brief retorna um ponteiro para os dados da classe (a,b,c,d)
	\return ponteiro para os dados da classe
	**/
	float *			toFloatPtr();
	/**
	\brief retorna uma string com os dados do plano
	**/
	const char *	toString( int precision = 2 ) const;

	//=========Distances====================
	/// Returns the distance of this plane to the given object.
	/** If the given object intersects or lies in this plane, then the returned distance is zero.
		\note This function always returns a positive distance value, even when the given object lies on the negative side
			of this plane. See the signedDistance() function to produce a distance value that differentiates between the
			front and back sides of this plane.
		\see signedDistance(), intersects(), contains(). */
	float distance(const CVec3D &point) const;
	float distance(const CLineSegment &lineSegment) const;
	float distance(const CSphere &sphere) const;
		/// Returns the signed distance of this plane to the given point.
	/** If this function returns a negative value, the given point lies in the negative halfspace of this plane.
		Conversely, if a positive value is returned, then the given point lies in the positive halfspace of this plane.
		\see distance(), isOnPositiveSide(), areOnSameSide(). */
	float signedDistance(const CVec3D &point) const;

	float signedDistance(const CAABBox &aabb) const;
	float signedDistance(const COBBox &obb) const;
//	float signedDistance(const CCircle &circle) const;
	float signedDistance(const CLine &line) const;
	float signedDistance(const CLineSegment &lineSegment) const;
	float signedDistance(const CRay &ray) const;
//	float signedDistance(const CPlane &plane) const;
	float signedDistance(const CPolygon &polygon) const;
	float signedDistance(const CPolyhedron &polyhedron) const;
	float signedDistance(const CSphere &sphere) const;
	float signedDistance(const CTriangle &triangle) const;

		/// Projects the given object onto this plane orthographically.
	/** \note This mapping can be expressed as a CMatJoint3x4 matrix operation. See the OrthoProjection() function.
		\see OrthoProjection(). */
	CVec3D project(const CVec3D &point) const;
	CLineSegment project(const CLineSegment &lineSegment) const;
	/** \param nonDegenerate [out] If the line or ray is perpendicular to the plane, the projection is
		a single point. In that case, the .pos parameter of the returned object will specify the point
		location, the .dir parameter of the object will be undefined and the nonDegenerate pointer will be
		set to false. This pointer may be null. */
	CLine project(const CLine &line, bool *nonDegenerate) const;
	CRay project(const CRay &ray, bool *nonDegenerate) const;

	CTriangle project(const CTriangle &triangle) const;
	CPolygon project(const CPolygon &polygon) const;


		/// Computes the intersection of three planes.
	/** This function computes the intersection of this plane, and the given two planes.
		\param outLine [out] If the three planes are configured in such a manner that intersection is a line,
			this parameter receives the line of intersection. This pointer may be null.
		\param outPoint [out] If the three planes are configured in such a manner that the interesction is a point,
			this parameter receives the point of intersection. This pointer may be null.
		\bug This function never outputs outLine.
		\return True if the intersection was a point, in which case outPoint is written to.
		\see contains(), distance(), closestPoint(). */
	bool intersects(const CPlane &plane, const CPlane &plane2, CLine *outLine = 0, CVec3D *outPoint = 0) const;

	/**
	\brief Tests whether this plane and the given object intersect.
	\param outLine [out] The intersection of two planes forms a line. If an intersection occurs, this parameter will receive
		the line of intersection. This pointer may be null.
	\return True if the given object intersects with this plane. */
	bool intersects(const CPlane &plane, CLine *outLine = 0) const;
	/** \param d [out] If specified, this parameter will receive the parametric distance of
			the intersection point along the line object. Use the getPoint(d) function of the line class
			to get the actual point of intersection. This pointer may be null. */
	bool intersects(const CRay &ray, float *d = 0) const;
	bool intersects(const CLine &line, float *d = 0) const;
	bool intersects(const CLineSegment &lineSegment, float *d = 0) const;
	bool intersects(const CSphere &sphere) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox  &obb) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	/// \todo add a version of Plane-Triangle intersection which returns the line segment of intersection.
	bool intersects(const CTriangle &triangle) const;
	/// Tests if this plane intersects with the given circle.
	/** \param pt1 [out] If specified, receives the first point of intersection. This pointer may be null.
		\param pt2 [out] If specified, receives the second point of intersection. This pointer may be null.
		\return The number of intersections that occurred: 0, 1 or 2. */
	int intersects(const CCircle &circle, CVec3D *pt1, CVec3D *pt2) const;
	int intersects(const CCircle &circle) const;
		/// Computes the intersection of a line and a plane.
	/** \param planeNormal The plane normal direction vector. This vector can be unnormalized.
		\param planeD The distance parameter of the plane equation.
		\param linePos The starting point of the line.
		\param lineDir The line direction vector. This vector does not need to be normalized.
		\param t [out] If this function returns true, this parameter will receive the distance along the line where intersection occurs.
					That is, the point lineStart + t * lineDir will be the intersection point. Note that if |lineDir| != 1,
					then t will not contain the real distance, but one scaled to the units of lineDir.
		\return If an intersection occurs, this function returns true. */
	static bool intersectLinePlane(const CVec3D &planeNormal, float planeD, const CVec3D &linePos, const CVec3D &lineDir, float &t);

	/// Tests if two planes are parallel.
	/** \see SetEquals(), compare(), DihedralAngle(). */
	bool isParallel(const CPlane &plane, float epsilon =CMath::EPSILON_SuperLow) const;
	/** 
	\brief Tests if the given direction vector points towards the positive side of this plane.
	 \param directionVector The direction vector to compare with the normal of this plane. This vector
	may be unnormalized.
	\see isOnPositiveSide. */
	bool isInPositiveDirection(const CVec3D &directionVector) const;

	/// Tests if the given point lies on the positive side of this plane.
	/** A plane divides the space in three sets: the negative halfspace, the plane itself, and the positive halfspace.
		The normal vector of the plane points towards the positive halfspace.
		\return This function returns true if the given point lies either on this plane itself, or in the positive
			halfspace of this plane.
		\see isInPositiveDirection, areOnSameSide(), distance(), signedDistance(). */
	bool isOnPositiveSide(const CVec3D &point) const;

	/// Performs a Triangle-Plane intersection test.
	/** \return This function returns the value 1 if the whole triangle is on the positive side of this plane, the
			value -1 if the whole triangle lies in the negative side of this plane, and 0 if the triangle intersects this plane.
		\see intersects(), areOnSameSide(), distance(), signedDistance(), contains(). */
	int examineSide(const CTriangle &triangle) const;

	/// Tests if two points are on the same side of this plane.
	/** \return This function returns true if both p1 and p2 are on the positive side or this plane, or if both p1 and p2
			are on the negative side of this plane.
		\see isOnPositiveSide(), distance(), signedDistance(). */
	bool areOnSameSide(const CVec3D &p1, const CVec3D &p2) const;

	/// Returns the distance of this plane to the given object.
	/** If the given object intersects or lies in this plane, then the returned distance is zero.
		\note This function always returns a positive distance value, even when the given object lies on the negative side
			of this plane. See the signedDistance() function to produce a distance value that differentiates between the
			front and back sides of this plane.
		\see signedDistance(), intersects(), contains(). */


public:
	float			a;
	float			b;
	float			c;
	float			d;
};

extern CPlane plane_origin;
#define plane_zero plane_origin

SMF_INLINE_FORCED CPlane::CPlane() {
}

SMF_INLINE_FORCED CPlane::CPlane( float a, float b, float c, float d ) {
	this->a = a;
	this->b = b;
	this->c = c;
	this->d = d;
}

SMF_INLINE_FORCED CPlane::CPlane( const CVec3D &normal, const float dist ) {
	this->a = normal.x;
	this->b = normal.y;
	this->c = normal.z;
	this->d = -dist;
}
SMF_INLINE_FORCED CPlane::CPlane( const CVec3D &normal_, const CVec3D &point)
{
	this->a = normal_.x;
	this->b = normal_.y;
	this->c = normal_.z;

	SMF_ASSERT(getNormal().isNormalized());
	d = point* this->getNormal();

	//mathassert(equalsAbs(signedDistance(point), 0.f, 0.01f));
	//mathassert(equalsAbs(signedDistance(point + normal_), 1.f, 0.01f));

}




SMF_INLINE_FORCED CPlane::CPlane(const CVec3D &v1, const CVec3D &v2, const CVec3D &v3){
	set(v1,v2,v3);
}


SMF_INLINE_FORCED void CPlane::set(const CVec3D &v1, const CVec3D &v2, const CVec3D &v3)
{
	setNormal(((v2-v1).cross(v3-v1)).normalized());
	d = (v1* getNormal());

#ifdef MATH_ASSERT_CORRECTNESS
	float d2 = dot(v2, normal);
	float d3 = dot(v3, normal);
	mathassert(equalsAbs(d, d2, 1e-2f));
	mathassert(equalsAbs(d, d3, 1e-2f));
#endif
}


SMF_INLINE_FORCED void CPlane::set(const CVec3D &point, const CVec3D &normal_)
{
	setNormal(normal_);
	SMF_ASSERT(getNormal().isNormalized());
	d = (point* getNormal());
#ifdef TEST_FOR_CORRECTNESS
	SMF_ASSERT(CMath::equalsAbs(signedDistance(point), 0.f, 0.01f));
	SMF_ASSERT(CMath::equalsAbs(signedDistance(point + normal_), 1.f, 0.01f));
#endif
}

SMF_INLINE_FORCED float CPlane::operator[]( int index ) const {
	return ( &a )[ index ];
}

SMF_INLINE_FORCED float& CPlane::operator[]( int index ) {
	return ( &a )[ index ];
}

SMF_INLINE_FORCED CPlane CPlane::operator-() const {
	return CPlane( -a, -b, -c, -d );
}

SMF_INLINE_FORCED CPlane &CPlane::operator=( const CVec3D &v ) {
	a = v.x;
	b = v.y;
	c = v.z;
	d = 0;
	return *this;
}

SMF_INLINE_FORCED CPlane CPlane::operator+( const CPlane &p ) const {
	return CPlane( a + p.a, b + p.b, c + p.c, d + p.d );
}

SMF_INLINE_FORCED CPlane CPlane::operator-( const CPlane &p ) const {
	return CPlane( a - p.a, b - p.b, c - p.c, d - p.d );
}

SMF_INLINE_FORCED CPlane &CPlane::operator*=( const CMat3D &m ) {
	getNormal() *= m;
	return *this;
}

SMF_INLINE_FORCED bool CPlane::compare( const CPlane &p ) const {
	return ( a == p.a && b == p.b && c == p.c && d == p.d );
}

SMF_INLINE_FORCED bool CPlane::compare( const CPlane &p, const float epsilon ) const {
	if ( CMath::fabs( a - p.a ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( b - p.b ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( c - p.c ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( d - p.d ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED bool CPlane::compare( const CPlane &p, const float normalEps, const float distEps ) const {
	if ( CMath::fabs( d - p.d ) > distEps ) {
		return false;
	}
	if ( !getNormal().compare( p.getNormal(), normalEps ) ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CPlane::operator==( const CPlane &p ) const {
	return compare( p );
}

SMF_INLINE_FORCED bool CPlane::operator!=( const CPlane &p ) const {
	return !compare( p );
}

SMF_INLINE_FORCED void CPlane::toZero() {
	a = b = c = d = 0.0f;
}

SMF_INLINE_FORCED void CPlane::setNormal( const CVec3D &normal ) {
	a = normal.x;
	b = normal.y;
	c = normal.z;
}

SMF_INLINE_FORCED const CVec3D &CPlane::getNormal() const {
	return *reinterpret_cast<const CVec3D *>(&a);
}

SMF_INLINE_FORCED CVec3D &CPlane::getNormal() {
	return *reinterpret_cast<CVec3D *>(&a);
}

SMF_INLINE_FORCED float CPlane::normalize( bool fixDegenerate ) {
	float length = reinterpret_cast<CVec3D *>(&a)->getLenght();
	reinterpret_cast<CVec3D *>(&a)->toNormal();
	if ( fixDegenerate ) {
		fixDegenerateNormal();
	}
	return length;
}

SMF_INLINE_FORCED bool CPlane::fixDegenerateNormal() {
	return getNormal().fixDegenerateNormal();
}

SMF_INLINE_FORCED bool CPlane::fixDegeneracies( float distEpsilon ) {
	bool fixedNormal = fixDegenerateNormal();
	// only fix dist if the normal was degenerate
	if ( fixedNormal ) {
		if ( CMath::fabs( d - CMath::rint( d ) ) < distEpsilon ) {
			d = CMath::rint( d );
		}
	}
	return fixedNormal;
}

SMF_INLINE_FORCED float CPlane::getDist() const {
	return -d;
}

SMF_INLINE_FORCED void CPlane::setDist( const float dist ) {
	d = -dist;
}

SMF_INLINE_FORCED bool CPlane::fromPoints( const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, bool fixDegenerate ) {

	getNormal() = (p1 - p2).cross( p3 - p2 ); //calcula o vetor normal aos pontos dados
	if ( normalize( fixDegenerate ) == 0.0f ) {
		return false;
	}
	d = -( getNormal() * p2 );
	return true;
}

SMF_INLINE_FORCED bool CPlane::fromVecs( const CVec3D &dir1, const CVec3D &dir2, const CVec3D &p, bool fixDegenerate ) {
	getNormal() = dir1.cross( dir2 );
	if ( normalize( fixDegenerate ) == 0.0f ) {
		return false;
	}
	d = -( getNormal() * p );
	return true;
}
SMF_INLINE_FORCED bool CPlane::fromPointNormal( const CVec3D &normal, const CVec3D &p, bool fixDegenerate ) {
	getNormal() = normal;
	if ( normalize( fixDegenerate ) == 0.0f ) {
		return false;
	}
	d = -( getNormal() * p );
	return true;
}
SMF_INLINE_FORCED void CPlane::fitThroughPoint( const CVec3D &p ) {
	d = -( getNormal() * p );
}

SMF_INLINE_FORCED CPlane CPlane::translate( const CVec3D &translation ) const {
	return CPlane( a, b, c, d - translation * getNormal() );
}

SMF_INLINE_FORCED CPlane &CPlane::translateSelf( const CVec3D &translation ) {
	d -= translation * getNormal();
	return *this;
}

SMF_INLINE_FORCED CPlane CPlane::rotate( const CVec3D &origin, const CMat3D &axis ) const {
	CPlane p;
	p.getNormal() = getNormal() * axis;
	p.d = d + origin * getNormal() - origin * p.getNormal();
	return p;
}

SMF_INLINE_FORCED CPlane &CPlane::rotateSelf( const CVec3D &origin, const CMat3D &axis ) {
	d += origin * getNormal();
	getNormal() *= axis;
	d -= origin * getNormal();
	return *this;
}

SMF_INLINE_FORCED float CPlane::getDistance( const CVec3D &v ) const {
	return a * v.x + b * v.y + c * v.z + d;
}
SMF_INLINE_FORCED float CPlane::getDistance( const CVec4D &v ) const {
	return a * v.x + b * v.y + c * v.z + d*v.w;
}
SMF_INLINE_FORCED float CPlane::getDistanceNormal( const CVec3D &v ) const{
	return a * v.x + b * v.y + c * v.z;
}
SMF_INLINE_FORCED int CPlane::side( const CVec3D &v, const float epsilon ) const {
	float dist = getDistance( v );
	if ( dist > epsilon ) {
		return PLANESIDE_FRONT;
	}
	else if ( dist < -epsilon ) {
		return PLANESIDE_BACK;
	}
	else {
		return PLANESIDE_ON;
	}
}

SMF_INLINE_FORCED bool CPlane::lineIntersection( const CVec3D &start, const CVec3D &end ) const {
	float d1, d2, fraction;

	d1 = getNormal() * start + d;
	d2 = getNormal() * end + d;
	if ( d1 == d2 ) {
		return false;
	}
	if ( d1 > 0.0f && d2 > 0.0f ) {
		return false;
	}
	if ( d1 < 0.0f && d2 < 0.0f ) {
		return false;
	}
	fraction = ( d1 / ( d1 - d2 ) );
	return ( fraction >= 0.0f && fraction <= 1.0f );
}

SMF_INLINE_FORCED bool CPlane::rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale ) const {
	float d1, d2;

	d1 = getNormal() * start + d;
	d2 = getNormal() * dir;
	if ( d2 == 0.0f ) {
		return false;
	}
	scale = -( d1 / d2 );
	return true;
}

SMF_INLINE_FORCED int CPlane::getDimension() const {
	return 4;
}
/*
é possível usar este método de conversão devido à forma como o compilador
armazena as classes na memória, visto que tem a mesma dimensão. veja esses papers:
http://www.openrce.org/articles/full_view/23
http://www.openrce.org/articles/files/jangrayhood.pdf
*/
SMF_INLINE_FORCED const CVec4D &CPlane::ToVec4() const {
	return *reinterpret_cast<const CVec4D *>(&a);
}

SMF_INLINE_FORCED CVec4D &CPlane::ToVec4() {
	return *reinterpret_cast<CVec4D *>(&a);
}

SMF_INLINE_FORCED const float *CPlane::toFloatPtr() const {
	return reinterpret_cast<const float *>(&a);
}

SMF_INLINE_FORCED float *CPlane::toFloatPtr() {
	return reinterpret_cast<float *>(&a);
}

} //end MATH
} //end SMF
#endif /* !__MATH_PLANE_H__ */
