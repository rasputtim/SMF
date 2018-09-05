

#ifndef __SMF_SPHERE_H__
#define __SMF_SPHERE_H__

#include "../SMF_Config.h"
#include "SMF_GeoDefs.h"
#include "SMF_Plane.h"

#include "../math/SMF_Vector.h"

namespace SMF{
namespace MATH {
class CMat3D;

}
using namespace MATH;
namespace GEO {
class CPlane;

/**
 * \class CSphere
 * \image html pics/sphere.png
 * \ingroup SMF_Geometric
 *
 * \brief Implementa uma Esfera de 3  dimensões
 *
 * \brief A 3D CSphere
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CSphere {
public:

	/** 
	\brief Usando este construtor os pontos min e max ficam indefinidos.
	       Você deve se lembrar de definí-los antes de usar a classe.
	\brief This means that the values of the members minPoint and maxPoint are undefined after creating a new CAABBox using this
        default constructor. Remember to assign to them before use.
    \see origin,radius 
	*/
	CSphere();
	
	explicit CSphere( const CVec3D &point );
	/** 
	\brief Cria um objeto CSphere com a dada posição (centro) e raio.

	\brief Constructs a sphere with a given position and radius.
	\param radius O raio da esfera
	\note radius A value > 0 constructs a sphere with positive volume. A value of <= 0 is valid, and constructs a degenerate sphere.
	\param centerpoint ponto central da esfera
	\see origin, radius, isFinite(), isDegenerate() 
	*/
	
	explicit CSphere( const CVec3D &centerpoint, const float radius );

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	/// returns tranlated sphere
	CSphere		operator+( const CVec3D &t ) const;				
	/// translate the sphere
	CSphere &		operator+=( const CVec3D &t );					
	CSphere		operator+( const CSphere &s ) const;
	CSphere &		operator+=( const CSphere &s );
	/// exact compare, no epsilon
	bool			compare( const CSphere &a ) const;							
	/// compare with epsilon
	bool			compare( const CSphere &a, const float epsilon ) const;	
	/// exact compare, no epsilon
	bool			operator==(	const CSphere &a ) const;						
	/// exact compare, no epsilon
	bool			operator!=(	const CSphere &a ) const;						
	/// inside out sphere
	void			clear();									
	/// single point at origin
	void			toZero();									
	/// set origin of sphere
	void			setOrigin( const CVec3D &o );					
	/// set square radius
	void			SetRadius( const float r );						
	/// returns origin of sphere
	const CVec3D &	getOrigin() const;						
	/// returns sphere radius
	float			getRadius() const;						
	/// returns true if sphere is inside out
	bool			isCleared() const;						
	/// add the point, returns true if the sphere expanded
	bool			addPoint( const CVec3D &p );					
	/// add the sphere, returns true if the sphere expanded
	bool			AddSphere( const CSphere &s );					
	/// return bounds expanded in all directions with the given value
	CSphere			expand( const float d ) const;					
	/// expand bounds in all directions with the given value
	CSphere &		expandSelf( const float d );					
	CSphere			translate( const CVec3D &translation ) const;
	CSphere &		translateSelf( const CVec3D &translation );

	float			planeDistance( const CPlane &plane ) const;
	int				planeSide( const CPlane &plane, const float epsilon = ON_EPSILON ) const;
	/// includes touching
	bool			containsPoint( const CVec3D &p ) const;		
	/// includes touching
	bool			intersectsSphere( const CSphere &s ) const;	
	bool			lineIntersection( const CVec3D &start, const CVec3D &end ) const;
	/**
	\brief  Computes the intersection of a line and a sphere.
	 This function solves the points of intersection between a line and a sphere.
		A line intersects sphere at 0, 1 or 2 points. When only one point of intersection is reported,
		the given line is tangent to the sphere.
		\param linePos The source position of the line.
		\param lineDir The direction of the line. This vector must be normalized in advance.
		\param sphereCenter The center position of the sphere to test.
		\param sphereRadius The radius of the sphere, which must be >= 0.
		\param t1 [out] This parameter receives the parametric distance along the line to the first point of intersection.
			If the sphere and the line do not intersect, this variable is not written to. To receive the actual
			world space point corresponding to this point of intersection, compute the vector 'linePos + t1 * lineDir'.
		\param t2 [out] This parameter receives the parametric distance along the line to the second point of intersection.
			If the sphere and the line do not intersect, this variable is not written to. If the line is tangent to this
			sphere (one point of intersection), this variable will be set equal to t1, so that the line segment
			[t1, t2] always forms a line segment completely enclosed inside the sphere. To receive the actual world space
			point corresponding to this point of intersection, compute the vector 'linePos + t2 * lineDir'.
		\return The number of intersection points: 0, 1 or 2. In case of 0, the line and sphere do not intersect. In
			case of 1, the line is tangent to the sphere. If the value of 2 is returned, the line properly intersects the
			sphere.
		\note The outputted variables t1 and t2 always satisfy t1 < t2. This allows distinguishing between the "enter"
			and "exit" positions of the line, if the line is interpreted more like a ray starting at linePos, and extending
			towards lineDir. */
	static int intersectLine(const CVec3D &linePos, const CVec3D &lineDir, const CVec3D &sphereCenter,
	                         float sphereRadius, float &t1, float &t2);

	
	/// intersection points are (start + dir * scale1) and (start + dir * scale2)
	bool			rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale1, float &scale2 ) const;

					/// Tight sphere for a point set.
	void			fromPoints( const CVec3D *points, const int numPoints );
					/// Most tight sphere for a translation.
	void			fromPointTranslation( const CVec3D &point, const CVec3D &translation );
	void			fromSphereTranslation( const CSphere &sphere, const CVec3D &start, const CVec3D &translation );
					/// Most tight sphere for a rotation.
//s	void			fromPointRotation( const CVec3D &point, const CRotation &rotation );
//s	void			FromSphereRotation( const CSphere &sphere, const CVec3D &start, const CRotation &rotation );

//outros
	
	/** 
	\if pt_br
	\brief Testa se a esfera é \b finita \b
	\note A esfera é \b finita \b se seus membros origin e radius não possuem NaNs (not a number http://en.wikipedia.org/wiki/NaN) ou números infinitos +/-. \n
	\elseif us_en
	\brief Tests if this CSphere is finite.
	\note A sphere is \b finite \b if its members origin and radius do not contain floating-point NaNs or +/-infs
		in them.
	\return True if the members pos and r both have finite floating-point values.
	\endif
	\see origin, radius, isDegenerate(), isFinite(), isInf(), IsNan(), isFinite(), CMath::INFINITE_FLOAT, CMath::NEG_INFINITE_FLOAT, CMath::NAN_FLOAT, CVec3D::nan, CVec3D::infinity. */
	bool isFinite() const;

	/** 
	\if pt_br
	\brief Retorna verdaddeiro se a esfera é \b degenerada \b.
	\note Um esfera é degenerada se ela não for finita, ou se seu raio for menor ou igual a zero \n
	\elseif us_en
	\brief Returns true if this CSphere is \b degenerate \b.
	\note A sphere is \b degenerate \b if it is not finite, or if the radius of the sphere is less or equal to 0. \p
	\endif
	\see origin, radius, isFinite() 
	**/
	bool isDegenerate() const;

	
	/** 
	\if pt_br
	\brief Especifica uma constante CSphere com valor (Vector::zero, 0.0f).
	\warning Devido a ordem de inicialização de veriáveis estáticas pelo compilador ser indefinida em C++,
	       não use esta variável para inicializar outras variáveis estáticas em outras unidades de compilação.

	\elseif us_en
	\brief Specifies a compile-time constant CSphere with value (Vector::zero, 0.0f). 
	\warning Due to static data initialization order being undefined in C++, do NOT use this
			member to initialize other static data in other compilation units! 
	\endif
	**/
	static const CSphere zero;
	//==============Distances======================

	/**
	\brief Returns the distance between this sphere and the given object.
	 This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
	\see contains(), intersects(), closestPoint().
	\todo add Sphere::distance(Polygon/Circle/Disc/Frustum/Polyhedron). */
	float distance(const CVec3D &point) const;
	float distance(const CSphere &sphere) const;
	float distance(const CAABBox &aabb) const;
	float distance(const COBBox &obb) const;
	float distance(const CPlane &plane) const;
	float distance(const CTriangle &triangle) const;
	float distance(const CRay &ray) const;
	float distance(const CLine &line) const;
	float distance(const CLineSegment &lineSegment) const;

		/// Projects this Sphere onto the given 1D axis direction vector.
	/** This function collapses this Sphere onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
		\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
		\param outMin [out] Returns the minimum extent of this object along the projection axis.
		\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;


		/// Tests whether this sphere and the given object intersect.
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this sphere,
		or this sphere is contained inside a polyhedron, etc.)
		\param intersectionPoint [out] If specified, receives the actual point of intersection. This pointer may be null.
		\param intersectionNormal [out] If specified, receives the sphere normal at the point of intersection. This pointer may be null.
		\param d [out] If specified, receives the distance along the Line/LineSegment/Ray to the intersection. This pointer may be null.
		\return In the case of Ray/Line/LineSegment intersection tests, the number of intersections is returned (0, 1 or 2).
			For other functions, true is returned if an intersection occurs or one of the objects is contained inside the other, and false otherwise.
		\see contains(), distance(), closestPoint(), LineSegment::getPoint(). */
	int intersects(const CLineSegment &lineSegment, CVec3D *intersectionPoint = 0, CVec3D *intersectionNormal = 0, float *d = 0, float *d2 = 0) const;
	int intersects(const CLine &line, CVec3D *intersectionPoint = 0, CVec3D *intersectionNormal = 0, float *d = 0, float *d2 = 0) const;
	int intersects(const CRay &ray, CVec3D *intersectionPoint = 0, CVec3D *intersectionNormal = 0, float *d = 0, float *d2 = 0) const;
	bool intersects(const CPlane &plane) const;
	bool intersects(const CAABBox &aabb, CVec3D *closestPointOnAABB = 0) const;
	bool intersects(const COBBox &obb, CVec3D *closestPointOnOBB = 0) const;
	bool intersects(const CTriangle &triangle, CVec3D *closestPointOnTriangle = 0) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	bool intersects(const CSphere &sphere) const;
	/// Returns the smallest AABB that encloses this sphere.
	/** The returned AABB is a cube, with a center position coincident with this sphere, and a side length of 2*r.
		\see maximalContainedAABB(). */
	CAABBox minimalEnclosingAABB() const;

	/// Returns the largest AABB that fits inside this sphere.
	/** The returned AABB is a cube, with a center position coincident with this sphere, and a side length of 2*r/sqrt(3).
		\see minimalEnclosingAABB(). */
	CAABBox maximalContainedAABB() const;

	/// Sets pos = (0,0,0) and r = -inf.
	/** After a call to this function, both isFinite() and isDegenerate() will return true.
		\see isFinite(), isDegenerate(). */
	void toNegativeInfinity();

private:
	/**
	\if pt_br
	\brief O centro da esfera
	\elseif us_en
	\brief The center of the CSphere
	\endif
	**/
	CVec3D		origin;
	
	/**
	\if pt_br
	\brief o raio da esfera
	\elseif us_en
	\brief The CSphere radius
	\endif
	**/
	float			radius;
};



SMF_INLINE CSphere::CSphere() {
}

SMF_INLINE CSphere::CSphere( const CVec3D &point ) {
	origin = point;
	radius = 0.0f;
}

SMF_INLINE CSphere::CSphere( const CVec3D &point, const float r ) {
	origin = point;
	radius = r;
}

SMF_INLINE float CSphere::operator[]( const int index ) const {
	return ((float *) &origin)[index];
}

SMF_INLINE float &CSphere::operator[]( const int index ) {
	return ((float *) &origin)[index];
}

SMF_INLINE CSphere CSphere::operator+( const CVec3D &t ) const {
	return CSphere( origin + t, radius );
}

SMF_INLINE CSphere &CSphere::operator+=( const CVec3D &t ) {
	origin += t;
	return *this;
}

SMF_INLINE bool CSphere::compare( const CSphere &a ) const {
	return ( origin.compare( a.origin ) && radius == a.radius );
}

SMF_INLINE bool CSphere::compare( const CSphere &a, const float epsilon ) const {
	return ( origin.compare( a.origin, epsilon ) && CMath::fabs( radius - a.radius ) <= epsilon );
}

SMF_INLINE bool CSphere::operator==( const CSphere &a ) const {
	return compare( a );
}

SMF_INLINE bool CSphere::operator!=( const CSphere &a ) const {
	return !compare( a );
}

SMF_INLINE void CSphere::clear() {
	origin.toZero();
	radius = -1.0f;
}

SMF_INLINE void CSphere::toZero() {
	origin.toZero();
	radius = 0.0f;
}

SMF_INLINE void CSphere::setOrigin( const CVec3D &o ) {
	origin = o;
}

SMF_INLINE void CSphere::SetRadius( const float r ) {
	radius = r;
}

SMF_INLINE const CVec3D &CSphere::getOrigin() const {
	return origin;
}

SMF_INLINE float CSphere::getRadius() const {
	return radius;
}

SMF_INLINE bool CSphere::isCleared() const {
	return ( radius < 0.0f );
}

SMF_INLINE bool CSphere::addPoint( const CVec3D &p ) {
	if ( radius < 0.0f ) {
		origin = p;
		radius = 0.0f;
		return true;
	}
	else {
		float r = ( p - origin ).getLengthSqr();
		if ( r > radius * radius ) {
			r = CMath::sqrt( r );
			origin += ( p - origin ) * 0.5f * (1.0f - radius / r );
			radius += 0.5f * ( r - radius );
			return true;
		}
		return false;
	}
}

SMF_INLINE bool CSphere::AddSphere( const CSphere &s ) {
	if ( radius < 0.0f ) {
		origin = s.origin;
		radius = s.radius;
		return true;
	}
	else {
		float r = ( s.origin - origin ).getLengthSqr();
		if ( r > ( radius + s.radius ) * ( radius + s.radius ) ) {
			r = CMath::sqrt( r );
			origin += ( s.origin - origin ) * 0.5f * (1.0f - radius / ( r + s.radius ) );
			radius += 0.5f * ( ( r + s.radius ) - radius );
			return true;
		}
		return false;
	}
}

SMF_INLINE CSphere CSphere::expand( const float d ) const {
	return CSphere( origin, radius + d );
}

SMF_INLINE CSphere &CSphere::expandSelf( const float d ) {
	radius += d;
	return *this;
}

SMF_INLINE CSphere CSphere::translate( const CVec3D &translation ) const {
	return CSphere( origin + translation, radius );
}

SMF_INLINE CSphere &CSphere::translateSelf( const CVec3D &translation ) {
	origin += translation;
	return *this;
}

SMF_INLINE bool CSphere::containsPoint( const CVec3D &p ) const {
	if ( ( p - origin ).getLengthSqr() > radius * radius ) {
		return false;
	}
	return true;
}

SMF_INLINE bool CSphere::intersectsSphere( const CSphere &s ) const {
	float r = s.radius + radius;
	if ( ( s.origin - origin ).getLengthSqr() > r * r ) {
		return false;
	}
	return true;
}

SMF_INLINE void CSphere::fromPointTranslation( const CVec3D &point, const CVec3D &translation ) {
	origin = point + 0.5f * translation;
	radius = CMath::sqrt( 0.5f * translation.getLengthSqr() );
}

SMF_INLINE void CSphere::fromSphereTranslation( const CSphere &sphere, const CVec3D &start, const CVec3D &translation ) {
	origin = start + sphere.origin + 0.5f * translation;
	radius = CMath::sqrt( 0.5f * translation.getLengthSqr() ) + sphere.radius;
}
#if 0
SMF_INLINE void CSphere::fromPointRotation( const CVec3D &point, const CRotation &rotation ) {
	CVec3D end = rotation * point;
	origin = ( point + end ) * 0.5f;
	radius = CMath::sqrt( 0.5f * ( end - point ).getLengthSqr() );
}

SMF_INLINE void CSphere::FromSphereRotation( const CSphere &sphere, const CVec3D &start, const CRotation &rotation ) {
	CVec3D end = rotation * sphere.origin;
	origin = start + ( sphere.origin + end ) * 0.5f;
	radius = CMath::sqrt( 0.5f * ( end - sphere.origin ).getLengthSqr() ) + sphere.radius;
}
#endif
SMF_INLINE void CSphere::projectToAxis( const CVec3D &dir, float &min, float &max ) const {
	float d;
	d = dir * origin;
	min = d - radius;
	max = d + radius;
}


} //end GEO
} //end SMF
#endif /* !__SMF_SPHERE_H__ */
