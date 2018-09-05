#ifndef __SMF_RAY_2D_
#define __SMF_RAY_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{


#ifdef MATH_OGRE_INTEROP
#include <OgreRay.h>
#endif

class CLine2D;
class CLineSegment2D;
class CPolygon2D;
class CPolyhedron2D;
class CCircle2D;
class CTriangle2D;
/**
 * \class CRay2D
 *
 * \ingroup SMF_Geometric
 *
 * \image html pics\line.png
 * \if pt_br
 * \brief   Representa um Raio em espaço cartesiano 2D
    \note  Um raio é definido por um ponto (origem) e um vetor (direção).
		   \n Os parâmetros são acondicionados num vetor RAY = [x0 y0 dx dy];
		   \n O raio comtém todos os pontos (x,y) de tal forma que \f$ x = x0 + t*dx y = y0 + t*dy \f$;
		   \n para todos t>0 . Ao contrário a uma reta, os pontos localizados antes da origem não pertemcem ao raio.
		   \n contudo, todos os raios e retas possuem a mesma representação; possuem  algumas funções em comum, como (like transformLine).
 * \elseif us_en
 * \brief 	A ray in 2D space is a line that starts from an origin point and extends to infinity in one direction.
 *  \note A ray is defined by a point (its origin), and a vector (its direction). 
         \n The different parameters are bundled into a row vector: RAY = [x0 y0 dx dy]; 
		 \n The ray contains all the points (x,y) such that: x = x0 + t*dx y = y0 + t*dy; 
		 \n for all t>0 Contrary to a (straight) line, the points located before the origin do not belong to the ray. 
		 \n However, as rays and lines have the same representation, some functions working on lines are also working on rays (like transformLine).

 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CRay2D
{
public:
	/// Specifies the origin of this ray.
	CPoint2D pos;

	/** \brief The normalized direction vector of this ray. [similarOverload: pos]
	 \note For proper functionality, this direction vector needs to always be normalized. If you set to this
		member manually, remember to make sure you only assign normalized direction vectors. */
	CVec2D dir;

	/** \brief The default constructor does not initialize any members of this class.
	 This means that the values of the members pos and dir are undefined after creating a new CRay2D using this
		default constructor. Remember to assign to them before use.
	\see pos, dir. */
	CRay2D() {}

	/** 
	\brief Constructs a new ray by explicitly specifying the member variables.
	\param pos The origin position of the ray.
	\param dir The direction of the ray. This vector must be normalized, this function will not normalize
			the vector for you (for performance reasons).
	\see pos, dir. */
	CRay2D(const CPoint2D &pos, const CVec2D &dir);
	/** 
	\brief Constructs a new ray by explicitly specifying the member variables.
	\param pos The origin position of the ray.
	\param p2 an arbitrary point where the ray must pass through
	\see pos, dir. */
	CRay2D(const CPoint2D &pos, const CPoint2D &p2);


	/// Converts a CLineSegment2D to a CRay2D.
	/** This constructor sets pos = lineSegment.point1 (default), or point2, and dir = (lineSegment.b - lineSegment.a).normalized().
		\param poitToUse point to use from the line segment to be the begin of this ray: 0-> point1, 1-> point2
		\see class CLineSegment2D, toLineSegment(). */
	explicit CRay2D(const CLineSegment2D &lineSegment, int pointToUse=0);
	/**
	\brief return the point origin of the ray
	**/
	CPoint2D origin()const { return pos;}
	bool isFinite() const;

	/** \brief Gets a point along the ray at the given distance.
	 Use this function to convert a 1D parametric point along the CRay2D to a 3D point in the linear space.
	\param distance The point to compute. getPoint(0) will return pos. getPoint(t) will return a point
			at distance |t| from pos. Passing in negative values is allowed, but in that case, the
			returned point does not actually lie on this CRay2D.
	\return pos + distance * dir.
	\see pos, dir. */
	CVec2D getPoint(float distance) const;


	/**
      \brief check whether p is on the direction of this Ray
      \param point considered point
      \param thr threshold angle buffer
      \return true or false
    **/
    bool inRightDir( const CPoint2D & point, const float & thr = 10.0 ) const
      {
          CVec2D vec = (point - origin());
		  float ang1 = CMath::atan( vec.y, vec.x );
		  float ang2 = CMath::atan( dir.y, dir.x );
		  return CMath::fabs( ang1 - ang2 ) < thr;
      }
	/** 
	\brief Translates this CRay2D in world space.
	\param offset The amount of displacement to apply to this CRay2D, in world space coordinates.
	\see transform(). */
	void translate(const CVec2D &offset);

	/** 
	\brief Applies a transformation to this CRay2D, in-place.
	\See translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	//void transform(const CMat3D &transform);
	//void transform(const CMatJoint3x4 &transform);
	//void transform(const CMat4D &transform);
	//void transform(const CQuaternion &transform);

	/** 
	\brief Tests if the given object is fully contained on this ray.
	\param distanceThreshold The magnitude of the epsilon test threshold to use. Since a CRay2D
		is a 1D object in a 3D space, an epsilon threshold is used to allow errors caused by floating-point
		inaccuracies.
	\return True if this ray contains the given object, up to the given distance threshold.
	\see class CLineSegment2D, distance(), closestPoint(), intersects(). */
	bool contains(const CPoint2D &point, float distanceThreshold = 1e-3f) const;
	bool contains(const CLineSegment2D &lineSegment, float distanceThreshold = 1e-3f) const;

	/** 
	\brief Tests if two rays are equal.
	\return True if this and the given CRay2D represent the same set of points, up to the given epsilon. */
	bool compare(const CRay2D &otherRay, float epsilon =CMath::EPSILON_SuperLow) const;

	/** 
	\brief Computes the distance between this ray and the given object.
	This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
	\param d [out] If specified, receives the parametric distance along this ray that
			specifies the closest point on this ray to the given object. The value returned here can be negative.
			This pointer may be null.
	\see contains(), intersects(), closestPoint(), getPoint(). */
	float distance(const CVec2D &point, float *d) const;
	float distance(const CVec2D &point) const;

	/** 
	\param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this ray. The value returned here can be negative. This pointer may
		be null. */
//	float distance(const CRay2D &other, float *d, float *d2 = 0) const;
//	float distance(const CRay2D &other) const;
	float distance(const CLine2D &other, float *d, float *d2 = 0) const;
//	float distance(const CLine2D &other) const;
//	float distance(const CLineSegment2D &other, float *d, float *d2 = 0) const;
//	float distance(const CLineSegment2D &other) const;
	float distance(const CCircle2D &sphere) const;
//	float distance(const Capsule &capsule) const;

	/** 
	\brief  Computes the closest point on this ray to the given object.
	If the other object intersects this ray, this function will return an arbitrary point inside
		the region of intersection.
	\param d [out] If specified, receives the parametric distance along this ray that
			specifies the closest point on this ray to the given object. The value returned here can be negative.
			This pointer may be null.
	\see contains(), distance(), intersects(), getPoint(). */
	CVec2D closestPoint(const CVec2D &targetPoint, float *d = 0) const;
	/** 
	\param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this ray. The value returned here can be negative. This pointer may
		be null. */
	//CVec2D closestPoint(const CRay2D &other, float *d = 0, float *d2 = 0) const;
	CVec2D closestPoint(const CLine2D &other, float *d = 0, float *d2 = 0) const;
	CVec2D closestPoint(const CLineSegment2D &other, float *d = 0, float *d2 = 0) const;
	
	CPoint2D closestPoint(CLineSegment2D const &l1)const;
    CPoint2D closestPoint(CPoint2D const  &p)const;
    CPoint2D closestPoint(CRay2D const &r1)const;

  /**
  \brief check intersection between one ray and other line,segment or ray
  \param [out] intersection if there in no intersection return NAN_FLOAT else return the intersection point
  \param [out] intType the type of intersection
  \return true if there is any kind of intersection
  \see CLineSegment2D::intersects_m2
  **/
  static bool  intersects_m2(const CRay2D& this_line, const CLine2D& other_line, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
  static bool  intersects_m2(const CRay2D& this_line, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
  static bool  intersects_m2(const CRay2D& this_line, const CRay2D& other_ray, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);


	/// Tests whether this ray and the given object intersect.	

    /*!
      \brief calculate the intersection points with circle
      \param circle considered
      \param [out] ICnt number of intersections
      \param [out] I1 pointer to the 1st solution variable
      \param [out] I2 pointer to the 2nd solution variable
      \return true if there is an intersection
     */
    bool intersection( const CCircle2D &circle,int *ICnt=NULL,CPoint2D * I1=NULL, CPoint2D * I2=NULL ) const;
	/*!
      \brief calculate the intersection points with triangle
      \param triangle considered 
      \param [out] ICnt number of intersections
      \param [out] I1 pointer to the 1st solution variable
      \param [out] I2 pointer to the 2nd solution variable
      \return true if there is an intersection
     */
    bool intersection(const CTriangle2D &triangle,  int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const;
    bool intersection(const CLine2D &line,  int *ICnt=NULL, CPoint2D *I1=NULL) const;
    bool intersection(const CLineSegment2D &segment,  int *ICnt=NULL, CPoint2D *I1=NULL) const;
    bool intersection(const CRay2D &ray,  int *ICnt=NULL, CPoint2D *I1=NULL) const;


	/** \param dNear [out] If specified, receives the distance along this ray to where the ray enters
		the bounding box.
		\param dFar [out] If specified, receives the distance along this ray to where the ray exits
		the bounding box. */
	bool intersects(const CAABBox2D &aabb, float &dNear, float &dFar) const;
	bool intersects(const CAABBox2D &aabb) const;
//	bool intersects(const COBBox2D &obb, float &dNear, float &dFar) const;
//	bool intersects(const COBBox2D &obb) const;
//	bool intersects(const CPolygon2D &polygon) const;
//	bool intersects(const CPolyhedron2D &polyhedron) const;
 /** \brief Check if the point intersects the line segment
  \see http://lucidarme.me/?p=1952
  **/
  bool intersects(CPoint2D const &p)const;
  /**
  \brief check if this Ray intersects the other Line
  **/
  bool intersects(CLine2D const &p)const;
  /**
  \brief check if this Ray intersects the other line  segment
  **/
  bool intersects(CLineSegment2D const &p)const;
  /**
  \brief check if this Ray intersects the other ray
  **/
  bool intersects(CRay2D const &p)const;

  bool intersects(const CCircle2D &s) const;

	bool intersects(const CTriangle2D &s) const;


  /*!
 * \brief isPointOnSegment check if a point is inside the current segment
 * \param point coordinates of the point to test
 * \param [out] intersecResult return the type of intersection found
 *   NOT_INTERESECTING if the point doesn't lay with the segment
 *   INTERESECTING_EXTREMITY_P1 if the point is merged with P1
 *   INTERESECTING_EXTREMITY_P2 if the point is merged with P2
 *   INTERESECTING if the point belongs to the segment (extremity no included)
 *  \return true if it has any intersection
 */
  bool isPointOnRay(CPoint2D const &p, int *intersecResult=NULL)const;
  



	/** 
	\brief  Tests if this ray intersects the given disc.
	\todo This signature will be moved to bool intersects(const Disc &disc) const;*/
	bool intersectsDisc(const CCircle2D &disc) const;

	/** 
	\brief  Converts this CRay2D to a CLine2D.
	 The pos and dir members of the returned CLine2D will be equal to this CRay2D. The only difference is
		that a CLine2D extends to infinity in two directions, whereas the CRay2D spans only in the positive
		direction.
	\see dir, CRay2D::CRay2D, class CLine2D, toLineSegment(). */
	CLine2D toLine2D() const;


	/**
	\if pt_br
	\brief Converte este raio num segmento 2D
	\note o ponto ininical do segmento é o centro do raio e o ponot final do segmento está na distância dada
	\param dist a distancia do centro aé o ponto final do segmento a ser criado
	\elseif us_en
	\brief Converts this CRay2D into a LineSegment2D
	\param d Specifies the position of the other endpoint along this CRay2D. This parameter may be negative,
		in which case the returned CLineSegment2D does not lie inside this CRay2D.
		\return A CLineSegment2D with point a at pos, and point b at pos + d * dir.
		\see pos, dir, CRay2D::CRay2D, class CLineSegment2D, toLine().
	endif
	**/
	CLineSegment2D toLineSegment2D(float d) const;

	/**
	\brief Converts this CRay2D to a CLineSegment2D.
	 \param dStart Specifies the position of the first endpoint along this CRay2D. This parameter may be negative,
		in which case the starting point lies outside this CRay2D to the opposite direction of the CRay2D.
		\param dEnd Specifies the position of the second endpoint along this CRay2D. This parameter may also be negative.
		\return A CLineSegment2D with point a at pos + dStart * dir, and point b at pos + dEnd * dir.
		\see pos, dir, CRay2D::CRay2D, class CLineSegment2D, toLine(). */
	CLineSegment2D toLineSegment2D(float dStart, float dEnd) const;

	/** 
	\brief  Projects this CRay2D onto the given 1D axis direction vector.
	 This function collapses this CRay2D onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
	\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
	\param outMin [out] Returns the minimum extent of this object along the projection axis.
	\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec2D &direction, float &outMin, float &outMax) const;


	
	/** 
	\brief  Returns a human-readable representation of this CRay2D.
	 The returned string specifies the position and direction of this CRay2D. */
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif

#ifdef MATH_OGRE_INTEROP
	CRay2D(const Ogre::CRay2D &other) { pos = other.getOrigin(); dir = other.getDirection(); }
	operator Ogre::CRay2D() const { return Ogre::CRay2D(pos, dir); }
#endif

};


#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CRay2D)
Q_DECLARE_METATYPE(CRay2D*)
#endif

namespace _2D{
std::ostream &operator <<(std::ostream &o, const CRay2D &ray);
} //end _2D


} //end GEO
}  //end SMF

#endif // __SMF_RAY_2D_
