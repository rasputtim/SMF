


#ifndef __SMF_AABB2D_H__
#define __SMF_AABB2D_H__

#include "../SMF_Config.h"
#include "../math/SMF_Vector.h"
#include "SMF_2DPoint.h"

namespace SMF{
namespace MATH{
class CMat2D;
}
namespace GEO{
using namespace MATH;
class CPoint2D;
class CLine2D;
class CTriangle2D;
class CCircle2D;
class CPolygon2D;
/**
 * \class CAABBox2D
 *
 * \ingroup SMF_Geometric
 * \image html pics/aabb.png
 * \if pt_br 
 * \brief (AABB) Um retângulo alinhado com os eixos \a x e \y em espaço 2D
 * \elseif us_en
 * \brief 	(AABB 2D) 2D Axis Aligned Bounding Box
 *  Represents a cartesian 2D box
 *  A 2d box with sides aligned with the \a x and \a y axes.
 *  Also supports operations required of a bounding box for geometric region
 *  tests.
 *
 *  A box can be empty; this is what the default constructor creates, or what
 *  is left after applying the empty() method.  Use the add() methods to enlarge
 *  a box, and use the contains() methods to check for inclusion of a point or
 *  an other box.
 *
 *  To make the convex union of two boxes, use box1.add(box2).
 *  \verbatim
 *                                  MaxPosition
 *                    O-------------O
 *                    |             |
 *                    |             |
 *                    |  centroid   |
 *                    |      o      |
 *                    |             |
 *        Y           |             |
 *        |           |             |
 *        |           O-------------O
 *        |       MinPosition
 *        O------X
 * \endverbatim
 * If you are using a CAABBox2D<int> to indicate a window on an image, do not forget
 * that your axes will be flipped. You could think of the window as follows.
 *  \verbatim
 *        O------X
 *        |       MinPosition
 *        |             O-------------O
 *        |             |             |
 *        Y             |             |
 *                      |  centroid   |
 *                      |      o      |
 *                      |             |
 *                      |             |
 *                      |             |
 *                      O-------------O
 *                               MaxPosition
 * \endverbatim
 * 
 *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */

class SMF_API CAABBox2D
{
 public:

  /** \brief Default constructor (creates empty box)*/
  CAABBox2D();

  /** \brief Construct using two corner points*/
  CAABBox2D(float const corner1[2],
             float const corner2[2]);

  /** \brief Construct using two corner points*/
  CAABBox2D(CPoint2D const& corner1,
             CPoint2D const& corner2);

  /** \brief Construct using ranges in \a x (first two args) and \a y (last two)*/
  CAABBox2D(float xmin, float xmax, float ymin, float ymax);

  enum point_Type { centre=0, min_pos, max_pos };

  /** \brief Construct a box sized width x height at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 4th argument.*/
  CAABBox2D(float const ref_point[2],
             float width, float height,
             point_Type);

  /** \brief Construct a box sized width x height at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 4th argument.*/
  CAABBox2D(CPoint2D const& ref_point,
             float width, float height,
             point_Type);

  /** \brief Equality test*/
  SMF_INLINE bool operator==(CAABBox2D const& b) const {
    // All empty boxes are equal:
    if (b.isEmpty()) return isEmpty();
    return  minX() == b.minX() && minY() == b.minY()
         && maxX() == b.maxX() && maxY() == b.maxY();
  }

  // Data Access---------------------------------------------------------------

  /** \brief Get width of this box (= \a x dimension)*/
  float width() const;
  /** \brief Get height of this box (= \a y dimension)*/
  float height() const;

  /** \brief Get "volume" (=area) of this box*/
  float area() const { return width()*height(); }
  /** \brief Get "volume" (=area) of this box*/
  float volume() const { return width()*height(); }

  /** \brief Get min \a x*/
  SMF_INLINE float minX() const { return min_pos_[0]; }
  /** \brief Get min \a y*/
  SMF_INLINE float minY() const { return min_pos_[1]; }
  /** \brief Get max \a x*/
  SMF_INLINE float maxX() const { return max_pos_[0]; }
  /** \brief Get max \a y*/
  SMF_INLINE float maxY() const { return max_pos_[1]; }

  /** \brief Get the centroid point*/
  CPoint2D centroid() const;
  /** \brief Get \a x component of centroid*/
  float getCentroidX() const;
  /** \brief Get \a y component of centroid*/
  float getCentroidY() const;

  /** \brief Return lower left corner of box*/
  CPoint2D minPoint() const;

  /** \brief Return upper right corner of box*/
  CPoint2D maxPoint() const;

  // Data Control--------------------------------------------------------------

  /** \brief Return true if this box is empty*/
  SMF_INLINE bool isEmpty() const {
    return minX() > maxX() || minY() > maxY();
  }

  /** \brief add a point to this box.
  // Do this by possibly enlarging the box so that the point just falls within the box.
  // Adding a point to an empty box makes it a size zero box only containing p.*/
  void add(CPoint2D const& p);

  /** \brief Make the convex union of two boxes.
  // Do this by possibly enlarging this box so that the corner points of the
  // given box just fall within the box.
  // Adding an empty box does not change the current box.*/
  void add(CAABBox2D const& b);

  /** \brief Return true if the point p is inside this box*/
  bool contains(CPoint2D const& p) const;

  /** \brief Return true if the corner points of b are inside this box*/
  bool contains(CAABBox2D const& b) const;

  /** \brief Return true if \a (x,y) inside box, ie \a x_min <= \a x <= \a x_max etc*/
  SMF_INLINE bool contains(float const& x, float const& y) const {
    return x >= minX() && x <= maxX() && y >= minY() && y <= maxY();
  }

/** \brief Return true if the point lies inside the box
// \see  CPoint2D
// \see  CAABBox2D*/

SMF_INLINE bool intersects(CPoint2D const& p) { return contains(p); }

    /**
  /if pt_br
  \brief Retorna Verdadeiro se a linha passada intersepta este AABBox2D. Se sim, computa os pontos de intersecção
  \elseif us_en
  \brief Return true if line intersects box. If so, compute intersection points.
  \endif
  **/
  bool intersects(const CLine2D& line)const;
  bool intersects(const CLineSegment2D &segment)const;
  bool intersects(const CAABBox2D &rhs) const;
  bool intersects(const CTriangle2D &tri) const;
  bool intersects(const CCircle2D &circle)const;

  /** \brief Return the intersection of two boxes (which is itself either a box, or empty)*/
  CAABBox2D intersection(CAABBox2D const&);


/** \brief Return true if line intersects box. If so, compute intersection points.
// \see  CLine2D*/

bool intersection(CLine2D const& line, int *ICnt=NULL, CPoint2D* p0=NULL, CPoint2D* p1=NULL)const ; 

/** Returns true if the box and line regions intersect
\param [out] ICnt the number of intersections of a line segment with a box, up to two are returned in p0 and p1.
// \see  CLineSegment2D
**/
bool intersection(CLineSegment2D const& segment, int *ICnt=NULL, CPoint2D* p0=NULL, CPoint2D* p1=NULL);

/** \brief Return true if the box and polygon regions intersect, regions include boundaries
// \see  CPolygont2D
// \see  CAABBox2D*/
bool intersection(CPolygon2D const& poly);

/** \brief Return the points from the list that lie inside the box
// \see  CPoint2D
// \see  CAABBox2D
*/
std::vector< CPoint2D > intersection(std::vector< CPoint2D > const& p);

  /** \brief Make the box empty*/
  void empty();

  /** \brief set left side of box (other side ordinates unchanged)*/
  SMF_INLINE void setMinX(float m) { min_pos_[0]=m; }
  /** \brief set bottom of box (other side ordinates unchanged)*/
  SMF_INLINE void setMinY(float m) { min_pos_[1]=m; }
  /** \brief set right side (other side ordinates unchanged)*/
  SMF_INLINE void setMaxX(float m) { max_pos_[0]=m; }
  /** \brief set top (other side ordinates unchanged)*/
  SMF_INLINE void setMaxY(float m) { max_pos_[1]=m; }

  /** \brief Move box so centroid lies at cx (width and height unchanged)*/
  void setCentroidX(float cx);
  /** \brief Move box so centroid lies at cy (width and height unchanged)*/
  void setCentroidY(float cy);

  /** \brief Modify width, retaining centroid at current position*/
  void setWidth(float width);
  /** \brief Modify height, retaining centroid at current position*/
  void setHeight(float height);

  /** \brief add to width and height, centroid unchanged.
  // Will move each side by \p expand / 2.*/
  void expandAboutCentroid(float expand);
  /** \brief scale width and height, centroid unchanged.*/
  void scaleAboutCentroid(double s);
  /** \brief scale width and height, keeping scaled position of origin unchanged.*/
  void scaleAboutOrigin(double s);

  /** \brief Modify bottom left. Top right only changed if necessary to avoid empty box*/
  void setMinPosition(float const min_position[2]);
  /** \brief Modify top right. Bottom left only changed if necessary to avoid empty box*/
  void setMaxPosition(float const max_position[2]);
  /** \brief Modify bottom left. Top right only changed if necessary to avoid empty box*/
  void setMinPoint(CPoint2D const& min_pt);
  /** \brief Modify top right. Bottom left only changed if necessary to avoid empty box*/
  void setMaxPoint(CPoint2D const& max_pt);

  /** \brief Move box so centroid lies at c (width, height unchanged)*/
  SMF_INLINE void setCentroid(float const c[2]) { setCentroidX(c[0]); setCentroidY(c[1]); }
  /** \brief Move box so centroid lies at c (width, height unchanged)*/
  SMF_INLINE void setCentroid(CPoint2D const& c) { setCentroidX(c.x()); setCentroidY(c.y()); }



	SMF_INLINE bool hasNegativeVolume() const
	{
		return max_pos_[0] < min_pos_[0] || max_pos_[1] < min_pos_[1];
	}

	SMF_INLINE bool isFinite() const
	{
		return CVec2D(min_pos_[0],min_pos_[1]).isFinite() && CVec2D(max_pos_[0],max_pos_[1]).isFinite() && CVec2D(min_pos_[0],min_pos_[1]).minElement() > -1e5f && CVec2D(max_pos_[0],max_pos_[1]).maxElement() < 1e5f;
	}
	SMF_INLINE 	bool isDegenerate() const
	{
		return min_pos_[0] >= max_pos_[0] || min_pos_[1] >= max_pos_[1];
	}

//	float distanceSq(const CVec2D &pt) const
//	{
//		CVec2D cp = pt.clamp(CVec2D(min_pos_[0],min_pos_[1]), CVec2D(max_pos_[0],max_pos_[1]));
//		return cp.distanceSq(pt);
//	}

	SMF_INLINE void toNegativeInfinity()
	{
		max_pos_[0]=(CMath::INFINITY_FLOAT);
		max_pos_[1]=(CMath::INFINITY_FLOAT);
		min_pos_[0]=(-CMath::INFINITY_FLOAT);
		min_pos_[1]=(-CMath::INFINITY_FLOAT);
	}

	SMF_INLINE void enclose(const CPoint2D &point)
	{
		min_pos_[0] = MIN(min_pos_[0], point.x_);
		min_pos_[1] = MIN(min_pos_[1], point.y_);
		max_pos_[0] = MAX(max_pos_[0], point.x_);
		max_pos_[1] = MAX(max_pos_[1], point.y_);
	}

CPoint2D ClosestPointOnRectangleFromPoint(const CPoint2D &point)const;
  static void test_box_2d();

  // I/O-----------------------------------------------------------------------

  /** \brief Write "<CAABBox2D x0,y0 to x1,y1>" to stream*/
  ostream& print(ostream&) const;

  /** \brief Write "x0 y0 x1 y1(endl)" to stream*/
  ostream& write(ostream&) const;

  /** \brief Read x0,y0,x1,y1 from stream*/
  istream& read(istream&);

  // INTERNALS-----------------------------------------------------------------
 protected:
  // Data Members--------------------------------------------------------------
  float min_pos_[2];
  float max_pos_[2];
};

/** \brief Write box to stream
// \see  CAABBox2D*/

ostream&  operator<<(ostream& s, CAABBox2D const& p);

/** \brief Read box from stream
// \see  CAABBox2D*/

istream&  operator>>(istream& is,  CAABBox2D & p);
#if 0
/** \brief Calculate the bounding box of a sequence of points or boxes.*/
void CAABBox2D_bounds(ITER begin, ITER end, CAABBox2D<T>& bounding_box)
{
  for (; begin != end; ++begin)
    bounding_box.add(*begin);
}

#endif



} // end GEO
} //end SMF

#endif /* !__SMF_BOUNDS_H__ */
