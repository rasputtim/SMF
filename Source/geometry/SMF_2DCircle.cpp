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

#include "SMF_Config.h"
#include "geometry/SMF_2DCircleSegment.h"
#include "geometry/SMF_2DLineSegment.h"
#include "geometry/SMF_2DLine.h"
#include "geometry/SMF_2DConics.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{



float distance(const CPoint2D &p1,const CPoint2D &p2)
{
	return CMath::sqrt((p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y()));
}
 /**
\brief return true if it is possible to get a circle from the points and radius
\param [out] ICnt a quantidade de possiveis circles [0,1,2]
\param [out] center1 the center of the first possible circle
\param [out] center2 the center of the second possible circle
\see http://rosettacode.org/wiki/Circles_of_given_radius_through_two_points
**/
static bool findCircles(CPoint2D p1, CPoint2D p2, float radius, int &ICnt, CPoint2D &center1, CPoint2D &center2)
{
	center1.toInfinite();
	center2.toInfinite();
	ICnt=0;
	float separation = distance(p1,p2), mirrorDistance;

	if(separation == 0.0)
	{

		//printf("\nNo circles can be drawn through (%.4f,%.4f)",p1.x(),p1.y())
	    return false;
	}//else printf("\nInfinitely many circles can be drawn through (%.4f,%.4f)",p1.x(),p1.y());
	if(separation == 2*radius)
	{
		ICnt=1;
		center1.set((p1.x()+p2.x())/2,(p1.y()+p2.y())/2);
		//printf("\nGiven points are opposite ends of a diameter of the circle with center (%.4f,%.4f) and radius %.4f",(p1.x()+p2.x())/2,(p1.y()+p2.y())/2,radius);
	    return true;
	}

	else if(separation > 2*radius)
	{
		//printf("\nGiven points are farther away from each other than a diameter of a circle with radius %.4f",radius);
		return false;
	}

	else
	{
		ICnt = 2;
		mirrorDistance =CMath::sqrt(CMath::pow(radius,2) - CMath::pow(separation/2,2));
		CPoint2D p1p2sum= p1+p2;
		CPoint2D p1p2dif= p1-p2;
		CPoint2D p1p2sum_h = p1p2sum /2;
		center1.set((p1p2sum_h.x() + mirrorDistance*(p1p2dif.y())/separation),(p1p2sum_h.y() + mirrorDistance*(p1p2dif.x())/separation));
		center2.set((p1p2sum_h.x() - mirrorDistance*(p1p2dif.y())/separation),(p1p2sum_h.y() - mirrorDistance*(p1p2dif.x())/separation));
		//printf("\nTwo circles are possible.");
		//printf("\nCircle C1 with center (%.4f,%.4f), radius %.4f and Circle C2 with center (%.4f,%.4f), radius %.4f",(p1.x()+p2.x())/2 + mirrorDistance*(p1.y()-p2.y())/separation,(p1.y()+p2.y())/2 + mirrorDistance*(p2.x()-p1.x())/separation,radius,(p1.x()+p2.x())/2 - mirrorDistance*(p1.y()-p2.y())/separation,(p1.y()+p2.y())/2 - mirrorDistance*(p2.x()-p1.x())/separation,radius);
	return true;
	}
}

CCircle2D::CCircle2D():CConic2D(){}

CCircle2D::CCircle2D(CPoint2D center,  float radius)
:CConic2D()
{
	set(center,radius);
}

CCircle2D::CCircle2D(CPoint2D p1,CPoint2D p2,CPoint2D p3)
:CConic2D()
{
	CTriangle2D tri(p1,p2,p3);
	CCircle2D circ=tri.circunCircle();
	CPoint2D pc=circ.getOrigin();
	set(pc,circ.getRadius());
}
CCircle2D::CCircle2D(CPoint2D p1, CPoint2D p2, float radius, int orientation)
:CConic2D()
{
	CPoint2D center1,center2;
	int count;
	if(findCircles(p1,p2,radius,count,center1,center2)) {
		if (count==1) set(center1,radius);
		else{ //count=2
			if (orientation==RightHandSide){
				set(center1,radius);
			}else set(center2,radius);
		}
	}else {//return a degenerate cirle???
		set(center1,0);
	}
}


bool CCircle2D::isFinite() const
{
	return center.isFinite() && MATH::isFinite(r);
}

bool CCircle2D::isDegenerate() const
{
	return r < 0.f;
}

void CCircle2D::set(CPoint2D &center_, float radius){
center=center_;
r=radius;
	 float a=1;
     float b=0;
	 float c=1;
	 float d=-center.x_;
	 float e=-center.y_;
	 float f=(d*d)+(e*e)-r;
//set CConic2D params
setFromPol(a,b,c,2*d,2*e,f);
//CCurveRoc::set_roc(radius);

}


float CCircle2D::sagitta(float dist) const{
	if( dist>= r) return r;
	if( dist<= 0) return 0;

	return r-dist;
}


CPoint2D CCircle2D::getPoint(float angleRadians) const
{
  float	x = center.x_ + r * CMath::cos(angleRadians);
  float y = center.y_ + r * CMath::sin(angleRadians);
  return CPoint2D(x,y);
}

CPoint2D CCircle2D::getPoint(float angleRadians, float d) const
{
  float	x = center.x_ + r * d * CMath::cos(angleRadians);
  float y = center.y_ + r * d * CMath::sin(angleRadians);
  return CPoint2D(x,y);
}



void CCircle2D::translate(const CVec2D &offset)
{
	center = center+offset;
}

void CCircle2D::transform(const CMat2D &transform)
{
//s	SMF_ASSERT(transform.hasUniformScale());
//s	SMF_ASSERT(transform.isColOrthogonal());
	center = transform.mul(center.toVec2DPos());
	r *= transform.Col(0).getLenght(); // scale the radius of the circle.
}



bool CCircle2D::edgeContains(const CPoint2D &point, float maxDistance) const
{
	return distanceToEdge(point) <= maxDistance;
}
/*
bool CCircle2D::DiscContains(const CVec2D &point, float maxDistance) const
{
	return distanceToDisc(point) <= maxDistance;
}

*/
float CCircle2D::distanceToEdge(const CPoint2D &point) const
{
	float dist = CPoint2D(center.x_,center.y_).distance(point);
	float dist2 = r-dist;
	// dist2 > 0 =>point outside circle
	// dist2 < 0 =>point inside circle
	// dist2 == 0 => point on circle edge
	return CMath::fabs(dist2);
}

bool CCircle2D::contains(const CVec2D &point) const{
	float dist = CPoint2D(center.x_,center.y_).distance(CPoint2D(point.x,point.y));
	float dist2 = r-dist;
	return dist2 <= r;
	// dist2 > 0 =>point outside circle
	// dist2 < 0 =>point inside circle
	// dist2 == 0 => point on circle edge

	}
bool CCircle2D::contains(const CPoint2D &point) const{
	bool Result;
	Result = (center.distanceSqr(point) <= (r * r));
	return Result;

}



/*
float CCircle2D::distanceToEdge(const CRay &ray, float *d, CVec2D *closestPoint) const
{
	float t;
	CVec2D cp = closestPointToEdge(ray, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(ray.getPoint(t));
}

float CCircle2D::distanceToEdge(const CLineSegment &lineSegment, float *d, CVec2D *closestPoint) const
{
	float t;
	CVec2D cp = closestPointToEdge(lineSegment, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(lineSegment.getPoint(t));
}

float CCircle2D::distanceToEdge(const CLine &line, float *d, CVec2D *closestPoint) const
{
	float t;
	CVec2D cp = closestPointToEdge(line, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(line.getPoint(t));
}
*/


bool CCircle2D::intersects(const CLine2D &line) const
{
	CPoint2D point(center.x_,center.y_);
	float dist = line.distance(point);
	return  dist <= r;
}

bool CCircle2D::intersects(const CLineSegment2D &lineSegment) const
{
	CPoint2D point(center.x_,center.y_);
	float dist = lineSegment.distance(point);
	return  dist <= r;

}

bool CCircle2D::intersects(const CRay2D &ray) const
{
	CPoint2D point(center.x_,center.y_);
	float dist;
	dist= ray.distance(point.toVec2DPos());
	return  dist <= r;
}

bool CCircle2D::intersects(const CCircle2D &c) const{
	CPoint2D p1(center);
	CPoint2D p2= c.center;
	float dist = p1.distance(p2);
	return dist <= r;
}

bool CCircle2D::intersects(const CPoint2D &p) const
{
	return contains(p);
}
bool CCircle2D::intersects(const CAABBox2D &rhs) const
{
	return rhs.intersects(*this);
}

#define USE_IT(x,val) if(x!=NULL) { *x=val;}
bool CCircle2D::intersection( const CLine2D & line,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const
{
    if(I1!=NULL) I1->toInfinite();
	if(I2!=NULL) I2->toInfinite();
	if(ICnt!=NULL) *ICnt=0;
	CPoint2D p1,p2;
	if ( CMath::fabs( line.a() ) < CMath::FLT_EPSILON )
    {
        SMF_ASSERT (! CMath::fabs( line.b() ) < CMath::FLT_EPSILON )

        // Line:    By + C = 0  ---> y = -C/B
        // Circle:  (x - cx)^2 + (y - cy)^2 = r^2
		// mx^2+nx+z=0
		//m = 1.0
		//n= -2 cx
		//z= cx^2+(-C/B)^2-(2cy*(-C/B))+cy^2-r^2
        // --->
        float x1 = 0.0, x2 = 0.0;
		float m,n,z;
		float c_b = line.c() / line.b();
		m=1.0f;
		n= -2.0 * origin().x();
		z= CMath::square(origin().x())+CMath::square( origin().y())+CMath::square(c_b)+(2*origin().y()*c_b)-CMath::square( radius() );
        int n_sol = CPolynomial::solveQuadratic( m,n,z, x1, x2 );

        if ( n_sol > 0 )
        {
            float y1 = -line.c() / line.b();

            if ( I1 )
            {
                p1.set(x1,y1);
				USE_IT(I1,p1)

            }

            if ( n_sol > 1 && I2 )
            {
                p2.set( x2, y1 );
				USE_IT(I2,p2)
            }
        }
        USE_IT(ICnt,n_sol)
		if (n_sol > 0) return true;
		else return false;
    }
    else
    {
        // include (fabs(l.b()) < EPSILON) case
        // use line & circle formula
        //   Ax + By + C = 0
        //   (x - cx)^2 + (y - cy)^2 = r^2
        // make y's quadratic formula using these fomula.
        float m = line.b() / line.a();
        float d = line.c() / line.a();

        float a = 1.0 + m * m;
        float b = 2.0 * ( -origin().y() + ( d + origin().x() ) * m );
        float c = CMath::square( d + origin().x() )
            + CMath::square( origin().y() )
            - CMath::square( radius() );

        float y1 = 0.0, y2 = 0.0;
        int n_sol = CPolynomial::solveQuadratic( a, b, c, y1, y2 );

        if ( n_sol > 0  )
        {
            p1.set( line.getX( y1 ), y1 );
			USE_IT(I1,p1)
        }

        if ( n_sol > 1  )
        {
            p2.set( line.getX( y2 ), y2 );
			USE_IT(I2,p2)
        }
		USE_IT(ICnt,n_sol)
		if (n_sol > 0) return true;
		else return false;

    }
}

/*-------------------------------------------------------------------*/
/*!

 */
bool CCircle2D::intersection( const CRay2D & ray,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const
{
    if(I1!=NULL) I1->toInfinite();
	if(I2!=NULL) I2->toInfinite();
	if(ICnt!=NULL) *ICnt=0;
	CLine2D line( ray.origin(), ray.dir );
    CPoint2D tsol1, tsol2;

    int n_sol;
	intersection( line, &n_sol,&tsol1, &tsol2 );

    if ( n_sol > 1
         && ! ray.inRightDir( tsol2, 1.0 ) )
    {
        --n_sol;
    }

    if ( n_sol > 0
         && ! ray.inRightDir( tsol1, 1.0 ) )
    {
        tsol1 = tsol2; // substituted by second solution
        --n_sol;
    }

    if ( n_sol > 0  )
    {
        USE_IT(I1,tsol1)
    }

    if ( n_sol > 1  )
    {
        USE_IT(I2,tsol2)
    }

    USE_IT(ICnt,n_sol)
	if (n_sol > 0) return true;
	else return false;
}

/*-------------------------------------------------------------------*/
/*!

 */
bool CCircle2D::intersection( const CLineSegment2D & segment,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const
{
    if(I1!=NULL) I1->toInfinite();
	if(I2!=NULL) I2->toInfinite();
	if(ICnt!=NULL) *ICnt=0;
	CLine2D line = segment.toLine2D();
    CPoint2D tsol1, tsol2;

    int n_sol;
	intersection( line, &n_sol,&tsol1, &tsol2 );

    if ( n_sol > 1
         && ! segment.recContains( tsol2 ) )
    {
        --n_sol;
    }

    if ( n_sol > 0
         && ! segment.recContains( tsol1 ) )
    {
        tsol1 = tsol2; // substituted by second solution
        --n_sol;
    }

    if ( n_sol > 0)
    {
        USE_IT(I1,tsol1)
    }

    if ( n_sol > 1)
    {
        USE_IT(I2,tsol2)
    }

    USE_IT(ICnt,n_sol)
	if (n_sol > 0) return true;
	else return false;
}

/*-------------------------------------------------------------------*/
/*!

 */
bool CCircle2D::intersection( const CCircle2D & circle,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const
{
    float rel_x = circle.origin().x() - this->origin().x();
    float rel_y = circle.origin().y() - this->origin().y();

    float center_dist2 = rel_x * rel_x + rel_y * rel_y;
    float center_dist = CMath::sqrt( center_dist2 );

    if ( center_dist < CMath::fabs( this->radius() - circle.radius() )
         || this->radius() + circle.radius() < center_dist )
    {
        return false;
    }

    //std::cerr << "must exist intersection C1: " << this->origin() << this->radius()
    //        << " C2: " << circle.origin() << circle.radius()
    //        << std::endl;
    // line that passes through the intersection points
    CLine2D line( -2.0 * rel_x,
                 -2.0 * rel_y,
                 circle.origin().toVec2DPos().getLengthSqr()
                 - circle.radius() * circle.radius()
                 - this->origin().toVec2DPos().getLengthSqr()
                 + this->radius() * this->radius() );

    return this->intersection( line,ICnt, I1, I2 );
}
#undef USE_IT(x,val)

std::string CCircle2D::toString() const
{
	char str[256];
	std::sprintf(str, "CCircle2D(center:(%.2f, %.2f) , r:%.2f)",
		center.x_, center.y_,   r);
	return str;
}
namespace _2D{
std::ostream &operator <<(std::ostream &o, const CCircle2D &circle)
{
	o << circle.toString();
	return o;
}
} //end _2D

#if 0
CCircle2D operator *(const CMat2D &transform, const CCircle2D &circle)
{
	CCircle2D c(circle);
	c.transform(transform);
	return c;
}
#endif

} //end GEO
}  //end SMF
