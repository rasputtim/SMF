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

float calcSagitta(CPoint2D &start, CPoint2D &end, float radius){
	float l2 = CMath::square(start.distance(end)*0.5), sag;
	float r2 = CMath::square(radius);
	sag = radius - (CMath::sqrt(r2-l2));
	return sag;
}

CCircleSegment2D::CCircleSegment2D(CPoint2D &start, CPoint2D &end, CCircle2D &circ){
	set(start, end,circ);
}

CCircle2D CCircleSegment2D::getSupCircle() const{CPoint2D pc= getOrigin();return CCircle2D(pc,getRadius());};


CCircleSegment2D::CCircleSegment2D(CPoint2D &start, CPoint2D &end, float radius, int orientation){
	CCircle2D circ(start,end,radius,orientation);
	set(start,end,circ);
}


void CCircleSegment2D::set(CPoint2D &start, CPoint2D &end, CCircle2D &circ){

startPoint=start;
endPoint=end;
if (circ.isDegenerate()){ //create a degenerate segment
	_sag=0;
	set_roc(0);
	circ.origin().toInfinite();
	CPoint2D pc=circ.origin();
	CCircle2D::set(pc,0);

}else{
	if (circ.edgeContains(start) && circ.edgeContains(end)){
		_sag=calcSagitta(start,end, circ.r);
		set_roc(circ.r);
		CPoint2D pc=circ.getOrigin();
		CCircle2D::set(pc,circ.getRadius());
	}else{ // circle is not possible.create a degenerate segment

		_sag=0;
		set_roc(0);
		getOrigin().toInfinite();
		CPoint2D pc=getOrigin();
		CCircle2D::set(pc,0);
		startPoint.toInfinite();
		endPoint.toInfinite();
	}
}
}

//s = r - (dist+square(r^2-l^2))
float CCircleSegment2D::sagitta(float dist) const{
	if( dist>= _sag) return _sag;
	if( dist<= 0) return 0;
	//the origin of the segment
	float r2= CMath::square(r);
	float l2= CMath::square(startPoint.distance(endPoint)*0.5);
	float sag = r - (dist-CMath::sqrt(r2-l2));
	return sag;
}

   /**A= r* Arccos((r-h)/r)-(r-h)*sqrt(2rh-h^2)
	r  is the radius of the circle of which the segment is a part.
	h  is the height of the segment (sagitta).
	**/
float CCircleSegment2D::area()const{
  float h2=CMath::square(_sag);
  float r2=CMath::square(getRadius());
  float r_sag=r-_sag;
  float arccos=CMath::acos(((r_sag)/r));
  float result= (r * arccos)-(r_sag*CMath::sqrt((2*r*_sag)-h2));
  return result;
}

float CCircleSegment2D::centralAngle()const{

	CVec2D p1=(startPoint-getOrigin());
	CVec2D p2=(endPoint-getOrigin());
	return p1.angleBetween(p2);

}
float CCircleSegment2D::getLenght()const{
	return CCircle2D::arcLenght(centralAngle());

}
/**
	\brief Calculating Height of an Arc at Any Point
	\note h = s + sqrt(r^2-x^2)-r
	where:
	h = the height of the arc;
	s = the sagitta of the arc;
	r = the radius of the arc;
	x = the horizontal offset from the center to the point where you want the height;
	**/
float CCircleSegment2D::height(CPoint2D &point)const{
	//todo: check if te point is on the arc

	CPoint2D center= _2D::midpoint(startPoint,endPoint);
	float x = CMath::fabs((point-center).getLenght());
	float x2 = CMath::square(x);
	float r2=CMath::square(r);
	float h = _sag + CMath::sqrt((r2-x2))-r;
	return h;
}

std::string CCircleSegment2D::toString() const
{
	char str[256];
	std::sprintf(str, "CCircleSegment2D(center:(%.2f, %.2f) , r:%.2f, startPoint( %.2f,%.2f), endoint(%.2f,%2.f))",
		center.x_, center.y_,   r,startPoint.x(),startPoint.y(),endPoint.x(),endPoint.y());
	return str;
}
namespace _2D{
std::ostream &operator <<(std::ostream &o, const CCircleSegment2D &segment)
{
	o << segment.toString();
	return o;
}
} //end _2D

} //end GEO
}  //end SMF
