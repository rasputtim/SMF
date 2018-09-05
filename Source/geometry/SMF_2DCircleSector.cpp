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

CCircSector2D::CCircSector2D(CPoint2D &start, CPoint2D &end, CCircle2D &circ){
	set(start, end,circ);
}

CCircSector2D::CCircSector2D(CPoint2D &start, CPoint2D &end, float radius, int orientation){
	CCircle2D circ(start,end,radius,orientation);
	set(start,end,circ);
}

CCircle2D CCircSector2D::getSupCircle() const{ CPoint2D pc=getOrigin(); CCircle2D circ; circ.set(pc,getRadius());
	    return circ;}; //CCircle2D(cent,);};

float CCircSector2D::sagitta(float dist) const{
	if( dist>= r) return r;
	if( dist<= 0) return 0;

	if( dist<= 0) return 0;
	return r-dist;
}


void CCircSector2D::set(CPoint2D &start, CPoint2D &end, CCircle2D &circ){

startPoint=start;
endPoint=end;
if (circ.isDegenerate()){ //create a degenerate segment
	set_roc(0);
	circ.origin().toInfinite();
	CPoint2D pc = circ.origin();
	CCircle2D::set(pc,0);

}else{
	if (circ.edgeContains(start) && circ.edgeContains(end)){
		set_roc(circ.r);
		CPoint2D pc = circ.getOrigin();
		CCircle2D::set(pc,circ.getRadius());
	}else{ // circle is not possible.create a degenerate segment


		set_roc(0);
		getOrigin().toInfinite();
		CPoint2D pc=getOrigin();
		CCircle2D::set(pc,0);
		startPoint.toInfinite();
		endPoint.toInfinite();
	}
}
}

float CCircSector2D::centralAngle()const{

	CVec2D p1=(startPoint-getOrigin());
	CVec2D p2=(endPoint-getOrigin());
	return p1.angleBetween(p2);

}
float CCircSector2D::getLenght()const{
	return CCircle2D::arcLenght(centralAngle());

}

bool areClockwise(const CPoint2D &v1, const CPoint2D &v2) {
  return -v1.x()*v2.y() + v1.y()*v2.x() > 0;
}
bool CCircSector2D::isInsideSector(const CPoint2D &point, const CPoint2D &center, const CPoint2D &sectorStart, const CPoint2D &sectorEnd)const {

CPoint2D relPoint= point - center;

  return !areClockwise(sectorStart, relPoint) &&
         areClockwise(sectorEnd, relPoint) &&
         CCircle2D::contains(point);
}

bool CCircSector2D::contains(const CVec2D &point) const{
	return isInsideSector(CPoint2D(point),getOrigin(),startPoint,endPoint);
}
bool CCircSector2D::contains(const CPoint2D &point) const{
	return isInsideSector(point,getOrigin(),startPoint,endPoint);
}


std::string CCircSector2D::toString() const
{
	char str[256];
	std::sprintf(str, "CCircSector2D(center:(%.2f, %.2f) , r:%.2f, startPoint( %.2f,%.2f ): endPoint(%2.f,%2.f)",
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
