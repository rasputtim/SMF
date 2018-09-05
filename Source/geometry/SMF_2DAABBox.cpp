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
#include "geometry/SMF_2DAABBox.h"
#include "math/all.h"
#include "geometry/all.h"
#include <vector>

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{
// Constructors/Destructor---------------------------------------------------


CAABBox2D::CAABBox2D()
{
  min_pos_[0]=min_pos_[1]=(float)1;
  max_pos_[0]=max_pos_[1]=(float)0; // empty box
}


CAABBox2D::CAABBox2D(float const corner1[2],
                             float const corner2[2])
{
  min_pos_[0]=max_pos_[0]=corner1[0];
  min_pos_[1]=max_pos_[1]=corner1[1];
  this->add(corner2);
}


CAABBox2D::CAABBox2D(CPoint2D const& corner1,
                             CPoint2D const& corner2)
{
  min_pos_[0]=max_pos_[0]=corner1.x();
  min_pos_[1]=max_pos_[1]=corner1.y();
  this->add(corner2);
}


CAABBox2D::CAABBox2D(float xmin, float xmax, float ymin, float ymax)
{
  min_pos_[0]=max_pos_[0]=xmin;
  min_pos_[1]=max_pos_[1]=ymin;
  this->add(CPoint2D(xmax,ymax));
  if (xmin > xmax || ymin > ymax) this->empty();
}


CAABBox2D::CAABBox2D(float const ref_point[2],
                             float w, float h,
                             CAABBox2D::point_Type t)
{
  if (t == CAABBox2D::centre)
  {
    min_pos_[0]=float(ref_point[0]-0.5*w);
    min_pos_[1]=float(ref_point[1]-0.5*h);
    max_pos_[0]=float(ref_point[0]+0.5*w);
    max_pos_[1]=float(ref_point[1]+0.5*h);
  }
  else if (t == CAABBox2D::min_pos)
  {
    min_pos_[0]=ref_point[0];
    min_pos_[1]=ref_point[1];
    max_pos_[0]=ref_point[0]+w;
    max_pos_[1]=ref_point[1]+h;
  }
  else if (t == CAABBox2D::max_pos)
  {
    min_pos_[0]=ref_point[0]-w;
    min_pos_[1]=ref_point[1]-h;
    max_pos_[0]=ref_point[0];
    max_pos_[1]=ref_point[1];
  }
  else
    SMF_ASSERT(!"point_type should be one of: centre, min_pos, max_pos");
}


CAABBox2D::CAABBox2D(CPoint2D const& ref_point,
                             float w, float h,
                             CAABBox2D::point_Type t)
{
  if (t == CAABBox2D::centre)
  {
    min_pos_[0]=float(ref_point.x()-0.5*w);
    min_pos_[1]=float(ref_point.y()-0.5*h);
    max_pos_[0]=float(ref_point.x()+0.5*w);
    max_pos_[1]=float(ref_point.y()+0.5*h);
  }
  else if (t == CAABBox2D::min_pos)
  {
    min_pos_[0]=ref_point.x();
    min_pos_[1]=ref_point.y();
    max_pos_[0]=ref_point.x()+w;
    max_pos_[1]=ref_point.y()+h;
  }
  else if (t == CAABBox2D::max_pos)
  {
    min_pos_[0]=ref_point.x()-w;
    min_pos_[1]=ref_point.y()-h;
    max_pos_[0]=ref_point.x();
    max_pos_[1]=ref_point.y();
  }
  else
    SMF_ASSERT(!"point_type should be one of: centre, min_pos, max_pos");
}


float CAABBox2D::getCentroidX() const
{
  SMF_ASSERT(!isEmpty());
  float centroid= float(0.5*(min_pos_[0] + max_pos_[0]));
  return centroid;
}


float CAABBox2D::getCentroidY() const
{
  SMF_ASSERT(!isEmpty());
  float centroid=float(0.5*(min_pos_[1] + max_pos_[1]));
  return centroid;
}


float CAABBox2D::width() const
{
  return (max_pos_[0] > min_pos_[0]) ? max_pos_[0] - min_pos_[0] : 0;
}


float CAABBox2D::height() const
{
  return (max_pos_[1] > min_pos_[1]) ? max_pos_[1] - min_pos_[1] : 0;
}


CPoint2D CAABBox2D::minPoint() const
{
  SMF_ASSERT(!isEmpty());
  return CPoint2D(min_pos_[0],min_pos_[1]);
}


CPoint2D CAABBox2D::maxPoint() const
{
  SMF_ASSERT(!isEmpty());
  return CPoint2D(max_pos_[0],max_pos_[1]);
}


CPoint2D CAABBox2D::centroid() const
{
  SMF_ASSERT(!isEmpty());
  return CPoint2D(getCentroidX(),getCentroidY());
}


void CAABBox2D::setCentroidX(float cent_x)
{
  SMF_ASSERT(!isEmpty());
  float delta = cent_x - getCentroidX();
  min_pos_[0]= min_pos_[0] + delta;
  max_pos_[0]= max_pos_[0] + delta;
}


void CAABBox2D::setCentroidY(float cent_y)
{
  SMF_ASSERT(!isEmpty());
  float delta = cent_y - getCentroidY();
  min_pos_[1]= min_pos_[1] + delta;
  max_pos_[1]= max_pos_[1] + delta;
}


// All this code is to avoid drift in the centroid.
// int version
SMF_INLINE void set_dim(int & minv, int& maxv, int spread)
{
  int sum = minv + maxv;
  sum = sum | (spread & 1); // if width is odd, then make sum odd
  minv = float(CMath::floor((sum-spread)/2.0));
  maxv = minv+spread;
}
//float version
SMF_INLINE void set_dim(float  & minv, float& maxv, float spread)
{
  float x = minv + maxv;
  minv = float( (x-spread)*0.5 );
  maxv = minv + spread;
}
#if 0
SMF_INLINE void set_dim(float & minv, float& maxv, float spread)
{
  float sum = minv + maxv;
  minv = float( (sum-spread)*0.5 );
  maxv = minv + spread;
}
#endif
//: Modify width, retaining centroid at current position
// For integer types, centroid might change slightly, but
// repeat calls to setHeight will not cause centroid drift.

void CAABBox2D::setWidth(float w)
{
  SMF_ASSERT(!isEmpty());
  set_dim(min_pos_[0], max_pos_[0], w);
}

//: Modify height, retaining centroid at current position
// For integer types, centroid might change slightly, but
// repeat calls to setHeight will not cause centroid drift.

void CAABBox2D::setHeight(float h)
{
  SMF_ASSERT(!isEmpty());
  set_dim(min_pos_[1], max_pos_[1], h);
}


//: add to width and height, centroid unchanged.
// Will move each side by \p expand / 2.

void CAABBox2D::expandAboutCentroid(float expand)
{
  SMF_ASSERT(!isEmpty());
  set_dim(min_pos_[0], max_pos_[0], width() + expand );
  set_dim(min_pos_[1], max_pos_[1], height() + expand );
}

//: scale width and height, centroid unchanged.

void CAABBox2D::scaleAboutCentroid(double s)
{
  SMF_ASSERT(!isEmpty());
  set_dim(min_pos_[0], max_pos_[0], static_cast<float>(width()*s));
  set_dim(min_pos_[1], max_pos_[1], static_cast<float>(height()*s));
}


//: scale width and height, keeping scaled position of origin unchanged.

void CAABBox2D::scaleAboutOrigin(double s)
{
  min_pos_[0] = static_cast<float>(min_pos_[0] * s);
  min_pos_[1] = static_cast<float>(min_pos_[1] * s);
  max_pos_[0] = static_cast<float>(max_pos_[0] * s);
  max_pos_[1] = static_cast<float>(max_pos_[1] * s);
}


void CAABBox2D::setMinPosition(float const min_position[2])
{
  min_pos_[0]=min_position[0];
  min_pos_[1]=min_position[1];
  if (max_pos_[0] < min_pos_[0]) {
    max_pos_[0]=min_pos_[0];
  }
  if (max_pos_[1] < min_pos_[1]) {
    max_pos_[1]=min_pos_[1];
  }
}


void CAABBox2D::setMaxPosition(float const max_position[2])
{
  max_pos_[0]=max_position[0];
  max_pos_[1]=max_position[1];
  if (max_pos_[0] < min_pos_[0])
    min_pos_[0]=max_pos_[0];
  if (max_pos_[1] < min_pos_[1])
    min_pos_[1]=max_pos_[1];
}


void CAABBox2D::setMinPoint(CPoint2D const& min_pt)
{
  min_pos_[0]=min_pt.x(); if (max_pos_[0]<min_pos_[0]) max_pos_[0]=min_pos_[0];
  min_pos_[1]=min_pt.y(); if (max_pos_[1]<min_pos_[1]) max_pos_[1]=min_pos_[1];
}


void CAABBox2D::setMaxPoint(CPoint2D const& max_pt)
{
  max_pos_[0]=max_pt.x(); if (max_pos_[0]<min_pos_[0]) min_pos_[0]=max_pos_[0];
  max_pos_[1]=max_pt.y(); if (max_pos_[1]<min_pos_[1]) min_pos_[1]=max_pos_[1];
}


bool CAABBox2D::intersection(CLine2D const& line, int *ICnt, CPoint2D* p0, CPoint2D* p1)const
{
  return line.intersection(*this,ICnt,p0,p1);
}

bool CAABBox2D::intersects(const CLine2D& line)const 
{
	return line.intersects(*this); 
};


bool CAABBox2D::intersects(const CAABBox2D &rhs) const
{
		return max_pos_[0] >= rhs.min_pos_[0] &&
		       max_pos_[1] >= rhs.min_pos_[1] &&
		       rhs.max_pos_[0] >= min_pos_[0] &&
		       rhs.max_pos_[1] >= min_pos_[1];
}



bool CAABBox2D::intersects(const CLineSegment2D &segment)const{

  bool Result = intersects(segment.toCAABBox());
  return Result;
}

bool CAABBox2D::intersects(const CTriangle2D &tri) const
{

  bool Result = intersects(tri.edge(tri.side1)) ||
                intersects(tri.edge(tri.side2)) ||
                intersects(tri.edge(tri.side3));

  return Result;
}


bool CAABBox2D::intersects(const CCircle2D &circle)const {
bool Result;
  Result = circle.contains(ClosestPointOnRectangleFromPoint(circle.center));
return Result;
} 

CPoint2D CAABBox2D::ClosestPointOnRectangleFromPoint(const CPoint2D &point)const {

	float Nx,Ny,x1,y1,x2,y2,Px,Py;
	Px=point.x();
	Py=point.y();
	x1= minX();
	x2= maxX();
	y1= minY();
	y2= maxY();

	if (Px < MIN(x1,x2)) Nx = MIN(x1,x2);
  else if (Px > MAX(x1,x2)) Nx = MAX(x1,x2);
  else Nx = Px;

  if (Py < MIN(y1,y2)) Ny = MIN(y1,y2);
  else if (Py > MAX(y1,y2)) Ny = MAX(y1,y2);
  else Ny = Py;
  return CPoint2D(Nx,Ny);
}  


bool CAABBox2D::intersection(const CLineSegment2D& line_seg, int *ICnt, CPoint2D * p0, CPoint2D* p1)
{
#define USE_IT(x,val) if(x!=NULL) { *x=val;}
	if (p0!=NULL) p0->toInfinite();
	if (p1!=NULL) p1->toInfinite();
	if (ICnt!=NULL) *ICnt=0;
	int nint2,nint = 0;
	CPoint2D pi0, pi1;
  
	CLine2D line(line_seg.a(), line_seg.b(), line_seg.c());
	// if no intersection just return
	if (!intersection(line, &nint2,&pi0, &pi1)){
		USE_IT(ICnt,nint)
    return false;
	}
  // check if intersection points are interior to the line segment
  if (line_seg.isPointOnSegment(pi0)) {
    USE_IT(p0,pi0)
    nint++;

  }
  if (line_seg.isPointOnSegment(pi1)) {
    USE_IT(p1,pi1)
    nint++;
  }
  USE_IT(ICnt,nint)
  return true;
#undef USE_IT(x,val)
}


bool CAABBox2D::intersection(const CPolygon2D& poly)
{
	const CAABBox2D& b = *this; 
  // easy checks first
  // check if any poly vertices are inside the box
  unsigned int ns = poly.numVertices();
  bool hit = false;
    for (unsigned int i = 0; i<ns&&!hit; ++i) {
      CPoint2D p = poly[i];
      hit = b.contains(p.x(), p.y());
    }
  
  if (hit) return true;
  // check if any box vertices are inside the polygon
  float minx = b.minX(), maxx = b.maxX();
  float miny = b.minY(), maxy = b.maxY();
  hit = poly.contains(CPoint2D(minx, miny)) || poly.contains(CPoint2D(maxx, maxy)) ||
    poly.contains(CPoint2D(minx, maxy)) || poly.contains(CPoint2D(maxx, miny));
  if (hit) return true;
  // check if any polygon edges intersect the box
  
    CPoint2D ia, ib;
    CPoint2D last = poly[0];
    for (unsigned int i = 1; i<ns&&!hit; ++i)
    {
      CPoint2D p = poly[i];
      CLineSegment2D l(last, p);
	  int count;
      intersection(l, &count,&ia, &ib);
	  hit = count > 0;
      last = p;
    }
    if (!hit) {
      CPoint2D start = poly[0];
	  int count;
      CLineSegment2D ll(last, start);
	  intersection(ll, &count, &ia, &ib);
      hit = count >0;
    }
  
  return hit;
}

std::vector< CPoint2D > CAABBox2D::intersection( std::vector< CPoint2D > const& p)
{ 

  std::vector< CPoint2D > r;
  std::vector< CPoint2D >::const_iterator i;
  for (i = p.begin(); i != p.end(); ++i)
    if (intersects(*i))
      r.push_back(*i);
  return r;
}









ostream& CAABBox2D::print(ostream& s) const
{
  if (isEmpty())
    return s << "<CAABBox2D (empty)>";
  else
    return s << "<CAABBox2D "
             << min_pos_[0] << ',' << min_pos_[1] << " to "
             << max_pos_[0] << ',' << max_pos_[1] << '>';
}


ostream& CAABBox2D::write(ostream& s) const
{
  return s << min_pos_[0] << ' ' << min_pos_[1] << ' '
           << max_pos_[0] << ' ' << max_pos_[1] << '\n';
}


istream& CAABBox2D::read(istream& s)
{
  return s >> min_pos_[0] >> min_pos_[1]
           >> max_pos_[0] >> max_pos_[1];
}
//: add a point to this box.
// Do this by possibly enlarging the box so that the point just falls within the box.
// Adding a point to an empty box makes it a size zero box only containing p.

void CAABBox2D::add(CPoint2D const& p)
{
  if (isEmpty())
  {
    min_pos_[0] = max_pos_[0] = p.x();
    min_pos_[1] = max_pos_[1] = p.y();
  }
  else
  {
    if (p.x() > max_pos_[0]) max_pos_[0] = p.x();
    if (p.x() < min_pos_[0]) min_pos_[0] = p.x();
    if (p.y() > max_pos_[1]) max_pos_[1] = p.y();
    if (p.y() < min_pos_[1]) min_pos_[1] = p.y();
  }
}

//: Make the convex union of two boxes
// Do this by possibly enlarging this box so that the corner points of the
// given box just fall within the box.
// Adding an empty box does not change the current box.

void CAABBox2D::add(CAABBox2D const& b)
{
  if (b.isEmpty()) return;
  add(b.minPoint());
  add(b.maxPoint());
}

//: Return true iff the point p is inside this box

bool CAABBox2D::contains(CPoint2D const& p) const
{
  return contains(p.x(), p.y());
}

//: Return true if the corner points of b are inside this box

bool CAABBox2D::contains(CAABBox2D const& b) const
{
  return
    contains(b.min_pos_[0], b.min_pos_[1]) &&
    contains(b.max_pos_[0], b.max_pos_[1]);
}

//: Make the box empty

void CAABBox2D::empty()
{
  min_pos_[0]=min_pos_[1]=(float)1;
  max_pos_[0]=max_pos_[1]=(float)0;
}

//: Return the intersection of two boxes (which is itself is a box, possibly the empty box)
CAABBox2D CAABBox2D::intersection(CAABBox2D const& b2)
{
	CAABBox2D const& b1=*this;
  float xmin = b1.min_pos_[0] > b2.min_pos_[0] ? b1.min_pos_[0] : b2.min_pos_[0];
  float ymin = b1.min_pos_[1] > b2.min_pos_[1] ? b1.min_pos_[1] : b2.min_pos_[1];
  float xmax = b1.max_pos_[0] < b2.max_pos_[0] ? b1.max_pos_[0] : b2.max_pos_[0];
  float ymax = b1.max_pos_[1] < b2.max_pos_[1] ? b1.max_pos_[1] : b2.max_pos_[1];
  return CAABBox2D(xmin,xmax,ymin,ymax);
}



COBBox CAABBox::toOBBox() const
{
	return COBBox(*this);
}

//: Print to stream

ostream& operator<<(ostream& s, CAABBox2D const& p)
{
  return p.print(s);
}

//: Read from stream

istream& operator>>(istream& is, CAABBox2D& p)
{
  return p.read(is);
}

namespace _2D{

} //end 2D
} //end GEO
}  //end SMF
