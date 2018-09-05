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
#include "geometry/SMF_2DEllipse.h"
#include "geometry/SMF_2DLineSegment.h"
#include "geometry/SMF_2DLine.h"
#include "geometry/SMF_2DConics.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{

CEllipse2D::CEllipse2D():CConic2D(){}

CEllipse2D::CEllipse2D(CPoint2D center,  float radius1,  float radius2, float angle)
:CConic2D()
{
	set(center,radius1,radius2,angle);
}



void CEllipse2D::set(CPoint2D center_,float rx_,  float ry_, float anglePhi){
center=center_;
rx= rx_;
ry= ry_;
set_ang(anglePhi);
xy_ratio = rx / ry;
//- Schwarzschild constant
//_e2 = CMath::square64(CMath::sqrt64(CMath::fabs((rx * rx) - (ry*ry))) / MAX(rx,ry));

//angle= angle_;

/*
phi=0
center(h,k)
A = 1 / rx^2
B = 0   .'.  b = 2 * B=0
C = 1 / ry^2;
D = - h/rx^2   .'.  d = 2 * D
E = - k/ry^2   .'.  e = 2 * E
F = (h^2/rx^2)+(k^2/ry^2)-1
*/
     float a,b,c,d,e,f,h,k;
	  h = center.x();
	  k = center.y();
	 float maxRadio2=CMath::square(this->a());
	 float minRadio2 = CMath::square(this->b());
	 if(!CMath::nearZero(Phi())) {
			a = 1/maxRadio2;
			b =0;
			c = 1/minRadio2;
			d = -2*(h/maxRadio2);
			e = -2*(k/minRadio2);

	  }else{
		  float cosPhi = CMath::cos(Phi());
		  float cosPhi2 = CMath::square(cosPhi);
		  float sinPhi = CMath::sin(Phi());
		  float sinPhi2 = CMath::square(sinPhi);

			a = (cosPhi2/maxRadio2)+(sinPhi2/minRadio2);
			c = (cosPhi2/minRadio2)+(sinPhi2/maxRadio2);
			b = cosPhi*sinPhi*((1/minRadio2)-(1/maxRadio2));
			d = 2*((-h*cosPhi)/maxRadio2)-((k*sinPhi)/minRadio2);
			e = 2*((-k*cosPhi)/minRadio2)+((h*sinPhi)/maxRadio2);
	 }


	 f= (((h*h)/maxRadio2)+((k*k)/minRadio2)-1);
	 //Set conic2D params
	 this->setFromPol(a,b,c,d,e,f);
	 //set roc


}

float CEllipse2D::getEccentricity()const{
	float ec=CMath::sqrt(CMath::fabs(CMath::square(a())-CMath::square(b())))/a();
	return ec;
}
//return the focci length (the distance from focus1 t focus2)
// focci = 2c. and c = Sqroot(rmajor^2-b()^2)
float CEllipse2D::getFocci() const{
	float ec=CMath::sqrt(CMath::fabs(CMath::square(a())-CMath::square(b())))*2;
	return ec;
}
// Return one of the focus of this elipse
CPoint2D CEllipse2D::getFocus1()const{
	return rx > ry ? CPoint2D(center.x()+c(), center.y()) : CPoint2D(center.x(),center.y()+c());
}
CPoint2D CEllipse2D::getFocus2()const{
	return rx > ry ? CPoint2D(center.x()-c(), center.y()) : CPoint2D(center.x(),center.y()-c());

}
float CEllipse2D::getY1(float x)const{
	float raiz1,raiz2;
	float invB2= 1/CMath::square(b());
	float invA2= 1/CMath::square(a());
	float h = center.x();
	float k = center.y();
	float A = 1 *invB2;
	float B = -(2*h*invB2);
	float C =(CMath::square(k)*invB2)+(CMath::square(x-h)*invA2)-1;
	int raizes = CPolynomial::solveQuadratic(A,B,C,raiz1,raiz2);
	if (raizes==0) return CMath::NAN_FLOAT;
	if (raizes==1) return raiz1;
	if (raizes==2) return raiz1;

}
float CEllipse2D::getY2(float x)const{
	float raiz1,raiz2;
	float invB2= 1/CMath::square(b());
	float invA2= 1/CMath::square(a());
	float h = center.x();
	float k = center.y();
	float A = 1 *invB2;
	float B = -(2*h*invB2);
	float C =(CMath::square(k)*invB2)+(CMath::square(x-h)*invA2)-1;
	int raizes = CPolynomial::solveQuadratic(A,B,C,raiz1,raiz2);
	if (raizes==0) return CMath::NAN_FLOAT;
	if (raizes==1) return raiz1;
	if (raizes==2) return raiz2;


}
double CEllipse2D::getRadius(const float angle) const
{
	double cosang=CMath::cos64(angle);
	double sinang=CMath::sin64(angle);
	double dist= a()*b()/(CMath::sqrt64((CMath::square64(a())*sinang)+(CMath::square64(b())*cosang)));
	return CMath::round(dist); //round to float to be compatible with other methods and classes
}

double CEllipse2D::getRadius(const CVec2D direction) const
{
	//todo verificar se já está normalizado???
	CVec2D dir = direction;
	//- Schwarzschild constant
	float _e2 = get_schwarzschild();
	if (!dir.isNormalized()) dir.toNormal();
     return rx > ry
        ? CMath::round(CMath::sqrt(CMath::square(ry) / (1. - _e2 * CMath::square(dir.x))))
        : CMath::round(CMath::sqrt(CMath::square(rx) / (1. - _e2 * CMath::square(dir.y))));
}
bool CEllipse2D::contains(const CVec2D &point) const{
	float num =(CMath::square(point.x-center.x())/CMath::square(a()))+(CMath::square(point.y-center.y())/CMath::square(b()));
	return num <=1;
	//return contains(CPoint2D(point.x,point.y));
}
bool CEllipse2D::contains(const CPoint2D &point) const{
	float num =(CMath::square(point.x()-center.x())/CMath::square(a()))+(CMath::square(point.y()-center.y())/CMath::square(b()));
	return num <=1;

	//return contains(point);
}
bool CEllipse2D::edgeContains(const CPoint2D &point, float maxDistance) const
{
		float num =(CMath::square(point.x()-center.x())/CMath::square(a()))+(CMath::square(point.y()-center.y())/CMath::square(b()));
		return num >= 1- maxDistance && num <= 1+ maxDistance;

}
/*
bool CEllipse2D::DiscContains(const CVec2D &point, float maxDistance) const
{
	return distanceToDisc(point) <= maxDistance;
}

*/
float CEllipse2D::distanceToEdge(const CPoint2D &point) const
{

	double dist = center.distance(point);
	double radius = getRadius(CVec2D(point.x(),point.y()));
	float dist2 = radius-dist;
	// dist2 > 0 =>point outside ellipse
	// dist2 < 0 =>point inside ellipse
	// dist2 == 0 => point on ellipse edge
	return CMath::fabs(dist2);
}

std::string CEllipse2D::toString() const
{
	char str[256];
	std::sprintf(str, "ELLIPSE> center:(%6.3f, %6.3f)\n"
		"\tMajor Radius:%6.3f \n"
		"\tMinor Radius:%6.3f \n"
		"\tRotation: %6.3f degrees\n",
		center.x(), center.y(),  a(),b(),180.0*Phi()/CMath::PI);
	return str;
}

CPoint2D CEllipse2D::getPoint(float angleRadians) const
{
  if (Phi()==0) {
  float	x = center.x() + getRadius(angleRadians) * CMath::cos(angleRadians);
  float y = center.y() + getRadius(angleRadians) * CMath::sin(angleRadians);
  return CPoint2D(x,y);
  }else {
	  return CPoint2D(CMath::NAN_FLOAT,CMath::NAN_FLOAT);
  }
}
CPoint2D CEllipse2D::getPoint(CVec2D direction) const
{
  float angleRadians = direction.getAngle(CVec2D(1,0));

  float	x = center.x() + getRadius(angleRadians) * CMath::cos(angleRadians);
  float y = center.y() + getRadius(angleRadians) * CMath::sin(angleRadians);
  return CPoint2D(x,y);
}
CPoint2D CEllipse2D::getPoint(CVec2D direction, float d) const
{
  //angulo entre a direção e o eixo x
  float angleRadians = direction.getAngle(CVec2D(getPoint(0.0f)));

  float	x = center.x() + getRadius(angleRadians) * d * CMath::cos(angleRadians);
  float y = center.y() + getRadius(angleRadians) * d * CMath::sin(angleRadians);
  return CPoint2D(x,y);
}

void CEllipse2D::translate(const CVec2D &offset)
{
	using _2D::operator+=;
	center += offset;
}

bool CEllipse2D::intersectsDisc(const CLine2D &line) const
{
	CPoint2D int_point,int_point2;
	bool res= intersection(*this,line,int_point,int_point2);
	return res;
}

IntersectResult	CEllipse2D::intersection(const CEllipse2D & elipse,const CLine2D &line,CPoint2D &int_point1,CPoint2D &int_point2){
	// calcular constantes
	CPoint2D p1,p2;
	line.get_two_points(p1,p2);
	float x1 = p1.x();
	float x2 = p2.x();
	float y1 = p1.y();
	float y2 = p2.y();
	float  h = elipse.center.x();
	float  k = elipse.center.y();
	//constants
	float invA= 1/elipse.a();
	float invB= 1/elipse.b();
	float x2_x1DivA=(x2-x1)*invA;
    float y2_y1DivB=(y2-y1)*invB;
	float A,B,C;
	A= (CMath::square(x2_x1DivA)+CMath::square(y2_y1DivB));
	B= ((2*x1*invA*x2_x1DivA)-(2*h*invA*x2_x1DivA)+(2*y1*invB*y2_y1DivB)-(2*k*invB*y2_y1DivB));
	C= ( (CMath::square((x1-h))*CMath::square(invA))+(CMath::square((y1-k))*CMath::square(invB)))-1;
    float raiz1,raiz2;
	int raizes= CPolynomial::solveQuadratic(A,B,C,raiz1,raiz2);
	if (raizes==0) return NOT_INTERESECTING;
	if (raizes==1) {
		float yy1 = line.getY(raiz1);
		if (MATH::isFinite(yy1)) {
		int_point1.set(raiz1,yy1);
		}else{
			yy1=elipse.getY1(raiz1);
		}
		return INTERESECTING_EXTREMITY_P1;
	}
	if (raizes==2) {
		float yy1 = line.getY(raiz1);
		if (MATH::isFinite(yy1)) {
		int_point1.set(raiz1,yy1);
		}else{

			yy1=elipse.getY1(raiz1+raiz2);
			int_point1.set(raiz1+raiz2,yy1);

		}
		float yy2 = line.getY(raiz2);
		if (MATH::isFinite(yy2)) {
		int_point2.set(raiz2,yy2);
		}else{
			yy2=elipse.getY2(raiz1+raiz2);
			int_point2.set(raiz1+raiz2,yy2);
		}
		return INTERESECTING;
	}
}


namespace _2D{
std::ostream &operator <<(std::ostream &o, const CEllipse2D &ellipse)
{
	o << ellipse.toString();
	return o;
}
}



#if 0







void CEllipse2D::transform(const CMat2D &transform)
{
//s	SMF_ASSERT(transform.hasUniformScale());
//s	SMF_ASSERT(transform.isColOrthogonal());
	center = transform.mul(center);
	r *= transform.Col(0).getLenght(); // scale the radius of the ellipse.
}








/*
float CEllipse2D::distanceToEdge(const CRay &ray, float *d, CVec2D *closestPoint) const
{
	float t;
	CVec2D cp = closestPointToEdge(ray, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(ray.getPoint(t));
}

float CEllipse2D::distanceToEdge(const CLineSegment &lineSegment, float *d, CVec2D *closestPoint) const
{
	float t;
	CVec2D cp = closestPointToEdge(lineSegment, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(lineSegment.getPoint(t));
}

float CEllipse2D::distanceToEdge(const CLine &line, float *d, CVec2D *closestPoint) const
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




bool CEllipse2D::intersectsDisc(const CLineSegment2D &lineSegment) const
{
	CPoint2D point(center.x,center.y);
	float dist = lineSegment.distance(point);
	return  dist <= r;

#if 0
bool CEllipse2D::intersectsDisc(const CRay2D &ray) const
{
	float d;
	return ray.getPoint(d).distanceSq(center) <= r*r;
}
#endif
}
#if 0
CEllipse2D operator *(const CMat2D &transform, const CEllipse2D &ellipse)
{
	CEllipse2D c(ellipse);
	c.transform(transform);
	return c;
}
#endif
#endif
} //end GEO
}  //end SMF
