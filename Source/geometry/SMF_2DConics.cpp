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
#include "geometry/SMF_2DConics.h"
#include "geometry/SMF_2DCircleSegment.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{
//===============CCurveBaseic===============
     /** Consider the ray as a plane and find
          intersection point with a line. Plane is defined by origin
          and normal vectors and line is defined by origin and
          direction vectors.

          \return scale factor of the line direction vector from line origin.

      */
static float pl_ln_intersect_scale(const CRay &ray, const CRay &line)
{
      return (ray.pos * ray.dir - ray.dir * line.pos) / (line.dir * ray.dir);
}



CCurveBase::~CCurveBase()
    {
    }

   // Default curve/ray intersection iterative method

    bool CCurveBase::intersect(CVec3D &point, const CRay &ray) const
    {
      CRay p;

      // initial intersection with z=0 plane
      {
        float  s = ray.dir.z;

        if (s == 0)
          return false;

        float  a = -ray.pos.z / s;

        if (a < 0)
          return false;

        p.pos = ray.pos + ray.dir * a;
      }

      unsigned int n = 32;      // avoid infinite loop

      while (n--)
        {
          float new_sag = sagitta(p.pos.project_xy()); //get saggita at point xy
          float old_sag = p.pos.z;

          // project previous intersection point on curve
          p.pos.z = new_sag;

          // stop if close enough
          if (fabs(old_sag - new_sag) < 1e-10)
            break;

          // get curve tangeante plane at intersection point
          normal(p.dir, p.pos);

          // intersect again with new tangeante plane
          float a = pl_ln_intersect_scale(p,ray);

          if (a < 0)
            return false;

          p.pos = ray.pos + ray.dir * a;
        }

      point = p.pos;

      return true;
    }

    // Default curve derivative use numerical differentiation

    struct curve_params_s
    {
      const CCurveBase *c;
      float x, y;
    };

    static float deriv_function_sagitta_x(float x, void *params)
    {
      struct curve_params_s *p = (struct curve_params_s*)params;

      return p->c->sagitta(CVec2D(x, p->y)); //get saggita at point xy
    }

    static float deriv_function_sagitta_y(float y, void *params)
    {
      struct curve_params_s *p = (struct curve_params_s*)params;

      return p->c->sagitta(CVec2D(p->x, y)); //get saggita at point xy
    }

    void CCurveBase::derivative(const CVec2D & xy, CVec2D & dxdy) const
    {
      float abserr;
      struct curve_params_s params;
      CNumDEqu deriv_function;

      deriv_function.params = &params;

      params.c = this;
      params.x = xy.x;
      params.y = xy.y;
	  deriv_function.function = deriv_function_sagitta_x;
      deriv_function.deriv_central( xy.x, 1e-6, &dxdy.x, &abserr);

      deriv_function.function = deriv_function_sagitta_y;
      deriv_function.deriv_central( xy.y, 1e-6, &dxdy.y, &abserr);
    }

    void CCurveBase::normal(CVec3D &normal, const CVec3D &point) const
    {
      CVec2D d;

      derivative(point.project_xy(), d);

      normal = CVec3D(d.x, d.y, -1.0);
      normal.normalize();
    }

	//=================CRotatConics=======================================

	   CRotationalCurv::CRotationalCurv()
    {
      deriv_function.function = deriv_function_sagitta;
      deriv_function.params = this;
    }

    void CRotationalCurv::normal(CVec3D &normal, const CVec3D &point) const
    {
      const float r = CMath::sqrt(CMath::square(point.x) + CMath::square(point.y));

      if (r == 0)
        normal = CVec3D(0, 0, -1);
      else
        {
          const float p = derivative(r);

          normal = CVec3D(point.x * p / r, point.y * p / r, -1.0);
          normal.normalize(); // FIXME simplify ?
        }
    }

    void CRotationalCurv::derivative(const CVec2D & xy, CVec2D & dxdy) const
    {
      const float r = xy.getLenght();

      if (r == 0)
        {
          dxdy.x = dxdy.y = 0.0;
          return;
        }

      const float p = derivative(r);

      dxdy = xy * (p / r);
    }

    float CRotationalCurv::deriv_function_sagitta(float x, void *params)
    {
      CRotationalCurv *c = static_cast<CRotationalCurv *>(params);

      return c->sagitta(x); //get saggita at distance x
    }

    float CRotationalCurv::derivative(float r) const
    {
      float result, abserr;
	  CNumDEqu func;
	  func=deriv_function;
	  func.deriv_central( r, 1e-4, &result, &abserr);
      return result;
    }


	float CRotationalCurv::sagitta(const CVec2D & xy) const
    {
      return sagitta(xy.getLenght()); //get saggita at point xy
    }

//===========CCurveRoc=================
CCurveRoc::CCurveRoc(float roc)
: _roc(roc)
{
}

void CCurveRoc::set_roc(float roc)
{
      _roc = roc;
}

float CCurveRoc::get_roc() const
{
      return _roc;
}

//================CConicBase==================
static float* getBuffer(int count){
	SMF_ASSERT (mem_isInitialized())
float *buffer= (float*) mem_Alloc(count*4);
return buffer;
}
static void freeBuffer(float *buffer){
	mem_Free(buffer);
}

static const char *conicName[] =
{
  "invalid conic",
  "real ellipse",
  "real circle",
  "imaginary ellipse",
  "imaginary circle",
  "hyperbola",
  "parabola",
  "real intersecting lines",
  "complex intersecting lines",
  "real parallel lines",
  "complex parallel lines",
  "coincident lines"
};


CMyString CConicBase::real_type() const { return conicName[(int)type_]; }

CConicBase::conicType_t CConicBase::type_by_name(CMyString const& name)
{
  for (int i = (int)no_type; i < num_conic_types; i++)
    if (name == conicName[i])
      return (conicType_t)i;
  return no_type; // should never reach this point
}


CMyString CConicBase::type_by_number(conicType_t type)
{
  if (type <= 0 || type >= num_conic_types) return conicName[no_type];
  return conicName[type];
}

CConicBase::CConicBase(float roc, float sc):
CCurveRoc(roc), _sh(sc + 1)
{
}

float CConicBase::get_eccentricity() const
{
      return CMath::sqrt(- _sh + 1.0);
}

float CConicBase::get_schwarzschild() const
{
      return _sh - 1.0;
}


float CConicBase::fit_roc(const CRotationalCurv &c, float radius, unsigned int count)
    {

      float *X, *Y;
	  X=getBuffer(count);
	  Y=getBuffer(count);
      float step = radius / (float)count;
      float y = step / 2.0;
      float c1, cov11, chisq;

      if (_sh != 0.0)
        {
          for (unsigned int i = 0; i < count; i++)
            {
              float x = c.sagitta(y); //get saggita at distance y

              Y[i] = CMath::square(y) + _sh * CMath::square(x);
              X[i] = 2.0 * _sh * x;

              y += step;
            }

          CLinRegression::fitMul(X, 1, Y, 1, count, &c1, &cov11, &chisq);

          _roc = c1 * _sh;
        }
      else
        { // Parabola special case
          for (unsigned int i = 0; i < count; i++)
            {
              float x = c.sagitta(y); //get saggita at distance y

              Y[i] = CMath::square(y);
              X[i] = 4.0 * x;

              y += step;
            }

          CLinRegression::fitMul(X, 1, Y, 1, count, &c1, &cov11, &chisq);

          _roc = 2.0 * c1;
        }
	  freeBuffer(X);
	  freeBuffer(Y);
      return CMath::sqrt(chisq / count); // FIXME bad rms error
    }

//========CConic================

void CConic::set_eccentricity(float e)
{
      _sh = - CMath::square(e) + 1.0;
}

void CConic::set_schwarzschild(float sc)
{
      _sh = sc + 1.0;
}

CConic::CConic(float roc, float sc):
CConicBase(roc, sc)
{
}

float CConic::derivative(float r) const
{
      // conical section differentiate (computed with Maxima)

      const float s2 = _sh * CMath::square(r);
      const float s3 = CMath::sqrt(1 - s2 / CMath::square(_roc));
      const float s4 = 2.0/(_roc * (s3+1)) + s2/(CMath::square(_roc) * _roc * s3 * CMath::square(s3 + 1));

      return r * s4;
    }

float CConic::sagitta(float r) const
    {
	//set radius of curvature on point with distance r

      return CMath::square(r) / (_roc * (CMath::sqrt( 1 - (_sh * CMath::square(r)) / CMath::square(_roc)) + 1));
    }

bool CConic::intersect(CVec3D &point, const CRay &ray) const
    {
      const float      ax = ray.pos.x;
      const float      ay = ray.pos.y;
      const float      az = ray.pos.z;
      const float      bx = ray.dir.x;
      const float      by = ray.dir.y;
      const float      bz = ray.dir.z;

      /*
        find intersection point between conical section and ray,
        Telescope optics, page 266
      */
      float a = (_sh * CMath::square(bz) + CMath::square(by) + CMath::square(bx));
      float b = ((_sh * bz * az + by * ay + bx * ax) / _roc - bz) * 2.0;
      float c = (_sh * CMath::square(az) + CMath::square(ay) + CMath::square(ax)) / _roc - 2.0 * az;

      float t;

      if (a == 0)
        {
          t = -c / b;
        }
      else
        {
          float d = CMath::square(b) - 4.0 * a * c / _roc;

          if (d < 0)
            return false;               // no intersection

          float s = CMath::sqrt(d);

          if (a * bz < 0)
            s = -s;

          if (_sh < 0)
            s = -s;

          t = (2 * c) / (s - b);
        }

      if (t <= 0)               // ignore intersection if before ray origin
        return false;

      point = ray.pos + ray.dir * t;

      return true;
    }


    float CConic::fit(const CRotationalCurv &c, float radius, unsigned int count)
    {
      float *X, *Y;
	  X=getBuffer(count);
	  Y=getBuffer(count);

      float step = radius / (float)count;
      float y = step / 2.0;

      for (unsigned int i = 0; i < count; i++)
        {
          float x = c.sagitta(y);

          Y[i] = CMath::square(y) / x;
          X[i] = x;

          y += step;
        }

      float c0, c1, cov00, cov01, cov11, chisq;

      CLinRegression::fit_linear(X, 1, Y, 1, count,
                     &c0, &c1, &cov00, &cov01, &cov11,
                     &chisq);

      _sh = -c1;
      _roc = c0 / 2.0;
	  freeBuffer(X);
	  freeBuffer(Y);

      return CMath::sqrt(chisq / count); // FIXME bad rms error
    }



void CConic2D::set_type_from_equation()
{
  float A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;

  /* determinant, subdeterminants and trace values */
  float det = A*(C*F - E*E) - B*(B*F - D*E) + D*(B*E - C*D); // determinant
  float J = A*C - B*B;  // upper 2x2 determinant
  float K = (C*F - E*E) + (A*F - D*D); // sum of two other 2x2 determinants
  float I = A + C; // trace of upper 2x2

  if (det != 0) {
    if (J > 0) {
      if (det*I < 0) {
        if (A==C && B==0)      type_ = real_circle;
        else                   type_ = real_ellipse;
      }
      else {
        if (A==C && B==0)      type_ = imaginary_circle;
        else                   type_ = imaginary_ellipse;
      }
    }
    else if (J < 0)            type_ = hyperbola;
    else /* J == 0 */          type_ = parabola;
  }
  else {    // limiting cases
    if (J < 0)                 type_ = real_intersecting_lines;
    else if (J > 0)            type_ = complex_intersecting_lines;
    else /* J == 0 */ {
      if ( A == 0 && B == 0 && C == 0 ) { // line at infinity is component
        if ( D !=0 || E != 0 ) type_ = real_intersecting_lines;
        else if (F != 0)       type_ = coincident_lines; // 2x w=0
        else                   type_ = no_type; // all coefficients are 0
      }
      else if (K < 0)          type_ = real_parallel_lines;
      else if (K > 0)          type_ = complex_parallel_lines;
      else                     type_ = coincident_lines;
    }
  }
}

CConic2D::CConic2D(float ta, float tb, float tc, float td, float te, float tf)
  :   a_(ta), b_(tb), c_(tc), d_(td), e_(te), f_(tf)
{
  set_type_from_equation();
  float eccentr, roc;
  eccentr= getEccentricity();
  set_eccentricity(eccentr);

}

CConic2D::CConic2D(float const co[])
  :  a_(co[0]), b_(co[1]), c_(co[2]), d_(co[3]), e_(co[4]), f_(co[5])
{
  set_type_from_equation();
  float eccentr;
  eccentr= getEccentricity();
  set_eccentricity(eccentr);
}
 CConic2D::CConic2D(CCircle2D const &circ){
	 float a=1;
     float b=0;
	 float c=1;
	 float d=-circ.center.x_;
	 float e=-circ.center.y_;
	 float f=(d*d)+(e*e)-circ.r;
	 setFromPol(a,b,c,2*d,2*e,f);
 }

int CConic2D::getPoint(float x,CPoint2D &p1,CPoint2D &p2)const{
p1.toInfinite();
p2.toInfinite();
//cy^2+(bx+e)y+(f+ax^2+dx)=0
float bx_e= (b()*x)+e();
float free = (f()+(a()*CMath::square(x))+(d()*x));
float root1,root2;
int num = CPolynomial::solveQuadratic(c(),bx_e,free,root1,root2);
if (num==0) return num;
else if(num>0) p1.set(x,root1);
else if (num >1) p2.set(x,root2);
return num;
}

bool CConic2D::operator==(CConic2D const& that) const
{
  if ( type() != that.type() ) return false;
  return   a()*that.b() == b()*that.a()
        && a()*that.c() == c()*that.a()
        && a()*that.d() == d()*that.a()
        && a()*that.e() == e()*that.a()
        && a()*that.f() == f()*that.a()
        && b()*that.c() == c()*that.b()
        && b()*that.d() == d()*that.b()
        && b()*that.e() == e()*that.b()
        && b()*that.f() == f()*that.b()
        && c()*that.d() == d()*that.c()
        && c()*that.e() == e()*that.c()
        && c()*that.f() == f()*that.c()
        && d()*that.e() == e()*that.d()
        && d()*that.f() == f()*that.d()
        && e()*that.f() == f()*that.e();
}

//-------------------------------------------------------------
// set values
void CConic2D::setFromPol(float ta, float tb, float tc, float td, float te, float tf)
{
  a_ = ta; b_ = tb; c_ = tc; d_ = td; e_ = te; f_ = tf;
  set_type_from_equation();
   float eccentr, roc;
  eccentr= getEccentricity();
  set_eccentricity(eccentr);

}

void CConic2D::setFromGen(float ta, float tb, float tc, float td, float te, float tf)
{
  a_ = ta; b_ = tb*2; c_ = tc; d_ = td*2; e_ = te*2; f_ = tf;
  set_type_from_equation();
   float eccentr, roc;
  eccentr= getEccentricity();
  set_eccentricity(eccentr);

}


//-------------------------------------------------------------


float CConic2D::getDiscriminator() const{
  float A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;

	float disc = A*A-(4*A*C);
return disc;
}
bool CConic2D::isDegenerate(float epsilon) const{
float A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
CMat3D mat(A,B/2,D/2,B/2,C,E/2,D/2,E/2,F);
float det= mat.determinant();
return det <= epsilon;
}

bool CConic2D::isCentral() const
{
  return type_ == real_ellipse|| type_ == imaginary_ellipse|| type_ == hyperbola
      || type_ == real_circle || type_ == imaginary_circle
      || type_ == real_intersecting_lines|| type_ == complex_intersecting_lines;
}
//--------------------------------------------------------------------------------
void CConic2D::translate_by(float x, float y)
{
  d_ += 2*a_*x + b_*y;
  f_ += c_ * y*y - a_ * x*x + d_ * x + e_ * y;
  e_ += 2*c_*y + b_*x;
  // This does not change the type, so no need to run set_type_from_equation()
}
CConic2D CConic2D::dual_conic() const
{
  float A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
  return CConic2D(E*E-C*F, 2*(B*F-D*E), D*D-A*F, 2*(C*D-B*E), 2*(A*E-B*D), B*B-A*C);
}

float CConic2D::curvature_at(CPoint2D const& p) const
{
  // Shorthands
  const float &a_xx = a_;
  const float &a_xy = b_;
  const float &a_yy = c_;
  const float &a_xw = d_;
  const float &a_yw = e_;

  const float x = p.x();
  const float y = p.y();

  float f_x  = 2*a_xx*x + a_xy*y + a_xw;
  float f_y  = 2*a_yy*y + a_xy*x + a_yw;
  float f_xy = a_xy;
  float f_xx = 2*a_xx;
  float f_yy = 2*a_yy;

  float f_x_2 = f_x*f_x;
  float f_y_2 = f_y*f_y;
  float denom = f_x_2 + f_y_2;
  denom = CMath::sqrt(denom*denom*denom);

  // Divergent of the unit normal grad f/|grad f|
  return (f_xx*f_y_2 - 2*f_x*f_y*f_xy + f_yy*f_x_2) / denom;
}

float CConic2D::getEccentricity()const{
if (type_==parabola) return 1;
if (type_== imaginary_ellipse || (type_==hyperbola && isDegenerate())||(type_==real_ellipse && isDegenerate()) ) {
return CMath::NAN_FLOAT;
}else {
float e,num,denom,B2,A_C2,n;
CMat3D mat=eccentMat();
float determ = mat.determinant();
determ <= 0 ? n=1 : n = -1;
B2=CMath::square(B());
A_C2=CMath::square(A()-C());
num = CMath::sqrt(A_C2+B2);
denom =(n* (A()+C()))+num;
e = CMath::sqrt(2*num/denom);
return e;
}

}

bool CConic2D::ellipse_geometry(float& xc, float& yc, float& major_axis_length,
                 float& minor_axis_length, float& angle_in_radians)
{
  if (type_!=real_ellipse && type_ != real_circle)
    return false;

  // Cast to float and half the non-diagonal (non-quadratic) entries B, D, E.
  float A = static_cast<float>(a_), B = static_cast<float>(b_)*0.5,
         C = static_cast<float>(c_), D = static_cast<float>(d_)*0.5,
         F = static_cast<float>(f_), E = static_cast<float>(e_)*0.5;
  if (A < 0)
    A=-A, B=-B, C=-C, D=-D, E=-E, F=-F;

  float det = A*(C*F - E*E) - B*(B*F- D*E) + D*(B*E-C*D);
  float D2 =  A*C - B*B;
  xc = (E*B - C*D)/D2;
  yc = (D*B - A*E)/D2;

  float trace = A + C;
  float disc = CMath::sqrt(trace*trace - 4.0*D2);
  float cmaj = (trace+disc)*D2/(2*det); if (cmaj < 0) cmaj = -cmaj;
  float cmin = (trace-disc)*D2/(2*det); if (cmin < 0) cmin = -cmin;
  minor_axis_length = 1.0/CMath::sqrt(cmaj>cmin?cmaj:cmin);
  major_axis_length = 1.0/CMath::sqrt(cmaj>cmin?cmin:cmaj);

  // find the angle that diagonalizes the upper 2x2 sub-matrix
  angle_in_radians  = -0.5 * CMath::atan64(2*B, C-A);
  //                  ^
  // and return the negative of this angle
  return true;
}


bool CConic2D::contains(CPoint2D const& p) const
{
  return p.x()*p.x()*a_+p.x()*p.y()*b_+p.y()*p.y()*c_+p.x()*d_+p.y()*e_+f_ <= 0;
}
bool CConic2D::edgeContains(CPoint2D const& p ,float epsilon) const
{
  float num =p.x()*p.x()*a_+p.x()*p.y()*b_+p.y()*p.y()*c_+p.x()*d_+p.y()*e_+f_ ;
 return num >= -epsilon && num <= epsilon;
}

bool CConic2D::intersection(const CConic2D &conic1, const CConic2D &conic2, float numIntersPoints,  CPoint2D &p1,  CPoint2D &p2, CPoint2D &p3, CPoint2D &p4,int prec){
	p1.toInfinite();
	p2.toInfinite();
	p3.toInfinite();
	p4.toInfinite();

	CVec2D pt1,pt2,pt3,pt4;
	CPolyConSystem sist;
	sist.setEqu1(conic1.a(),conic1.b(),conic1.c(),conic1.d(),conic1.e(),conic1.f());
	sist.setEqu2(conic2.a(),conic2.b(),conic2.c(),conic2.d(),conic2.e(),conic2.f());
	int soluc=sist.getSol(pt1,pt2,pt3,pt4,prec);
	numIntersPoints=soluc;
	if(soluc>0){
		p1=pt1;
		p2=pt2;
		p3=pt3;
		p4=pt4;
		return true;
	}
	else return false;
}

CMat3D CConic2D::characMat()const{
	CMat3D conMat;
  /* compute the transformation using matrix muliplication */
   conMat[0][0] = this->A();
   conMat[1][0] = conMat[0][1] = this->B();
   conMat[1][1] = this->C();
   conMat[2][0] = conMat[0][2] = this->D();
   conMat[2][1] = conMat[1][2] = this->E();
   conMat[2][2] = this->F();
   return conMat;
}
CMat3D CConic2D::eccentMat()const{
	CMat3D conMat;
  /* compute the transformation using matrix muliplication */
   conMat[0][0] = this->A();
   conMat[1][0] = conMat[0][1] = this->B()/2;
   conMat[1][1] = this->C();
   conMat[2][0] = conMat[0][2] = this->D()/2;
   conMat[2][1] = conMat[1][2] = this->E()/2;
   conMat[2][2] = this->F();
   return conMat;
}


//: Write "<CConic2D aX^2+bXY+cY^2+dX+eY+f=0>" to stream

ostream& CConic2D::operator<<(ostream& s)
{
  CConic2D const& co=*this;
  s << "<Cconic2D ";
  if (co.a() == 1) s << "X^2";
  else if (co.a() == -1) s << "-X^2";
  else if (co.a() != 0) s << co.a() << "X^2";
  if (co.b() > 0) s << '+';
  if (co.b() == 1) s << "XY";
  else if (co.b() == -1) s << "-XY";
  else if (co.b() != 0) s << co.b() << "XY";
  if (co.c() > 0) s << '+';
  if (co.c() == 1) s << "Y^2";
  else if (co.c() == -1) s << "-Y^2";
  else if (co.c() != 0) s << co.c() << "Y^2";
  if (co.d() > 0) s << '+';
  if (co.d() == 1) s << "X";
  else if (co.d() == -1) s << "-X";
  else if (co.d() != 0) s << co.d() << "X";
  if (co.e() > 0) s << '+';
  if (co.e() == 1) s << "Y";
  else if (co.e() == -1) s << "-Y";
  else if (co.e() != 0) s << co.e() << "Y";
  if (co.f() > 0) s << '+';
  if (co.f() == 1) s ;
  else if (co.f() == -1) s ;
  else if (co.f() != 0) s << co.f() ;
  return s << "=0 " << co.real_type() << "> ";
}

//: Read a b c d e f from stream
istream& CConic2D::operator>>(istream& is)
{
  CConic2D& co=*this;
  float ta, tb, tc, td, te, tf; is >> ta >> tb >> tc >> td >> te >> tf;
  co.setFromPol(ta,tb,tc,td,te,tf); return is;
}

//============ELipse intersection=================

/* transformation matrix type */
struct TMat
   {
	   TMat(){};
	    TMat(float A,float B,float C,float D,float M,float N):a(A),b(B),c(C),d(D),m(M),n(N){};
	   //remember CMat3D is column major
	CMat3D toMat3(){CMat3D trans(a,b,0,c,d,0,m,n,1); return trans;};

   float a,b,c,d;		 /* tranformation coefficients */
   float m,n;			 /* translation coefficients   */
   } ;


/* prototypes */

/* Functions with names beginning with a O:
 *
 * 1. Ignore CEllipse2D center coordinates.
 *
 * 2. Ignore any translation components in tranformation matrices.
 *
 * 3. Assume that conic coeficients were generated by "O-name" functions.
 *
 * This keeps the computations relatively simple. The CEllipse2D centers
 * can then be transformed separately as points.
 */

/* OTransformConic - transform conic coefficients about the origin */
void OTransformConic(CConic2D *ConicP,TMat *TMatP);

/* OGenMyEllipseCoefs - Generate conic coefficients of an CEllipse2D */
void OGenMyEllipseCoefs(CEllipse2D *MyEllipseP,CConic2D *ConicP);

/* OGenMyEllipseGeom - Generates CEllipse2D geometry from conic coefficients */
void OGenMyEllipseGeom(CConic2D *ConicP,CEllipse2D *MyEllipseP);

/* TransformPoint - transform a point using a tranformation matrix */
void TransformPoint(CPoint2D *PointP,TMat *TMatP);

/* TransformMyEllipse - transform an CEllipse2D using a tranformation matrix */
void TransformMyEllipse(CEllipse2D *MyEllipseP,TMat *TMatP);



/* identity matrix */
static const TMat IdentMat(1.0,0.0,0.0,1.0,0.0,0.0);

/* Transformation matrix routines */

/* translate a matrix by m,n */
void TranslateMat(TMat *Mat,float m,float n)
   {
   Mat->m += m;
   Mat->n += n;
   }

/* rotate a matrix by Phi */
void RotateMat(TMat *Mat,float Phi)
   {
   float SinPhi=sin(Phi);
   float CosPhi=cos(Phi);
   TMat temp=*Mat;		/* temporary copy of Mat */

   /* These are just the matrix operations written out long hand */
   Mat->a = temp.a*CosPhi - temp.b*SinPhi;
   Mat->b = temp.b*CosPhi + temp.a*SinPhi;
   Mat->c = temp.c*CosPhi - temp.d*SinPhi;
   Mat->d = temp.d*CosPhi + temp.c*SinPhi;
   Mat->m = temp.m*CosPhi - temp.n*SinPhi;
   Mat->n = temp.n*CosPhi + temp.m*SinPhi;
   }

/* scale a matrix by sx, sy */
void ScaleMat(TMat *Mat,float sx,float sy)
   {
   Mat->a *= sx;
   Mat->b *= sy;
   Mat->c *= sx;
   Mat->d *= sy;
   Mat->m *= sx;
   Mat->n *= sy;
   }

/* TransformPoint - transform a point using a tranformation matrix */
void TransformPoint(CPoint2D *PointP,TMat *TMatP)
   {
   CMat3D trans(TMatP->a,TMatP->b,0,TMatP->c,TMatP->d,0,TMatP->m,TMatP->n,1);
   *PointP=PointP->transform(TMatP->toMat3());
   }

/* CConic2D routines */

/* near zero test */
#define EPSILON 1e-9
#define isZero(x) (x > -EPSILON && x < EPSILON)

/* GenMyEllipseCoefs - Generate conic coefficients of an CEllipse2D */
static void GenMyEllipseCoefs(CEllipse2D *elip,CConic2D *M)
   {
   float sqr_r1,sqr_r2;
   float sint,cost,sin2t,sqr_sint,sqr_cost;
   float cenx,ceny,sqr_cenx,sqr_ceny,invsqr_r1,invsqr_r2;

   /* common coeficients */
   sqr_r1 = elip->a();
   sqr_r2 = elip->b();
   sqr_r1 *= sqr_r1;
   sqr_r2 *= sqr_r2;
   sint = sin(elip->Phi());
   cost = cos(elip->Phi());
   sin2t = 2.0*sint*cost;
   sqr_sint = sint*sint;
   sqr_cost = cost*cost;
   cenx = elip->center.x_;
   sqr_cenx = cenx*cenx;
   ceny = elip->center.y_;
   sqr_ceny = ceny*ceny;
   invsqr_r1 = 1.0/sqr_r1;
   invsqr_r2 = 1.0/sqr_r2;
   float A,B,C,D,E,F;
   /* Compute the coefficients. These formulae are the transformations
      on the unit circle written out long hand */
   A = sqr_cost/sqr_r1 + sqr_sint/sqr_r2;
   B = (sqr_r2-sqr_r1)*sin2t/(2.0*sqr_r1*sqr_r2);
   C = sqr_cost/sqr_r2 + sqr_sint/sqr_r1;
   D = -ceny*B-cenx*A;
   E = -cenx*B-ceny*C;
   F = -1.0 + (sqr_cenx + sqr_ceny)*(invsqr_r1 + invsqr_r2)/2.0 +
      (sqr_cost - sqr_sint)*(sqr_cenx - sqr_ceny)*(invsqr_r1 - invsqr_r2)/2.0 +
      cenx*ceny*(invsqr_r1 - invsqr_r2)*sin2t;
   M->setFromGen(A,B,C,D,E,F);
   }

/* Compute the transformation which turns an CEllipse2D into a circle */
void Elp2Cir(CEllipse2D *Elp,TMat *CirMat)
   {
   /* Start with identity matrix */
   *CirMat = IdentMat;
   /* translate to origin */
   TranslateMat(CirMat,-Elp->center.x_,-Elp->center.y_);
   /* rotate into standard position */
   RotateMat(CirMat,-Elp->Phi());
   /* scale into a circle. */
   ScaleMat(CirMat,1.0/Elp->a(),1.0/Elp->b());
   }

/* Compute the inverse of the transformation
   which turns an CEllipse2D into a circle */
void InvElp2Cir(CEllipse2D *Elp,TMat *InvMat)
   {
   /* Start with identity matrix */
   *InvMat = IdentMat;
   /* scale back into an CEllipse2D. */
   ScaleMat(InvMat,Elp->a(),Elp->b());
   /* rotate */
   RotateMat(InvMat,Elp->Phi());
   /* translate from origin */
   TranslateMat(InvMat,Elp->center.x_,Elp->center.y_);
   }

/* OTransformConic - transform conic coefficients about the origin	 */
/* This routine ignores the translation components of *TMatP and
   assumes the conic is "centered" at the origin (i.e. D,E=0, F=-1)	 */
/* The computations are just the matrix operations written out long hand */
/* This code assumes that the transformation is not degenerate		 */
void OTransformConic(CConic2D *ConicP,TMat *TMatP)
   {
   float A,B,C,D,E,F,Denom;

   /* common denominator for transformed cooefficients */
   Denom = TMatP->a*TMatP->d - TMatP->b*TMatP->c;
   Denom *= Denom;

   A = (ConicP->C()*TMatP->b*TMatP->b - 2.0*ConicP->B()*TMatP->b*TMatP->d +
      ConicP->A()*TMatP->d*TMatP->d)/Denom;

   B = (-ConicP->C()*TMatP->a*TMatP->b + ConicP->B()*TMatP->b*TMatP->c +
      ConicP->B()*TMatP->a*TMatP->d - ConicP->A()*TMatP->c*TMatP->d)/Denom;

   C = (ConicP->C()*TMatP->a*TMatP->a - 2.0*ConicP->B()*TMatP->a*TMatP->c +
      ConicP->A()*TMatP->c*TMatP->c)/Denom;
   D = ConicP->D();
   E = ConicP->E();
   F = ConicP->F();
   ConicP->setFromGen(A,B,C,D,E,F);
   }

/* OGenMyEllipseCoefs - Generate conic coefficients of an CEllipse2D */
/* The CEllipse2D is assumed to be centered at the origin. */
void OGenMyEllipseCoefs(CEllipse2D *MyEllipseP,CConic2D *ConicP)
   {
   float SinPhi = sin(MyEllipseP->Phi()); /* sine of CEllipse2D rotation   */
   float CosPhi = cos(MyEllipseP->Phi()); /* cosine of CEllipse2D rotation */
   float SqSinPhi = SinPhi*SinPhi;    /* square of sin(phi)	     */
   float SqCosPhi = CosPhi*CosPhi;    /* square of cos(phi)	     */
   float SqMaxRad = MyEllipseP->a()*MyEllipseP->a();
   float SqMinRad = MyEllipseP->b()*MyEllipseP->b();
   float A,B,C,D,E,F;
   /* compute coefficients for the CEllipse2D in standard position */
   A = SqCosPhi/SqMaxRad + SqSinPhi/SqMinRad;
   B = (1.0/SqMaxRad - 1.0/SqMinRad)*SinPhi*CosPhi;
   C = SqCosPhi/SqMinRad + SqSinPhi/SqMaxRad;
   D = E = 0.0;
   F = -1.0;
   ConicP->setFromGen(A,B,C,D,E,F);
   }

/* OGenMyEllipseGeom - Generates CEllipse2D geometry from conic coefficients */
/* This routine assumes the conic coefficients D=E=0, F=-1 */
void OGenMyEllipseGeom(CConic2D *ConicP,CEllipse2D *MyEllipseP)
   {
   float Numer,Denom,Temp;
   TMat DiagTransform;		 /* transform diagonalization */
   CConic2D ConicT = *ConicP;	 /* temporary copy of conic coefficients */
   float SinPhi,CosPhi;

   /* compute new CEllipse2D rotation */
   Numer = ConicT.B() + ConicT.B();
   Denom = ConicT.A() - ConicT.C();
   /* Phi = 1/2 atan(Numer/Denom) = 1/2 (pi/2 - atan(Denom/Numer)
      We use the form that keeps the argument to atan between -1 and 1 */

   float Phi = 0.5*(fabs(Numer) < fabs(Denom)?
      atan(Numer/Denom) :  CMath::HALF_PI-CMath::atan(Denom/Numer));

   /* diagonalize the conic */
   SinPhi = CMath::sin(MyEllipseP->Phi());
   CosPhi = CMath::cos(MyEllipseP->Phi());
   DiagTransform.a = CosPhi;	 /* rotate by -Phi */
   DiagTransform.b = -SinPhi;
   DiagTransform.c = SinPhi;
   DiagTransform.d = CosPhi;
   DiagTransform.m = DiagTransform.n = 0.0;
   OTransformConic(&ConicT,&DiagTransform);

   /* compute new radii from diagonalized coefficients */
   float NewMaxRadio = 1.0/CMath::sqrt(ConicT.A());
   float NewMinRadio = 1.0/CMath::sqrt(ConicT.C());

   /* be sure a() >= b() */
   if (NewMaxRadio < NewMinRadio)
      {
      Temp = MyEllipseP->a();	 /* exchange the radii */
      NewMaxRadio = NewMinRadio;
      NewMinRadio = Temp;
      Phi += CMath::HALF_PI;	 /* adjust the rotation */
      }
   MyEllipseP->set(CPoint2D(0,0),NewMaxRadio,NewMaxRadio,Phi);
   }

/* TransformMyEllipse - transform an CEllipse2D using a tranformation matrix */
void TransformMyEllipse(CEllipse2D *MyEllipseP,TMat *TMatP)
   {
   CConic2D MyEllipseCoefs;

   /* generate the CEllipse2D coefficients (using center=origin) */
   OGenMyEllipseCoefs(MyEllipseP,&MyEllipseCoefs);

   /* transform the coefficients */
   OTransformConic(&MyEllipseCoefs,TMatP);

   /* turn the transformed coefficients back into geometry */
   OGenMyEllipseGeom(&MyEllipseCoefs,MyEllipseP);

   /* translate the center */
   CMat3D trans(TMatP->a,TMatP->b,0,TMatP->c,TMatP->d,0,TMatP->m,TMatP->n,1);
   MyEllipseP->center=MyEllipseP->center.transform(trans);
   }


/* transform a conic by a tranformation matrix */
void TransformConic(CConic2D *ConicP,TMat *TMatP)
   {
   CMat3D invMat, conMat,tranInvMat,result1,result2;
   float D;
   int i,j;

   /* Compute M' = Inv(TMat).M.transpose(Inv(TMat))

   /* compute the transformation using matrix muliplication */
   //get the characteristic matriz
   conMat = ConicP->characMat();

   /* inverse transformation */
   D = TMatP->a*TMatP->d - TMatP->b*TMatP->c;
   invMat[0][0] = TMatP->d/D;
   invMat[1][0] = -TMatP->b/D;
   invMat[2][0] = 0.0;
   invMat[0][1] = -TMatP->c/D;
   invMat[1][1] = TMatP->a/D;
   invMat[2][1] = 0.0;
   invMat[0][2] = (TMatP->c*TMatP->n - TMatP->d*TMatP->m)/D;
   invMat[1][2] = (TMatP->b*TMatP->m - TMatP->a*TMatP->n)/D;
   invMat[2][2] = 1.0;

   /* compute transpose */
   tranInvMat = invMat.transpose();

   /* multiply the matrices */
   result1 = invMat*conMat;
   result2 = result1*tranInvMat;

   float AA,BA,CA,DA,EA,FA;
   AA = result2[0][0];	       /* return to conic form */
   BA = result2[1][0];
   CA = result2[1][1];
   DA = result2[2][0];
   EA = result2[2][1];
   FA = result2[2][2];

   ConicP->setFromGen(AA,BA,CA,DA,EA,FA);
   }

/* Compute the intersection of a circle and a line */
/* See Graphic Gems volume 1, page 5 for a description of this algorithm */
int IntCirLine(CPoint2D *IntPts,CCircle2D *Cir,CLineSegment2D *Ln)
   {
   CPoint2D G,V;
   float a,b,c,d,t,sqrt_d;

   G.x_ = Ln->point1_.x_ - Cir->center.x_;     	/* G = Ln->point1_ - Cir->center */
   G.y_ = Ln->point1_.y_ - Cir->center.y_;
   V.x_ = Ln->point2_.x_ - Ln->point1_.x_;           /* V = Ln->point2_ - Ln->point1_ */
   V.y_ = Ln->point2_.y_ - Ln->point1_.y_;
   a = V.x_*V.x_ + V.y_*V.y_;		/* a = V.V */
   b = V.x_*G.x_ + V.y_*G.y_;b += b; 	/* b = 2(V.G) */
   c = (G.x_*G.x_ + G.y_*G.y_) -    	/* c = G.G + CCircle2D->Radius^2 */
      Cir->r*Cir->r;
   d = b*b - 4.0*a*c;			/* discriminant */

   if (d <= 0.0)
      return 0;				/* no intersections */

   sqrt_d = CMath::sqrt(d);
   t = (-b + sqrt_d)/(a + a);           /* t = (-b +/- CMath::sqrt(d))/2a */
   IntPts[0].x_ = Ln->point1_.x_ + t*V.x_;      /* Pt = Ln->point1_ + t V */
   IntPts[0].y_ = Ln->point1_.y_ + t*V.y_;
   t = (-b - sqrt_d)/(a + a);
   IntPts[1].x_ = Ln->point1_.x_ + t*V.x_;
   IntPts[1].y_ = Ln->point1_.y_ + t*V.y_;
   return 2;
   }

/* compute all intersections of two MyEllipses */
/* E1 and E2 are the two MyEllipses */
/* IntPts points to an array of twelve points
    (some duplicates may be returned) */
/* The number of intersections found is returned */
/* Both MyEllipses are assumed to have non-zero radii */
int intersects(CEllipse2D &E1,CEllipse2D &E2, CPoint2D *IntPts)
   {
   TMat ElpCirMat1,ElpCirMat2,InvMat,TempMat;
   CConic2D Conic1,Conic2,Conic3,TempConic;
   float Roots[3],qRoots[2];
   CPoint2D pc(0.0,0.0);
   static CCircle2D TestCir(pc,1.0);
   CLineSegment2D TestLine[2];
   CPoint2D TestPoint;
   float PolyCoef[4];		/* coefficients of the polynomial */
   float D;			/* discriminant: B^2 - AC */
   float Phi;			/* CEllipse2D rotation */
   float m,n;			/* CEllipse2D translation */
   float Scl;			/* scaling factor */
   int NumRoots,NumLines;
   int CircleInts;		/* intersections between line and circle */
   int IntCount = 0;		/* number of intersections found */
   int i,j,k;

   /* compute the transformations which turn E1 and E2 into circles */
   Elp2Cir(&E1,&ElpCirMat1);
   Elp2Cir(&E2,&ElpCirMat2);

   /* compute the inverse transformation of ElpCirMat1 */
   InvElp2Cir(&E1,&InvMat);

   /* Compute the characteristic matrices */
   GenMyEllipseCoefs(&E1,&Conic1);
   GenMyEllipseCoefs(&E2,&Conic2);

   /* find x such that Det(Conic1 + x Conic2) = 0 */
   PolyCoef[0] = -Conic1.C()*Conic1.D()*Conic1.D() + 2.0*Conic1.B()*Conic1.D()*Conic1.E() -
      Conic1.A()*Conic1.E()*Conic1.E() - Conic1.B()*Conic1.B()*Conic1.F() +
      Conic1.A()*Conic1.C()*Conic1.F();
   PolyCoef[1] = -(Conic2.C()*Conic1.D()*Conic1.D()) -
      2.0*Conic1.C()*Conic1.D()*Conic2.D() + 2.0*Conic2.B()*Conic1.D()*Conic1.E() +
      2.0*Conic1.B()*Conic2.D()*Conic1.E() - Conic2.A()*Conic1.E()*Conic1.E() +
      2.0*Conic1.B()*Conic1.D()*Conic2.E() - 2.0*Conic1.A()*Conic1.E()*Conic2.E() -
      2.0*Conic1.B()*Conic2.B()*Conic1.F() + Conic2.A()*Conic1.C()*Conic1.F() +
      Conic1.A()*Conic2.C()*Conic1.F() - Conic1.B()*Conic1.B()*Conic2.F() +
      Conic1.A()*Conic1.C()*Conic2.F();
   PolyCoef[2] = -2.0*Conic2.C()*Conic1.D()*Conic2.D() - Conic1.C()*Conic2.D()*Conic2.D() +
      2.0*Conic2.B()*Conic2.D()*Conic1.E() + 2.0*Conic2.B()*Conic1.D()*Conic2.E() +
      2.0*Conic1.B()*Conic2.D()*Conic2.E() - 2.0*Conic2.A()*Conic1.E()*Conic2.E() -
      Conic1.A()*Conic2.E()*Conic2.E() - Conic2.B()*Conic2.B()*Conic1.F() +
      Conic2.A()*Conic2.C()*Conic1.F() - 2.0*Conic1.B()*Conic2.B()*Conic2.F() +
      Conic2.A()*Conic1.C()*Conic2.F() + Conic1.A()*Conic2.C()*Conic2.F();
   PolyCoef[3] = -Conic2.C()*Conic2.D()*Conic2.D() + 2.0*Conic2.B()*Conic2.D()*Conic2.E() -
      Conic2.A()*Conic2.E()*Conic2.E() - Conic2.B()*Conic2.B()*Conic2.F() +
      Conic2.A()*Conic2.C()*Conic2.F();
   NumRoots = SMF::MATH::CPolynomial::solveCubic(PolyCoef,Roots);

   if (NumRoots == 0)
      return 0;

   /* we try all the roots, even though it's redundant, so that we
      avoid some pathological situations */
   for (i=0;i<NumRoots;i++)
      {
      NumLines = 0;
	  float A,B,C,D,E,F;
      /* Conic3 = Conic1 + mu Conic2 */
      A = Conic1.A() + Roots[i]*Conic2.A();
      B = Conic1.B() + Roots[i]*Conic2.B();
      C = Conic1.C() + Roots[i]*Conic2.C();
      D = Conic1.D() + Roots[i]*Conic2.D();
      E = Conic1.E() + Roots[i]*Conic2.E();
      F = Conic1.F() + Roots[i]*Conic2.F();
	  Conic3.setFromGen(A,B,C,D,E,F);
      D = Conic3.B()*Conic3.B() - Conic3.A()*Conic3.C();
      if (isZero(Conic3.A()) && isZero(Conic3.B()) && isZero(Conic3.C()))
	 {
	 /* Case 1 - Single line */
	 NumLines = 1;
	 /* compute endpoints of the line, avoiding division by zero */
	 if (fabs(Conic3.D()) > fabs(Conic3.E()))
	    {
	    TestLine[0].point1_.y_ = 0.0;
	    TestLine[0].point1_.x_ = -Conic3.F()/(Conic3.D() + Conic3.D());
	    TestLine[0].point2_.y_ = 1.0;
	    TestLine[0].point2_.x_ = -(Conic3.E() + Conic3.E() + Conic3.F())/
	       (Conic3.D() + Conic3.D());
	    }
	 else
	    {
	    TestLine[0].point1_.x_ = 0.0;
	    TestLine[0].point1_.y_ = -Conic3.F()/(Conic3.E() + Conic3.E());
	    TestLine[0].point2_.x_ = 1.0;
	    TestLine[0].point2_.x_ = -(Conic3.D() + Conic3.D() + Conic3.F())/
	       (Conic3.E() + Conic3.E());
	    }
	 }
      else
	 {
	 /* use the espresion for Phi that takes atan of the
	    smallest argument */
	 Phi = (fabs(Conic3.B() + Conic3.B()) < fabs(Conic3.A()-Conic3.C())?
	    CMath::atan((Conic3.B() + Conic3.B())/(Conic3.A() - Conic3.C())) :
	    CMath::HALF_PI - CMath::atan((Conic3.A() - Conic3.C())/(Conic3.B() + Conic3.B())))/2.0;
	 if (isZero(D))
	    {
	    /* Case 2 - Parallel lines */
	    TempConic = Conic3;
	    TempMat = IdentMat;
	    RotateMat(&TempMat,-Phi);
	    TransformConic(&TempConic,&TempMat);
	    if (isZero(TempConic.C()))   /* vertical */
	       {
	       PolyCoef[0] = TempConic.F();
	       PolyCoef[1] = TempConic.D();
	       PolyCoef[2] = TempConic.A();
	       if ((NumLines=SMF::MATH::CPolynomial::solveQuadratic(PolyCoef,qRoots))!=0)
		  {
		  TestLine[0].point1_.x_ = qRoots[0];
		  TestLine[0].point1_.y_ = -1.0;
		  TestLine[0].point2_.x_ = qRoots[0];
		  TestLine[0].point2_.y_ = 1.0;
		  if (NumLines==2)
		     {
		     TestLine[1].point1_.x_ = qRoots[1];
		     TestLine[1].point1_.y_ = -1.0;
		     TestLine[1].point2_.x_ = qRoots[1];
		     TestLine[1].point2_.y_ = 1.0;
		     }
		  }
	       }
	    else		    /* horizontal */
	       {
	       PolyCoef[0] = TempConic.F();
	       PolyCoef[1] = TempConic.E();
	       PolyCoef[2] = TempConic.C();
	       if ((NumLines=SMF::MATH::CPolynomial::solveQuadratic(PolyCoef,qRoots))!=0)
		  {
		  TestLine[0].point1_.x_ = -1.0;
		  TestLine[0].point1_.y_ = qRoots[0];
		  TestLine[0].point2_.x_ = 1.0;
		  TestLine[0].point2_.y_ = qRoots[0];
		  if (NumLines==2)
		     {
		     TestLine[1].point1_.x_ = -1.0;
		     TestLine[1].point1_.y_ = qRoots[1];
		     TestLine[1].point2_.x_ = 1.0;
		     TestLine[1].point2_.y_ = qRoots[1];
		     }
		  }
	       }
	    TempMat = IdentMat;
	    RotateMat(&TempMat,Phi);

   TestLine[0].point1_.transform(TempMat.toMat3());
	    TransformPoint(&TestLine[0].point2_,&TempMat);
	    if (NumLines==2)
	       {
	       TransformPoint(&TestLine[1].point1_,&TempMat);
	       TransformPoint(&TestLine[1].point2_,&TempMat);
	       }
	    }
	 else
	    {
	    /* Case 3 - Crossing lines */
	    NumLines = 2;

	    /* translate the system so that the intersection of the lines
	       is at the origin */
	    TempConic = Conic3;
	    m = (Conic3.C()*Conic3.D() - Conic3.B()*Conic3.E())/D;
	    n = (Conic3.A()*Conic3.E() - Conic3.B()*Conic3.D())/D;
	    TempMat = IdentMat;
	    TranslateMat(&TempMat,-m,-n);
	    RotateMat(&TempMat,-Phi);
	    TransformConic(&TempConic,&TempMat);

	    /* Compute the line endpoints */
	    TestLine[0].point1_.x_ = CMath::sqrt(fabs(1.0/TempConic.A()));
	    TestLine[0].point1_.y_ = CMath::sqrt(fabs(1.0/TempConic.C()));
	    Scl = MAX(TestLine[0].point1_.x_,TestLine[0].point1_.y_);  /* adjust range */
	    TestLine[0].point1_.x_ /= Scl;
	    TestLine[0].point1_.y_ /= Scl;
	    TestLine[0].point2_.x_ = - TestLine[0].point1_.x_;
	    TestLine[0].point2_.y_ = - TestLine[0].point1_.y_;
	    TestLine[1].point1_.x_ = TestLine[0].point1_.x_;
	    TestLine[1].point1_.y_ = - TestLine[0].point1_.y_;
	    TestLine[1].point2_.x_ = - TestLine[1].point1_.x_;
	    TestLine[1].point2_.y_ = - TestLine[1].point1_.y_;

	    /* translate the lines back */
	    TempMat = IdentMat;
	    RotateMat(&TempMat,Phi);
	    TranslateMat(&TempMat,m,n);
	    TransformPoint(&TestLine[0].point1_,&TempMat);
	    TransformPoint(&TestLine[0].point2_,&TempMat);
	    TransformPoint(&TestLine[1].point1_,&TempMat);
	    TransformPoint(&TestLine[1].point2_,&TempMat);
	    }
	 }

      /* find the CEllipse2D line intersections */
      for (j = 0;j < NumLines;j++)
	 {
	 /* transform the line endpts into the circle space of the CEllipse2D */
	 TransformPoint(&TestLine[j].point1_,&ElpCirMat1);
	 TransformPoint(&TestLine[j].point2_,&ElpCirMat1);

	 /* compute the number of intersections of the transformed line
	    and test circle */
	 CircleInts = IntCirLine(&IntPts[IntCount],&TestCir,&TestLine[j]);
	 if (CircleInts>0)
	    {
	    /* transform the intersection points back into CEllipse2D space */
	    for (k = 0;k < CircleInts;k++)
	       TransformPoint(&IntPts[IntCount+k],&InvMat);
	    /* update the number of intersections found */
	    IntCount += CircleInts;
	    }
	 }
      }
   /* validate the points */
   j = IntCount;
   IntCount = 0;
   for (i = 0;i < j;i++)
      {
      TestPoint = IntPts[i];
      TransformPoint(&TestPoint,&ElpCirMat2);
      if (TestPoint.x_ < 2.0 && TestPoint.y_ < 2.0 &&
	 isZero(1.0 - CMath::sqrt(TestPoint.x_*TestPoint.x_ +
	 TestPoint.y_*TestPoint.y_)))
	 IntPts[IntCount++]=IntPts[i];
      }
   return IntCount;
   }

/* Test routines */

/* CEllipse2D with center at (1,2), major radius 2, minor radius 1,
   and rotation 0. */
   CPoint2D pc(2.0,1.0);
CEllipse2D TestMyEllipse(pc,2.0,1.0,0.0);

/* transform matrix for shear of 45 degrees from vertical */
TMat TestTransform(1.0f,0.0f,1.0f,1.0f,0.0f,0.0f);



/* test MyEllipses for intersection */
CEllipse2D Elp1(CPoint2D(5.0,4.0),1.0,0.5,CMath::PI/3.0);
CEllipse2D Elp2(CPoint2D(4.0,3.0),2.0,1.0,0.0);
CEllipse2D Elp3(CPoint2D(1.0,1.0),2.0,1.0,M_PI_2);
CEllipse2D Elp4(CPoint2D(1.0,1.0),2.0,0.5,0.0);
CEllipse2D Sal1(CPoint2D(0.0,1.0),4.0,2.0,0.0);
CEllipse2D Sal2(CPoint2D(0.0,2.0),CMath::sqrt(2),CMath::sqrt(12),0.0);

void CConic2D::testEllipseIntersection()
   {
  CPoint2D IntPts[12];
   int IntCount;
   int i;
   /* find CEllipse2D intersections */
   printf("Intersections of MyEllipses Sal1 & Sal2:\n\n");
   IntCount=intersects(Sal1,Sal2,IntPts);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].x_,IntPts[i].y_);

   /* find CEllipse2D intersections */
   printf("Intersections of MyEllipses 1 & 2:\n\n");
   IntCount=intersects(Elp1,Elp2,IntPts);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].x_,IntPts[i].y_);
   printf("Intersections of MyEllipses 3 & 4:\n");
   IntCount=intersects(Elp3,Elp4,IntPts);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].x_,IntPts[i].y_);

   /* transform MyEllipses */
   printf("\n\nBefore transformation:\n");
   printf(TestMyEllipse.toString().c_str());

   TransformMyEllipse(&TestMyEllipse,&TestTransform);

   printf("After transformation:\n");
   printf(TestMyEllipse.toString().c_str());

   }


} //end GEO
}  //end SMF
