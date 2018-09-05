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

#include "math/SMF_Polynomial.h"
//#include "util/SMF_Util.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_EulerAngles.h"
#include "util/SMF_StringUtils.h"
#ifndef FLT_EPSILON
#define FLT_EPSILON 1.0E-12f
#endif
#include <complex>
//#pragma hdrstop

namespace SMF {
namespace MATH{

const float EPSILON		= 1e-6f;

CPolyConSystem::CPolyConSystem():
sol1(CVec2D::nan),sol2(CVec2D::nan),sol3(CVec2D::nan),sol4(CVec2D::nan),
alpha(CMath::NAN_FLOAT), beta(CMath::NAN_FLOAT),gamma(CMath::NAN_FLOAT), delta(CMath::NAN_FLOAT), epsilon(CMath::NAN_FLOAT),
Y1(CMath::NAN_FLOAT),Y2(CMath::NAN_FLOAT), Y3(CMath::NAN_FLOAT), Y4(CMath::NAN_FLOAT),
X1A(CMath::NAN_FLOAT),X2A(CMath::NAN_FLOAT), X3A(CMath::NAN_FLOAT), X4A(CMath::NAN_FLOAT),
X1B(CMath::NAN_FLOAT),X2B(CMath::NAN_FLOAT), X3B(CMath::NAN_FLOAT), X4B(CMath::NAN_FLOAT),
numRealSol(0),solved(false){}

CPolyConSystem::CPolyConSystem(float a,float b, float c,float d, float e, float f,
							   float g,float h, float i,float j, float k, float l):
sol1(CVec2D::nan),sol2(CVec2D::nan),sol3(CVec2D::nan),sol4(CVec2D::nan),
A(a),
B(b),
C(c),
D(d),
E(e),
F(f),
G(g),
H(h),
I(i),
J(j),
K(k),
L(l),alpha(CMath::NAN_FLOAT), beta(CMath::NAN_FLOAT),gamma(CMath::NAN_FLOAT), delta(CMath::NAN_FLOAT), epsilon(CMath::NAN_FLOAT),
Y1(CMath::NAN_FLOAT),Y2(CMath::NAN_FLOAT), Y3(CMath::NAN_FLOAT), Y4(CMath::NAN_FLOAT),
X1A(CMath::NAN_FLOAT),X2A(CMath::NAN_FLOAT), X3A(CMath::NAN_FLOAT), X4A(CMath::NAN_FLOAT),
X1B(CMath::NAN_FLOAT),X2B(CMath::NAN_FLOAT), X3B(CMath::NAN_FLOAT), X4B(CMath::NAN_FLOAT),
numRealSol(0),solved(false){}

void CPolyConSystem::resetCoeff(){
sol1=CVec2D::nan;sol2=CVec2D::nan;sol3=CVec2D::nan;sol4=CVec2D::nan;
alpha=CMath::NAN_FLOAT; beta=CMath::NAN_FLOAT;gamma=CMath::NAN_FLOAT; delta=CMath::NAN_FLOAT; epsilon=CMath::NAN_FLOAT;
Y1=CMath::NAN_FLOAT;Y2=CMath::NAN_FLOAT; Y3=CMath::NAN_FLOAT; Y4=CMath::NAN_FLOAT;
X1A=CMath::NAN_FLOAT;X2A=CMath::NAN_FLOAT; X3A=CMath::NAN_FLOAT; X4A=CMath::NAN_FLOAT;
X1B=CMath::NAN_FLOAT;X2B=CMath::NAN_FLOAT; X3B=CMath::NAN_FLOAT; X4B=CMath::NAN_FLOAT;

numRealSol=0;
solved=false;
}
void CPolyConSystem::setEqu1(float a,float b, float c,float d, float e, float f){
A=a;
B=b;
C=c;
D=d;
E=e;
F=f;
resetCoeff();
}
void CPolyConSystem::setEqu2(float g,float h, float i,float j, float k, float l){
G=g;
H=h;
I=i;
J=j;
K=k;
L=l;
resetCoeff();
}

void CPolyConSystem::calcCoeff(){

	sf_double  alpha_d, beta_d, gamma_d,delta_d, epsilon_d;
	sf_double  Ad=static_cast<double>(A);sf_double  Bd=static_cast<double>(B);sf_double  Cd=static_cast<double>(C);sf_double  Dd=static_cast<double>(D);sf_double  Ed=static_cast<double>(E);sf_double  Fd=static_cast<double>(F);
	sf_double  Gd=static_cast<double>(G);sf_double  Hd=static_cast<double>(H);sf_double  Id=static_cast<double>(I);sf_double  Jd=static_cast<double>(J);sf_double  Kd=static_cast<double>(K);sf_double  Ld=static_cast<double>(L);
	if (false) { //A=0 \\todo
		
	}
	else {

		if (A == 0.0) {
			alpha_d = G * CMath::pow64(C, 2) - B * C * H + CMath::pow64(B, 2) * I;
			beta_d = 2.0 * C * E * G - C * D * H - B * C * J - B * E * H + CMath::pow64(B, 2) * K + 2.0 * B * D * I;
			gamma_d = G * CMath::pow64(E, 2) + 2.0 * C * F * G - C * D * J - D * E * H - B * F * H - B * E * J + I * CMath::pow64(D, 2) + CMath::pow64(B, 2) * L + 2.0 * B * D * K;
			delta_d = 2.0 * E * F * G - D * F * H - D * E * J - B * F * J + CMath::pow64(D, 2) * K + 2.0 * B * D * L;
			epsilon_d = G * CMath::pow64(F, 2) - D * F * J + CMath::pow64(D, 2) * L;
		}
		else {
			double A_2 = CMath::pow64(Ad, 2);
			double A_3 = CMath::pow64(Ad, 3);
			double Bd_2= CMath::pow64(Bd, 2);
			double D_2 = CMath::pow64(Dd, 2);
			double G_2 = CMath::pow64(Gd, 2);

			double AD_2 = CMath::pow64(Ad * Hd, 2);
			double BG_2 = CMath::pow64(Bd * Gd, 2);
			double AJ_2 = CMath::pow64(Ad * Jd, 2);
			double DG_2 = CMath::pow64(Dd * Gd, 2);
			double BJ_pl_DH = (Bd * Jd + Dd * Hd);
			double D2_4AF = (D_2 - 4.0 * Ad * Fd);
			double B2_4AC = (Bd_2 - 4.0 * Ad * Cd);
			double BD_2AE = (Bd * Dd - 2.0 * Ad * Ed);
			double BD_AE = (Bd * Dd - Ad * Ed);
			double A2K_K_ADH_ABJ = 2.0 * A_2 * Kd - Ad * Dd * Hd - Ad * Bd * Jd;
			double BD_AE_p_A2K_K_ADH_ABJ = (((BD_AE)* Gd) + A2K_K_ADH_ABJ);

			alpha_d = 4.0 * Gd * (Bd_2 - 2.0 * Ad * Cd) * (Gd * Bd_2 - 2.0 * Ad * Cd * Gd + 4.0 * A_2 * Id - 2.0 * Ad * Bd * Hd) + \
				16.0 * A_3 * (Ad * CMath::pow64(Id, 2) - Bd * Id * Hd) + \
				4.0 * CMath::pow64(Ad * Bd * Hd, 2) - 4.0 * B2_4AC * (AD_2 - 2.0 * Ad * Bd * Gd * Hd + BG_2);

			beta_d = 8.0 * Gd * ((BD_AE) * (Gd * Bd_2 - 2.0 * Ad * Cd * Gd + 4.0 * A_2 * Id - 2.0 * Ad * Bd * Hd) + \
				 BD_AE_p_A2K_K_ADH_ABJ * (Bd_2 - 2.0 * Ad * Cd)) + \
				16.0 * A_3 * (Ad * Id * Kd - Dd * Id * Hd - Bd * Id * Jd + Ad * Id * Kd - Bd * Hd * Kd) + 8.0 * A_2 * Bd * Hd * (Dd * Hd + Bd * Jd) - \
				8.0 * B2_4AC * (A_2 * Hd * Jd - Ad * Gd * BJ_pl_DH + Bd * Dd * G_2) - 8.0 * BD_2AE * (AD_2 - 2.0 * Ad * Bd * Gd * Hd + BG_2);

			gamma_d = 4.0 * Gd * ((Bd_2 - 2.0 * Ad * Cd) * (Gd * D_2 - 2.0 * Ad * Fd * Gd + 4.0 * A_2 * Ld - 2.0 * Ad * Dd * Jd) + \
				(Gd * Bd_2 - 2.0 * Ad * Cd * Gd + 4.0 * A_2 * Id - 2.0 * Ad * Bd * Hd) * (D_2 - 2.0 * Ad * Fd) + 4.0 * (BD_AE) * \
				 BD_AE_p_A2K_K_ADH_ABJ) + 16.0 * A_3 * (Ad * Id * Ld - Dd * Id * Jd + Ad * CMath::pow64(Kd, 2) - Dd * Hd * Kd - Bd * Jd * Kd + Ad * Id * Ld - Bd * Hd * Ld) + \
				4.0 * A_2 * (CMath::pow64((Dd * Hd + Bd * Jd), 2) + 2.0 * Bd * Dd * Hd * Jd) - 4.0 * B2_4AC * \
				(AJ_2 - 2.0 * Ad * Dd * Gd * Jd + DG_2) - 4.0 * D2_4AF * (AD_2 - 2.0 * Ad * Bd * Gd * Hd + \
				BG_2) - 16.0 * BD_2AE * (A_2 * Hd * Jd - Ad * Gd * BJ_pl_DH + Bd * Dd * G_2);

			delta_d = 8.0 * Gd * ((BD_AE) * (Gd * D_2 - 2.0 * Ad * Fd * Gd + 4.0 * A_2 * Ld - 2.0 * Ad * Dd * Jd) + \
				(D_2 - 2.0 * Ad * Fd) *  BD_AE_p_A2K_K_ADH_ABJ) + \
				16.0 * A_3 * (Ad * Kd * Ld - Dd * Jd * Kd + Ad * Kd * Ld - Dd * Hd * Ld - Bd * Jd * Ld) + 8.0 * A_2 * Dd * Jd * (Dd * Hd + Bd * Jd) - \
				8.0 * BD_2AE * (AJ_2 - 2.0 * Ad * Dd * Gd * Jd + DG_2) - 8.0 * D2_4AF * (A_2 * Hd * Jd - Ad * Gd * BJ_pl_DH + Bd * Dd * G_2);
			
			epsilon_d = 4.0 * Gd * (D_2 - 2.0 * Ad * Fd) * (Gd * D_2 - 2.0 * Ad * Fd * Gd + 4.0 * A_2 * Ld - 2.0 * Ad * Dd * Jd) + \
				16.0 * A_3 * (Ad * CMath::pow64(Ld, 2) - Dd * Jd * Ld) + 4.0 * CMath::pow64(Ad * Dd * Jd, 2) - 4.0 * D2_4AF * (AJ_2 - 2.0 * Ad * Dd * Gd * Jd + DG_2);

		}

alpha=static_cast<float>(alpha_d);
beta=static_cast<float>(beta_d);
gamma=static_cast<float>(gamma_d);
delta=static_cast<float>(delta_d);
epsilon=static_cast<float>(epsilon_d);
}
#if 0
	CVec5D vec(alpha,beta,gamma,delta,epsilon);
	vec.toNormal();
	alpha=vec.x;
	beta=vec.y;
	gamma=vec.z;
	delta=vec.s;
	epsilon=vec.t;
#endif

}

int CPolyConSystem::solveY(){
	calcCoeff();
	CComplex rt1=CComplex::nan;
	CComplex rt2=CComplex::nan;
	CComplex rt3=CComplex::nan;
	CComplex rt4=CComplex::nan;
	
	int sol=0;
	CPolynomial::solveQuarticArthur( alpha,beta,gamma,delta,epsilon,  rt1,rt2,rt3,rt4 );
	if(MATH::isFinite(rt1.real())){
		float val =rt1.real();
		if( !FLOAT_IS_DENORMAL(val)) Y1=val;
		sol++;
	}
	if(MATH::isFinite(rt2.real()) && !(rt2==rt1)){
		float val =rt2.real();
		if( !FLOAT_IS_DENORMAL(val)) Y2=val;
		sol++;
	}
	if(MATH::isFinite(rt3.real())&& !(rt3==rt1) && !(rt3==rt2)){
		float val =rt3.real();
		if( !FLOAT_IS_DENORMAL(val)) Y3=val;
		sol++;
	}
	if(MATH::isFinite(rt4.real())&& !(rt4==rt1) && !(rt4==rt2) && !(rt4==rt3)){
		float val =rt4.real();
		if( !FLOAT_IS_DENORMAL(val)) Y4=val;
		sol++;
	}
	return sol; //CPolynomial::solveQuartic(alpha,beta,gamma,delta,epsilon,Y1,Y2,Y3,Y4);

}

float CPolyConSystem::calcX_a0(float y) {
	float	x = -(C * CMath::pow(y, 2) + E * y + F) / (B * y + D);
	return x;
}
float CPolyConSystem::calcX1(float y) {
	float by= B * y;
	float by_d = by+D;
	float by_d2 = CMath::square(by_d);
	float delta_part1 = ( 4.0 * A * ((C *y*y) + (E * y) + F));
	float delta = by_d2 - delta_part1;
	float deltaSqr = CMath::sqrt(delta);
	float  x1 = (-(by_d) + deltaSqr ) / (2.0 * A);
				
	float	x = (-(by_d) + CMath::sqrt(by_d2 - (4.0 * A * ((C * CMath::pow(y, 2)) + (E * y) + F)))) / (2.0 * A);
	return x1;
}
float CPolyConSystem::calcX2(float y) {
	float by= B * y;
	float by_d = by+D;
	float by_d2 = CMath::square(by_d);
	float delta_part1 = ( 4.0 * A * ((C *y*y) + (E * y) + F));
	float delta = by_d2 - delta_part1;
	float deltaSqr = CMath::sqrt(delta);
	float  x1 = (-(by_d) - deltaSqr ) / (2.0 * A);

	float	x = (-(by_d) - CMath::sqrt(by_d2 - (4.0 * A * (C * CMath::pow(y, 2) + E * y + F)))) / (2.0 * A);
	return x1;
}

bool CPolyConSystem::testValues( float x1, float y){
	float val1 = CMath::fabs((A * CMath::square(x1)) + (B * x1 * y) + (C * CMath::square(y)) + (D * x1) +( E * y) + F);
	float val2 = CMath::fabs((G * CMath::square(x1)) + (H * x1 * y) + (I * CMath::square(y)) + (J * x1) +( K * y) + L);
	//bool t1 = CMath::nearZero(val1);
	//bool t2 = CMath::nearZero(val1);
	bool passed = (int)val1 <= 0 && (int)val2 <=0;
	return passed;
}

void CPolyConSystem::setSol(float x,float y){
	CVec2D sol(x,y);
	if(!sol1.isFinite()){
		sol1=sol;
		numRealSol++;
	}else if(!sol2.isFinite()){
		//verifica s é igual a solução 1
		if(sol==sol1) return;
		else {
			sol2=sol;
			numRealSol++;
		}

	}else if(!sol3.isFinite()){
		//verifica s é igual a uma solução anterior
		if(sol==sol1 || sol==sol2) return;
		else {
			sol3=sol;
			numRealSol++;
		}
	}else if(!sol4.isFinite()){
		//verifica se é igual a uma solução anterior
		if(sol==sol1 || sol==sol2 || sol==sol3) return;
		else {
			sol4=sol;
			numRealSol++;
		}
	}
}
void CPolyConSystem::solveSystem(int precision){
	int solucoes=solveY();
	float y1,y2,y3,y4;
	if (precision>0){
		y1=CMath::round(Y1,precision);
		y2=CMath::round(Y2,precision);
		y3=CMath::round(Y3,precision);
		y4=CMath::round(Y4,precision);
	}else{
		y1=Y1;
		y2=Y2;
		y3=Y3;
		y4=Y4;
	}
	float x1,x2;
	bool t1,t2;
	if(solucoes >0 && MATH::isFinite(y1)&& !FLOAT_IS_DENORMAL(y1)){ //calcula possiveis x
		if (A == 0.0) {
				if (false) {
					
				}
				else {
					x1 = x2 = calcX_a0(y1);
				}
			}else {
				x1 = calcX1(y1);
				x2 = calcX2(y1);
			}

			bool v1 = MATH::isValid(x1);
			bool v2 = MATH::isValid(x2);
			if (v1){t1=testValues(x1,y1);}else t1=false;
			if (v2){t2=testValues(x2,y1);}else t2=false;
			if (t1) {
				X1A = x1;
				setSol(x1,y1);
			}
			if (t2) {
				X1B = x2;
				setSol(x2,y1);
			}
	}
	if(solucoes >1 && MATH::isFinite(y2)&& !FLOAT_IS_DENORMAL(y2)){ //calcula possiveis x
		if (A == 0.0) {
				if (false) {
					
				}
				else {
					x1 = x2 = calcX_a0(y2);
				}
			}
			else {
				x1 = calcX1(y2);
				x2 = calcX2(y2);
			}
			bool v1 = MATH::isValid(x1);
			bool v2 = MATH::isValid(x2);
			if (v1){t1=testValues(x1,y2);}else t1=false;
			if (v2){t2=testValues(x2,y2);}else t2=false;
			if (t1) {
				X2A = x1;
				setSol(x1,y2);
			}
			if (t2) {
				X2B = x2;
				setSol(x2,y2);
			}
	}
	if(solucoes >2 && MATH::isFinite(y3) && !FLOAT_IS_DENORMAL(y3)){ //calcula possiveis x
		if (A == 0.0) {
				if (false) {
					
				}
				else {
					x1 = x2 = calcX_a0(y3);
				}
			}
			else {
				x1 = calcX1(y3);
				x2 = calcX2(y3);
			}
			bool v1 = MATH::isValid(x1);
			bool v2 = MATH::isValid(x2);
			if (v1){t1=testValues(x1,y3);}else t1=false;
			if (v2){t2=testValues(x2,y3);}else t2=false;
			if (t1) {
				X3A = x1;
				setSol(x1,y3);
			}
			if (t2) {
				X3B = x2;
				setSol(x2,y3);
			}
	}
	if(solucoes >3 &&  MATH::isFinite(y4)&& !FLOAT_IS_DENORMAL(y4)){ //calcula possiveis x
		if (A == 0.0) {
				if (false) {
					
				}
				else {
					x1 = x2 = calcX_a0(y4);
				}
			}
			else {
				x1 = calcX1(y4);
				x2 = calcX2(y4);
			}
			bool v1 = MATH::isValid(x1);
			bool v2 = MATH::isValid(x2);
			if (v1){t1=testValues(x1,y4);}else t1=false;
			if (v2){t2=testValues(x2,y4);}else t2=false;
			if (t1) {
				X4A = x1;
				setSol(x1,y4);
			}
			if (t2) {
				X4B = x2;
				setSol(x2,y4);
			}
	}


	solved=true;
}


int CPolyConSystem::getSol(CVec2D &so1, CVec2D &so2, CVec2D &so3, CVec2D &so4,int prec){
	if (!solved)solveSystem(prec);
	so1=sol1;
	so2=sol2;
	so3=sol3;
	so4=sol4;
	return numRealSol;
}

/*
=============
CPolynomial::Laguer
=============
*/
int CPolynomial::Laguer( const CComplex *coef, const int degree, CComplex &x ) {
	const int MT = 10, MAX_ITERATIONS = MT * 8;
	static const float frac[] = { 0.0f, 0.5f, 0.25f, 0.75f, 0.13f, 0.38f, 0.62f, 0.88f, 1.0f };
	int i, j;
	float abx, abp, abm, err;
	CComplex dx, cx, b, d, f, g, s, gps, gms, g2;

	for ( i = 1; i <= MAX_ITERATIONS; i++ ) {
		b = coef[degree];
		err = b.abs();
		d.toZero();
		f.toZero();
		abx = x.abs();
		for ( j = degree - 1; j >= 0; j-- ) {
			f = x * f + d;
			d = x * d + b;
			b = x * b + coef[j];
			err = b.abs() + abx * err;
		}
		if ( b.abs() < err * EPSILON ) {
			return i;
		}
		g = d / b;
		g2 = g * g;
		s = ( ( degree - 1 ) * ( degree * ( g2 - 2.0f * f / b ) - g2 ) ).sqrt();
		gps = g + s;
		gms = g - s;
		abp = gps.abs();
		abm = gms.abs();
		if ( abp < abm ) {
			gps = gms;
		}
		if ( MAX( abp, abm ) > 0.0f ) {
			dx = degree / gps;
		} else {
			dx = CMath::exp( CMath::log( 1.0f + abx ) ) * CComplex( CMath::cos( i ), CMath::sin( i ) );
		}
		cx = x - dx;
		if ( x == cx ) {
			return i;
		}
		if ( i % MT == 0 ) {
			x = cx;
		} else {
			x -= frac[i/MT] * dx;
		}
	}
	return i;
}

/*
=============
CPolynomial::getRoots
=============
*/
int CPolynomial::getRoots( CComplex *roots ) const {
	int i, j;
	CComplex x, b, c, *coef;

	coef = (CComplex *) _allocafloat16( ( degree + 1 ) * sizeof( CComplex ) );
	for ( i = 0; i <= degree; i++ ) {
		coef[i].set( coefficient[i], 0.0f );
	}

	for ( i = degree - 1; i >= 0; i-- ) {
		x.toZero();
		Laguer( coef, i + 1, x );
		if ( CMath::fabs( x.i ) < 2.0f * EPSILON * CMath::fabs( x.r ) ) {
			x.i = 0.0f;
		}
		roots[i] = x;
		b = coef[i+1];
		for ( j = i; j >= 0; j-- ) {
			c = coef[j];
			coef[j] = b;
			b = x * b + c;
		}
	}

	for ( i = 0; i <= degree; i++ ) {
		coef[i].set( coefficient[i], 0.0f );
	}
	for ( i = 0; i < degree; i++ ) {
		Laguer( coef, degree, roots[i] );
	}

	for ( i = 1; i < degree; i++ ) {
		x = roots[i];
		for ( j = i - 1; j >= 0; j-- ) {
			if ( roots[j].r <= x.r ) {
				break;
			}
			roots[j+1] = roots[j];
		}
		roots[j+1] = x;
	}

	return degree;
}

int CPolynomial::solveQuarticLaguer(float c1,float c2,float c3, float c4, float c5, CComplex *roots ) {
	int deg =4;
	int i, j;
	CComplex x, b, c, *coef;
	float coeffs[5];
	coeffs[4]=c1;
	coeffs[3]=c2;
	coeffs[2]=c3;
	coeffs[1]=c4;
	coeffs[0]=c5;
	coef = (CComplex *) _allocafloat16( ( deg + 1 ) * sizeof( CComplex ) );
	for ( i = 0; i <= deg; i++ ) {
		coef[i].set( coeffs[i], 0.0f );
	}

	for ( i = deg - 1; i >= 0; i-- ) {
		x.toZero();
		Laguer( coef, i + 1, x );
		if ( CMath::fabs( x.i ) < 2.0f * EPSILON * CMath::fabs( x.r ) ) {
			x.i = 0.0f;
		}
		roots[i] = x;
		b = coef[i+1];
		for ( j = i; j >= 0; j-- ) {
			c = coef[j];
			coef[j] = b;
			b = x * b + c;
		}
	}

	for ( i = 0; i <= deg; i++ ) {
		coef[i].set( coeffs[i], 0.0f );
	}
	for ( i = 0; i < deg; i++ ) {
		Laguer( coef, deg, roots[i] );
	}

	for ( i = 1; i < deg; i++ ) {
		x = roots[i];
		for ( j = i - 1; j >= 0; j-- ) {
			if ( roots[j].r <= x.r ) {
				break;
			}
			roots[j+1] = roots[j];
		}
		roots[j+1] = x;
	}

	return deg;
}

/*
=============
CPolynomial::getRoots
=============
*/
int CPolynomial::getRoots( float *roots ) const {
	int i, num;
	CComplex *complexRoots;

	switch( degree ) {
		case 0: return 0;
		case 1: return getRoots1( coefficient[1], coefficient[0], roots );
		case 2: return getRoots2( coefficient[2], coefficient[1], coefficient[0], roots );
		case 3: return getRoots3( coefficient[3], coefficient[2], coefficient[1], coefficient[0], roots );
		case 4: return getRoots4( coefficient[4], coefficient[3], coefficient[2], coefficient[1], coefficient[0], roots );
	}

	// The Abel-Ruffini theorem states that there is no general solution
	// in radicals to polynomial equations of degree five or higher.
	// A polynomial equation can be solved by radicals if and only if
	// its Galois group is a solvable group.

	complexRoots = (CComplex *) _allocafloat16( degree * sizeof( CComplex ) );

	getRoots( complexRoots );

	for ( num = i = 0; i < degree; i++ ) {
		if ( complexRoots[i].i == 0.0f ) {
			roots[i] = complexRoots[i].r;
			num++;
		}
	}
	return num;
}

/*
=============
CPolynomial::toString
=============
*/
const char *CPolynomial::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
=============
CPolynomial::Test
=============
*/
void CPolynomial::Test() {
	int i, num;
	float roots[4], value;
	CComplex complexRoots[4], complexValue;
	CPolynomial p;

	p = CPolynomial( -5.0f, 4.0f );
	num = p.getRoots( roots );
	for ( i = 0; i < num; i++ ) {
		value = p.getValue( roots[i] );
		SMF_ASSERT( CMath::fabs( value ) < 1e-4f );
	}

	p = CPolynomial( -5.0f, 4.0f, 3.0f );
	num = p.getRoots( roots );
	for ( i = 0; i < num; i++ ) {
		value = p.getValue( roots[i] );
		SMF_ASSERT( CMath::fabs( value ) < 1e-4f );
	}

	p = CPolynomial( 1.0f, 4.0f, 3.0f, -2.0f );
	num = p.getRoots( roots );
	for ( i = 0; i < num; i++ ) {
		value = p.getValue( roots[i] );
		SMF_ASSERT( CMath::fabs( value ) < 1e-4f );
	}

	p = CPolynomial( 5.0f, 4.0f, 3.0f, -2.0f );
	num = p.getRoots( roots );
	for ( i = 0; i < num; i++ ) {
		value = p.getValue( roots[i] );
		SMF_ASSERT( CMath::fabs( value ) < 1e-4f );
	}

	p = CPolynomial( -5.0f, 4.0f, 3.0f, 2.0f, 1.0f );
	num = p.getRoots( roots );
	for ( i = 0; i < num; i++ ) {
		value = p.getValue( roots[i] );
		SMF_ASSERT( CMath::fabs( value ) < 1e-4f );
	}

	p = CPolynomial( 1.0f, 4.0f, 3.0f, -2.0f );
	num = p.getRoots( complexRoots );
	for ( i = 0; i < num; i++ ) {
		complexValue = p.getValue( complexRoots[i] );
		SMF_ASSERT( CMath::fabs( complexValue.r ) < 1e-4f && CMath::fabs( complexValue.i ) < 1e-4f );
	}

	p = CPolynomial( 5.0f, 4.0f, 3.0f, -2.0f );
	num = p.getRoots( complexRoots );
	for ( i = 0; i < num; i++ ) {
		complexValue = p.getValue( complexRoots[i] );
		SMF_ASSERT( CMath::fabs( complexValue.r ) < 1e-4f && CMath::fabs( complexValue.i ) < 1e-4f );
	}
}


int CPolynomial::solveQuadratic(float a, float b, float c, float &root1, float &root2)
{
root1=root2=CMath::NAN_FLOAT;
	// ax^2 + bx + c == 0 => x = [ -b +/- sqrt(b^2 - 4ac) ] / 2a.

	///\todo numerical float issues: catastrophic cancellation can occur in the subtraction.
	float radicand = b*b - 4.f * a * c;
	if (radicand < -1e-6f) // add a small epsilon to allow the radicand to be slightly zero.
		return 0;
	float denom = 1.f / (2.f * a);
	if (radicand < 1e-6f) // Consider the radicand to be zero, and hence only one solution.
	{
		root1 = -b * denom;
		return 1;
	}
	radicand = CMath::sqrt(radicand);
	root1 = (-b + radicand) * denom;
	root2 = (-b - radicand) * denom;
	return 2;
}

int CPolynomial::solveQuadratic(float c[ 3 ], float s[ 2 ])
{
    float p, q, D;
	s[0]=s[1]=CMath::NAN_FLOAT;
    /* normal form: x^2 + px + q = 0 */

    p = c[ 1 ] / (2 * c[ 2 ]);
    q = c[ 0 ] / c[ 2 ];

    D = p * p - q;

    if (CMath::nearZero(D))
    {
	s[ 0 ] = - p;
	return 1;
    }
    else if (D < 0)
    {
	return 0;
    }
    else if (D > 0)
    {
	float sqrt_D = CMath::sqrt(D);

	s[ 0 ] =   sqrt_D - p;
	s[ 1 ] = - sqrt_D - p;
	return 2;
    }
}
int CPolynomial::solveCubic(float c3, float c2, float c1, float c0, float &root1, float &root2, float &root3)
{
	float entrada[4];
	float saida[3];
	entrada[0]= c0;
	entrada[1]=c1;
	entrada[2]=c2;
	entrada[3]=c3;
	saida[0]=saida[1]=saida[2]=CMath::NAN_FLOAT;
	int num = solveCubic(entrada,saida);
	root1=saida[0];
	root2=saida[1];
	root3=saida[2];
	return num;
}

int CPolynomial::solveCubic(float c[ 4 ], float s[ 3 ])
{
    int     i, num;
    float   sub;
    float  A, B, C;
    float  sq_A, p, q;
    float  cb_p, D;
	s[0]=s[1]=s[2]=CMath::NAN_FLOAT;
    /* divide by c3 to get normal form: x^3 + Ax^2 + Bx + C = 0 */

    A = c[ 2 ] / c[ 3 ];
    B = c[ 1 ] / c[ 3 ];
    C = c[ 0 ] / c[ 3 ];

    /*  substitute x = y - A/3 to eliminate quadric term:
	x^3 +px + q = 0 */

    sq_A = A * A;
    p = 1.0/3 * (- 1.0/3 * sq_A + B);
    q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);

    /* use Cardano's formula */

    cb_p = p * p * p;
    D = q * q + cb_p;

    if (CMath::nearZero(D))
    {
	if (CMath::nearZero(q)) /* one triple solution */
	{
	     s[ 0 ] = 0;
	    num = 1;
	}
	else /* one single and one double solution */
	{
	    double u = CMath::cbrt(-q);
	    s[ 0 ] = 2 * u;
	    s[ 1 ] = - u;
	     num = 2;
	}
    }
    else if (D < 0) /* Casus irreducibilis: three real solutions */
    {
	double phi = 1.0/3 * CMath::acos(-q / CMath::sqrt(-cb_p));
	double t = 2 * sqrt(-p);

	s[ 0 ] =   t * CMath::cos(phi);
	s[ 1 ] = - t * CMath::cos(phi + CMath::PI / 3);
	s[ 2 ] = - t * CMath::cos(phi - CMath::PI / 3);
	num = 3;
    }
    else /* one real solution */
    {
	double sqrt_D = sqrt(D);
	double u = CMath::cbrt(sqrt_D - q);
	double v = - CMath::cbrt(sqrt_D + q);

	s[ 0 ] = u + v;
	num = 1;
    }

    /* resubstitute */

    sub = 1.0/3 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;
    return num;
}

int CPolynomial::solveQuartic(float c4,float c3, float c2, float c1, float c0, float &root1, float &root2, float &root3, float &root4)
{
    float s[4];  //saida
	float  coeffs[ 4 ];
    float  z, u, v, sub;
    float  A, B, C, D;
    float  sq_A, p, q, r;
    int     i, num;
	s[0]=s[1]=s[2]=s[3]=CMath::NAN_FLOAT;
	root1=root2=root3=root4=CMath::NAN_FLOAT;
    /* divide by c4 to get the normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 */

    A = c3 / c4;
    B = c2 / c4;
    C = c1 / c4;
    D = c0 / c4;

    /*  substitute x = y - A/4 to eliminate cubic term:
	x^4 + px^2 + qx + r = 0 */

    sq_A = A * A;
    p = - 3.0/8 * sq_A + B;
    q = 1.0/8 * sq_A * A - 1.0/2 * A * B + C;
    r = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D;

    if (CMath::nearZero(r))
    {
	/* no absolute term: y(y^3 + py + q) = 0 */

	coeffs[ 0 ] = q;
	coeffs[ 1 ] = p;
	coeffs[ 2 ] = 0;
	coeffs[ 3 ] = 1;

	num = solveCubic(coeffs, s);

	s[ num++ ] = 0;
    }
    else
    {
	/* solve the resolvent cubic ... */

	coeffs[ 0 ] = 0.5f * r * p - 1.0/8 * q * q;
	coeffs[ 1 ] = - r;
	coeffs[ 2 ] = - 0.5f  * p;
	coeffs[ 3 ] = 1;

	(void) solveCubic(coeffs, s);

	/* ... and take the one real solution ... */

	z = s[ 0 ];

	/* ... to build two quadric equations */

	u = z * z - r;
	v = 2 * z - p;

	if (CMath::nearZero(u))
	    u = 0;
	else if (u > 0)
	    u = sqrt(u);
	else
	    return 0;

	if (CMath::nearZero(v))
	    v = 0;
	else if (v > 0)
	    v = sqrt(v);
	else
	    return 0;

	coeffs[ 0 ] = z - u;
	coeffs[ 1 ] = q < 0 ? -v : v;
	coeffs[ 2 ] = 1;

	num = solveQuadratic(coeffs, s);

	coeffs[ 0 ]= z + u;
	coeffs[ 1 ] = q < 0 ? v : -v;
	coeffs[ 2 ] = 1;

	num += solveQuadratic(coeffs, s + num);
    }

    /* resubstitute */

    sub = 1.0/4 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

	root1 = s[0];
	root2 = s[1];
	root3 = s[2];
	root4 = s[3];
    return num;
}

//==========ARTHUR================

using namespace std;


//método de Ferrari

void CPolynomial::solveSimple(const float& ALPHA, const float& BETA, CComplex &R1)
{
	R1 = CMath::NAN_FLOAT;

	if (ALPHA == 0.0f) {
		if (BETA == 0.0f) {
			cout << "Equacao possivel indeterminada";
		}
		else {
			cout << "Equacao Impossivel";
		}
	}
	else {
		R1 = -BETA / ALPHA;
	}
}

void CPolynomial::solveQuadraticArthur(const float& ALPHA, const float& BETA, const float& GAMMA, CComplex &R1, CComplex &R2)
{
	R1 = R2 = CMath::NAN_FLOAT;

	if (ALPHA == 0.0f) {
		solveSimple(BETA, GAMMA, R1);
	}
	else {
		CComplex D = sqrt(BETA * BETA - 4.0f * ALPHA * GAMMA);

		R1 = (-BETA + D) / (2.0f * ALPHA);
		R2 = (-BETA - D) / (2.0f * ALPHA);
	}
}

void CPolynomial::solveCubicArthur(const float& ALPHA, const float& BETA, const float& GAMMA, const float& DELTA, CComplex &RR1, CComplex &RR2, CComplex &RR3)
{
	RR1 = RR2 = RR3 = CMath::NAN_FLOAT;

	if (ALPHA == 0.0f) {
		solveQuadraticArthur(BETA, GAMMA, DELTA, RR1, RR2);
	}
	else {
		CComplex MI = CComplex::sqrt(static_cast<CComplex>((4.0f * BETA * BETA - 12.0f * ALPHA * GAMMA) / (9.0f * ALPHA * ALPHA)));

		if (MI == 0.0f) {
			RR1 = RR2 = RR3 = -BETA / (3.0f * ALPHA);
		}
		else {
			CComplex TETA = CComplex::arccos((36.0f * ALPHA * BETA * GAMMA - 8.0f * BETA * BETA * BETA - 108.0f * ALPHA * ALPHA * DELTA) / (27.0f * ALPHA * ALPHA * ALPHA * MI * MI * MI)) / 3.0f;
			//ccos aqui
			RR1 = MI * CComplex::cos(TETA) - BETA / (3.0f * ALPHA);
			TETA = TETA + 2.0f * 3.141592f / 3.0f;
			
			RR2 = MI * CComplex::cos(TETA) - BETA / (3.0f * ALPHA);
			TETA = TETA + 2.0f * 3.141592f / 3.0f;

			RR3 = MI * CComplex::cos(TETA) - BETA / (3.0f * ALPHA);
		}
	}
}

void CPolynomial::solveQuarticArthur(const float& ALPHA, const float& BETA, const float& GAMMA, const float& DELTA, const float& EPSILON, CComplex &RR1, CComplex &RR2, CComplex &RR3, CComplex &RR4) {
	RR1=CMath::NAN_FLOAT,RR2=CMath::NAN_FLOAT,RR3=CMath::NAN_FLOAT,RR4=CMath::NAN_FLOAT;
	//CComplex RRR1=CMath::NAN_FLOAT,RRR2=CMath::NAN_FLOAT,RRR3=CMath::NAN_FLOAT,RRR4=CMath::NAN_FLOAT;

	if (ALPHA == 0.0f) {
		solveCubicArthur(BETA, GAMMA, DELTA, EPSILON, RR1, RR2, RR3);
		
	}
	else {
		float P, Q, R;

		P = (8.0f * ALPHA * GAMMA - 3.0f * BETA * BETA) / (8.0f * ALPHA * ALPHA);
		Q = (BETA * BETA * BETA - 4.0f * ALPHA * BETA * GAMMA + 8.0f * ALPHA * ALPHA * DELTA) / (8.0f * ALPHA * ALPHA * ALPHA);
		R = (16.0f * ALPHA * BETA * BETA * GAMMA + 256.0f * ALPHA * ALPHA *ALPHA * EPSILON - 3.0f * BETA * BETA * BETA * BETA - 64.0f * ALPHA * ALPHA * BETA * DELTA) / (256.0f * ALPHA * ALPHA * ALPHA * ALPHA);

		
		CComplex MMI = 2.0f * CComplex::sqrt(static_cast<CComplex>(P * P + 12 * R)) / 3.0f;
		//CComplex MMMI = CComplex::sqrt(static_cast<CComplex>(4.0f * (P * P + 12 * R) / 9.0f));
		
		CComplex SSIGMA = MMI == 0.0f ? 0.0f : (3.0f * MMI * CComplex::cos(CComplex::arccos((8.0f * P * P * P - 288.0f * P * R + 108.0f * Q * Q) / (27.0f * MMI * MMI * MMI)) / 3.0f) - 2.0f * P) / 3.0f;
		//CComplex SSSIGMA = (MMMI == 0.0f ? 0.0f :MMMI * CComplex::cos(1.0f / 3.0f * CComplex::arccos(4 * (2.0f * P * P * P - 72.0f * P * R + 27.0f * Q * Q) / (27.0f * MMMI * MMMI * MMMI)))) - 2.0f * P / 3.0f;

		CComplex DD1= CComplex::sqrt(8.0f * (0.5f * SSIGMA * SSIGMA - SSIGMA * (P + SSIGMA) - CComplex::sqrt(SSIGMA) * Q));
		//CComplex DDD1 = CComplex::sqrt(4.0f * SSSIGMA * SSSIGMA - 8.0f * SSSIGMA * (P + SSSIGMA) - 8.0f * CComplex::sqrt(SSSIGMA) * Q);
		
		CComplex DD2 = CComplex::sqrt(8.0f * (0.5f * SSIGMA * SSIGMA - SSIGMA * (P + SSIGMA) + CComplex::sqrt(SSIGMA) * Q));
		//CComplex DDD2 = CComplex::sqrt(4.0f * SSSIGMA * SSSIGMA - 8.0f * SSSIGMA * (P + SSSIGMA) + 8.0f * CComplex::sqrt(SSSIGMA) * Q);


		RR1 = (4.0f * ALPHA * (2.0f * SSIGMA + DD1) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		RR2 = (4.0f * ALPHA * (2.0f * SSIGMA - DD1) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		RR3 = (4.0f * ALPHA * (-2.0f * SSIGMA + DD2) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		RR4 = (4.0f * ALPHA * (-2.0f * SSIGMA - DD2) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);

		//RRR1 = (4.0f * ALPHA * (2.0f * SSIGMA + DD1) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		//RRR2 = (4.0f * ALPHA * (2.0f * SSIGMA - DD1) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		//RRR3 = (4.0f * ALPHA * (-2.0f * SSIGMA + DD2) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);
		//RRR4 = (4.0f * ALPHA * (-2.0f * SSIGMA - DD2) / (4.0f * CComplex::sqrt(SSIGMA)) - BETA) / (4.0f * ALPHA);

	}

}
} // end MATH
} //end SMF
