
#include <io.h>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <values.h>
#include "SMF.h"

using namespace std;
using std::set_terminate;
using std::set_unexpected;
using std::set_new_handler;

/* CONMAT.C - MyEllipse tranformation and intersection functions */
/* Written by Kenneth J. Hill, June, 24, 1994  */


#ifndef	M_PI
#define M_PI	3.14159265358979323846
#endif

#ifndef	M_PI_2
#define M_PI_2	1.57079632679489661923
#endif

#ifndef max
#define max(a,b) ((a)>=(b) ? (a) : (b) )
#endif

#include    <math.h>
extern float   sqrt(), cbrt(), cos(), acos();

/* epsilon surrounding for near zero values */
#define NOCBRT
#define     EQN_EPS     1e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)

#ifdef NOCBRT
#define     cbrt(x)     ((x) > 0.0 ? pow((float)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -pow((float)-(x), 1.0/3.0) : 0.0))
#endif

int SolveQuadric(
    float c[ 3 ],
    float s[ 2 ])
{
    float p, q, D;

    /* normal form: x^2 + px + q = 0 */

    p = c[ 1 ] / (2 * c[ 2 ]);
    q = c[ 0 ] / c[ 2 ];

    D = p * p - q;

    if (IsZero(D))
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
	float sqrt_D = sqrt(D);

	s[ 0 ] =   sqrt_D - p;
	s[ 1 ] = - sqrt_D - p;
	return 2;
    }
}


int SolveCubic(
    float c[ 4 ],
    float s[ 3 ])
{
    int     i, num;
    float  sub;
    float  A, B, C;
    float  sq_A, p, q;
    float  cb_p, D;

    /* normal form: x^3 + Ax^2 + Bx + C = 0 */

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

    if (IsZero(D))
    {
	if (IsZero(q)) /* one triple solution */
	{
	    s[ 0 ] = 0;
	    num = 1;
	}
	else /* one single and one float solution */
	{
	    float u = cbrt(-q);
	    s[ 0 ] = 2 * u;
	    s[ 1 ] = - u;
	    num = 2;
	}
    }
    else if (D < 0) /* Casus irreducibilis: three real solutions */
    {
	float phi = 1.0/3 * acos(-q / sqrt(-cb_p));
	float t = 2 * sqrt(-p);

	s[ 0 ] =   t * cos(phi);
	s[ 1 ] = - t * cos(phi + M_PI / 3);
	s[ 2 ] = - t * cos(phi - M_PI / 3);
	num = 3;
    }
    else /* one real solution */
    {
	float sqrt_D = sqrt(D);
	float u = cbrt(sqrt_D - q);
	float v = - cbrt(sqrt_D + q);

	s[ 0 ] = u + v;
	num = 1;
    }

    /* resubstitute */

    sub = 1.0/3 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

    return num;
}


int SolveQuartic(
    float c[ 5 ],
    float s[ 4 ])
{
    float  coeffs[ 4 ];
    float  z, u, v, sub;
    float  A, B, C, D;
    float  sq_A, p, q, r;
    int     i, num;

    /* normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 */

    A = c[ 3 ] / c[ 4 ];
    B = c[ 2 ] / c[ 4 ];
    C = c[ 1 ] / c[ 4 ];
    D = c[ 0 ] / c[ 4 ];

    /*  substitute x = y - A/4 to eliminate cubic term:
	x^4 + px^2 + qx + r = 0 */

    sq_A = A * A;
    p = - 3.0/8 * sq_A + B;
    q = 1.0/8 * sq_A * A - 1.0/2 * A * B + C;
    r = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D;

    if (IsZero(r))
    {
	/* no absolute term: y(y^3 + py + q) = 0 */

	coeffs[ 0 ] = q;
	coeffs[ 1 ] = p;
	coeffs[ 2 ] = 0;
	coeffs[ 3 ] = 1;

	num = SolveCubic(coeffs, s);

	s[ num++ ] = 0;
    }
    else
    {
	/* solve the resolvent cubic ... */

	coeffs[ 0 ] = 1.0/2 * r * p - 1.0/8 * q * q;
	coeffs[ 1 ] = - r;
	coeffs[ 2 ] = - 1.0/2 * p;
	coeffs[ 3 ] = 1;

	(void) SolveCubic(coeffs, s);

	/* ... and take the one real solution ... */

	z = s[ 0 ];

	/* ... to build two quadric equations */

	u = z * z - r;
	v = 2 * z - p;

	if (IsZero(u))
	    u = 0;
	else if (u > 0)
	    u = sqrt(u);
	else
	    return 0;

	if (IsZero(v))
	    v = 0;
	else if (v > 0)
	    v = sqrt(v);
	else
	    return 0;

	coeffs[ 0 ] = z - u;
	coeffs[ 1 ] = q < 0 ? -v : v;
	coeffs[ 2 ] = 1;

	num = SolveQuadric(coeffs, s);

	coeffs[ 0 ]= z + u;
	coeffs[ 1 ] = q < 0 ? v : -v;
	coeffs[ 2 ] = 1;

	num += SolveQuadric(coeffs, s + num);
    }

    /* resubstitute */

    sub = 1.0/4 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

    return num;
}


/* Type definitions */

/* Point structure */
typedef struct PointTag
   {
   float X,Y;
   } Point;

/* Line structure */
typedef struct LineTag
   {
   Point P1,P2;
   } Line;

/* Circle structure */
typedef struct CircleTag
   {
   Point Center;
   float Radius;
   } Circle;

/* MyEllipse structure */
typedef struct EllipseTag
   {
   Point Center;		 /* ellipse center	 */
   float MaxRad,MinRad;	 /* major and minor axis */
   float Phi;			 /* major axis rotation  */
   } MyEllipse;

/* Conic coefficients structure */
typedef struct ConicTag
   {
   float A,B,C,D,E,F;
   } Conic;

/* transformation matrix type */
typedef struct TMatTag
   {
   float a,b,c,d;		 /* tranformation coefficients */
   float m,n;			 /* translation coefficients   */
   } TMat;

/* prototypes */

/* Functions with names beginning with a O:
 *
 * 1. Ignore ellipse center coordinates.
 *
 * 2. Ignore any translation components in tranformation matrices.
 *
 * 3. Assume that conic coeficients were generated by "O-name" functions.
 *
 * This keeps the computations relatively simple. The ellipse centers
 * can then be transformed separately as points.
 */

/* OTransformConic - Transform conic coefficients about the origin */
void OTransformConic(Conic *ConicP,TMat *TMatP);

/* OGenEllipseCoefs - Generate conic coefficients of an ellipse */
void OGenEllipseCoefs(MyEllipse *EllipseP,Conic *ConicP);

/* OGenEllipseGeom - Generates ellipse geometry from conic coefficients */
void OGenEllipseGeom(Conic *ConicP,MyEllipse *EllipseP);

/* TransformPoint - Transform a point using a tranformation matrix */
void TransformPoint(Point *PointP,TMat *TMatP);

/* TransformEllipse - Transform an ellipse using a tranformation matrix */
void TransformEllipse(MyEllipse *EllipseP,TMat *TMatP);



/* Identity matrix */
static const TMat IdentMat={1.0,0.0,0.0,1.0,0.0,0.0};

/* Transformation matrix routines */

/* Translate a matrix by m,n */
void TranslateMat(TMat *Mat,float m,float n)
   {
   Mat->m += m;
   Mat->n += n;
   }

/* Rotate a matrix by Phi */
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

/* Scale a matrix by sx, sy */
void ScaleMat(TMat *Mat,float sx,float sy)
   {
   Mat->a *= sx;
   Mat->b *= sy;
   Mat->c *= sx;
   Mat->d *= sy;
   Mat->m *= sx;
   Mat->n *= sy;
   }

/* TransformPoint - Transform a point using a tranformation matrix */
void TransformPoint(Point *PointP,TMat *TMatP)
   {
   Point TempPoint;

   TempPoint.X = PointP->X*TMatP->a + PointP->Y*TMatP->c + TMatP->m;
   TempPoint.Y = PointP->X*TMatP->b + PointP->Y*TMatP->d + TMatP->n;
   *PointP=TempPoint;
   }

/* Conic routines */

/* near zero test */
#define EPSILON 1e-9
#define IsZero(x) (x > -EPSILON && x < EPSILON)

/* GenEllipseCoefs - Generate conic coefficients of an ellipse */
static void GenEllipseCoefs(MyEllipse *elip,Conic *M)
   {
   float sqr_r1,sqr_r2;
   float sint,cost,sin2t,sqr_sint,sqr_cost;
   float cenx,ceny,sqr_cenx,sqr_ceny,invsqr_r1,invsqr_r2;

   /* common coeficients */
   sqr_r1 = elip->MaxRad;
   sqr_r2 = elip->MinRad;
   sqr_r1 *= sqr_r1;
   sqr_r2 *= sqr_r2;
   sint = sin(elip->Phi);
   cost = cos(elip->Phi);
   sin2t = 2.0*sint*cost;
   sqr_sint = sint*sint;
   sqr_cost = cost*cost;
   cenx = elip->Center.X;
   sqr_cenx = cenx*cenx;
   ceny = elip->Center.Y;
   sqr_ceny = ceny*ceny;
   invsqr_r1 = 1.0/sqr_r1;
   invsqr_r2 = 1.0/sqr_r2;

   /* Compute the coefficients. These formulae are the transformations
      on the unit circle written out long hand */
   M->A = sqr_cost/sqr_r1 + sqr_sint/sqr_r2;
   M->B = (sqr_r2-sqr_r1)*sin2t/(2.0*sqr_r1*sqr_r2);
   M->C = sqr_cost/sqr_r2 + sqr_sint/sqr_r1;
   M->D = -ceny*M->B-cenx*M->A;
   M->E = -cenx*M->B-ceny*M->C;
   M->F = -1.0 + (sqr_cenx + sqr_ceny)*(invsqr_r1 + invsqr_r2)/2.0 +
      (sqr_cost - sqr_sint)*(sqr_cenx - sqr_ceny)*(invsqr_r1 - invsqr_r2)/2.0 +
      cenx*ceny*(invsqr_r1 - invsqr_r2)*sin2t;
   }

/* Compute the transformation which turns an ellipse into a circle */
void Elp2Cir(MyEllipse *Elp,TMat *CirMat)
   {
   /* Start with identity matrix */
   *CirMat = IdentMat;
   /* Translate to origin */
   TranslateMat(CirMat,-Elp->Center.X,-Elp->Center.Y);
   /* Rotate into standard position */
   RotateMat(CirMat,-Elp->Phi);
   /* Scale into a circle. */
   ScaleMat(CirMat,1.0/Elp->MaxRad,1.0/Elp->MinRad);
   }

/* Compute the inverse of the transformation
   which turns an ellipse into a circle */
void InvElp2Cir(MyEllipse *Elp,TMat *InvMat)
   {
   /* Start with identity matrix */
   *InvMat = IdentMat;
   /* Scale back into an ellipse. */
   ScaleMat(InvMat,Elp->MaxRad,Elp->MinRad);
   /* Rotate */
   RotateMat(InvMat,Elp->Phi);
   /* Translate from origin */
   TranslateMat(InvMat,Elp->Center.X,Elp->Center.Y);
   }

/* OTransformConic - Transform conic coefficients about the origin	 */
/* This routine ignores the translation components of *TMatP and
   assumes the conic is "centered" at the origin (i.e. D,E=0, F=-1)	 */
/* The computations are just the matrix operations written out long hand */
/* This code assumes that the transformation is not degenerate		 */
void OTransformConic(Conic *ConicP,TMat *TMatP)
   {
   float A,B,C,Denom;

   /* common denominator for transformed cooefficients */
   Denom = TMatP->a*TMatP->d - TMatP->b*TMatP->c;
   Denom *= Denom;

   A = (ConicP->C*TMatP->b*TMatP->b - 2.0*ConicP->B*TMatP->b*TMatP->d +
      ConicP->A*TMatP->d*TMatP->d)/Denom;

   B = (-ConicP->C*TMatP->a*TMatP->b + ConicP->B*TMatP->b*TMatP->c +
      ConicP->B*TMatP->a*TMatP->d - ConicP->A*TMatP->c*TMatP->d)/Denom;

   C = (ConicP->C*TMatP->a*TMatP->a - 2.0*ConicP->B*TMatP->a*TMatP->c +
      ConicP->A*TMatP->c*TMatP->c)/Denom;

   ConicP->A=A;
   ConicP->B=B;
   ConicP->C=C;
   }

/* OGenEllipseCoefs - Generate conic coefficients of an ellipse */
/* The ellipse is assumed to be centered at the origin. */
void OGenEllipseCoefs(MyEllipse *EllipseP,Conic *ConicP)
   {
   float SinPhi = sin(EllipseP->Phi); /* sine of ellipse rotation   */
   float CosPhi = cos(EllipseP->Phi); /* cosine of ellipse rotation */
   float SqSinPhi = SinPhi*SinPhi;    /* square of sin(phi)	     */
   float SqCosPhi = CosPhi*CosPhi;    /* square of cos(phi)	     */
   float SqMaxRad = EllipseP->MaxRad*EllipseP->MaxRad;
   float SqMinRad = EllipseP->MinRad*EllipseP->MinRad;

   /* compute coefficients for the ellipse in standard position */
   ConicP->A = SqCosPhi/SqMaxRad + SqSinPhi/SqMinRad;
   ConicP->B = (1.0/SqMaxRad - 1.0/SqMinRad)*SinPhi*CosPhi;
   ConicP->C = SqCosPhi/SqMinRad + SqSinPhi/SqMaxRad;
   ConicP->D = ConicP->E = 0.0;
   ConicP->F = -1.0;
   }

/* OGenEllipseGeom - Generates ellipse geometry from conic coefficients */
/* This routine assumes the conic coefficients D=E=0, F=-1 */
void OGenEllipseGeom(Conic *ConicP,MyEllipse *EllipseP)
   {
   float Numer,Denom,Temp;
   TMat DiagTransform;		 /* transform diagonalization */
   Conic ConicT = *ConicP;	 /* temporary copy of conic coefficients */
   float SinPhi,CosPhi;

   /* compute new ellipse rotation */
   Numer = ConicT.B + ConicT.B;
   Denom = ConicT.A - ConicT.C;
   /* Phi = 1/2 atan(Numer/Denom) = 1/2 (pi/2 - atan(Denom/Numer)
      We use the form that keeps the argument to atan between -1 and 1 */

   EllipseP->Phi = 0.5*(fabs(Numer) < fabs(Denom)?
      atan(Numer/Denom):M_PI_2-atan(Denom/Numer));

   /* diagonalize the conic */
   SinPhi = sin(EllipseP->Phi);
   CosPhi = cos(EllipseP->Phi);
   DiagTransform.a = CosPhi;	 /* rotate by -Phi */
   DiagTransform.b = -SinPhi;
   DiagTransform.c = SinPhi;
   DiagTransform.d = CosPhi;
   DiagTransform.m = DiagTransform.n = 0.0;
   OTransformConic(&ConicT,&DiagTransform);

   /* compute new radii from diagonalized coefficients */
   EllipseP->MaxRad = 1.0/sqrt(ConicT.A);
   EllipseP->MinRad = 1.0/sqrt(ConicT.C);

   /* be sure MaxRad >= MinRad */
   if (EllipseP->MaxRad < EllipseP->MinRad)
      {
      Temp = EllipseP->MaxRad;	 /* exchange the radii */
      EllipseP->MaxRad = EllipseP->MinRad;
      EllipseP->MinRad = Temp;
      EllipseP->Phi += M_PI_2;	 /* adjust the rotation */
      }
   }

/* TransformEllipse - Transform an ellipse using a tranformation matrix */
void TransformEllipse(MyEllipse *EllipseP,TMat *TMatP)
   {
   Conic EllipseCoefs;

   /* generate the ellipse coefficients (using Center=origin) */
   OGenEllipseCoefs(EllipseP,&EllipseCoefs);

   /* transform the coefficients */
   OTransformConic(&EllipseCoefs,TMatP);

   /* turn the transformed coefficients back into geometry */
   OGenEllipseGeom(&EllipseCoefs,EllipseP);

   /* translate the center */
   TransformPoint(&EllipseP->Center,TMatP);
   }

/* MultMat3 - Multiply two 3x3 matrices */
void MultMat3(float *Mat1,float *Mat2, float *Result)
   {
   int i,j;

   for (i = 0;i < 3;i++)
      for (j = 0;j < 3;j++)
	 Result[i*3+j] = Mat1[i*3+0]*Mat2[0*3+j] +
	    Mat1[i*3+1]*Mat2[1*3+j] +
	    Mat1[i*3+2]*Mat2[2*3+j];
   }

/* Transform a conic by a tranformation matrix */
void TransformConic(Conic *ConicP,TMat *TMatP)
   {
   float InvMat[3][3],ConMat[3][3],TranInvMat[3][3];
   float Result1[3][3],Result2[3][3];
   float D;
   int i,j;

   /* Compute M' = Inv(TMat).M.Transpose(Inv(TMat))

   /* compute the transformation using matrix muliplication */
   ConMat[0][0] = ConicP->A;
   ConMat[0][1] = ConMat[1][0] = ConicP->B;
   ConMat[1][1] = ConicP->C;
   ConMat[0][2] = ConMat[2][0] = ConicP->D;
   ConMat[1][2] = ConMat[2][1] = ConicP->E;
   ConMat[2][2] = ConicP->F;

   /* inverse transformation */
   D = TMatP->a*TMatP->d - TMatP->b*TMatP->c;
   InvMat[0][0] = TMatP->d/D;
   InvMat[0][1] = -TMatP->b/D;
   InvMat[0][2] = 0.0;
   InvMat[1][0] = -TMatP->c/D;
   InvMat[1][1] = TMatP->a/D;
   InvMat[1][2] = 0.0;
   InvMat[2][0] = (TMatP->c*TMatP->n - TMatP->d*TMatP->m)/D;
   InvMat[2][1] = (TMatP->b*TMatP->m - TMatP->a*TMatP->n)/D;
   InvMat[2][2] = 1.0;

   /* compute transpose */
   for (i = 0;i < 3;i++)
      for (j = 0;j < 3;j++)
	 TranInvMat[j][i] = InvMat[i][j];

   /* multiply the matrices */
   MultMat3((float *)InvMat,(float *)ConMat,(float *)Result1);
   MultMat3((float *)Result1,(float *)TranInvMat,(float *)Result2);
   ConicP->A = Result2[0][0];	       /* return to conic form */
   ConicP->B = Result2[0][1];
   ConicP->C = Result2[1][1];
   ConicP->D = Result2[0][2];
   ConicP->E = Result2[1][2];
   ConicP->F = Result2[2][2];
   }

/* Compute the intersection of a circle and a line */
/* See Graphic Gems Volume 1, page 5 for a description of this algorithm */
int IntCirLine(Point *IntPts,Circle *Cir,Line *Ln)
   {
   Point G,V;
   float a,b,c,d,t,sqrt_d;

   G.X = Ln->P1.X - Cir->Center.X;     	/* G = Ln->P1 - Cir->Center */
   G.Y = Ln->P1.Y - Cir->Center.Y;
   V.X = Ln->P2.X - Ln->P1.X;           /* V = Ln->P2 - Ln->P1 */
   V.Y = Ln->P2.Y - Ln->P1.Y;
   a = V.X*V.X + V.Y*V.Y;		/* a = V.V */
   b = V.X*G.X + V.Y*G.Y;b += b; 	/* b = 2(V.G) */
   c = (G.X*G.X + G.Y*G.Y) -    	/* c = G.G + Circle->Radius^2 */
      Cir->Radius*Cir->Radius;
   d = b*b - 4.0*a*c;			/* discriminant */

   if (d <= 0.0)
      return 0;				/* no intersections */

   sqrt_d = sqrt(d);
   t = (-b + sqrt_d)/(a + a);           /* t = (-b +/- sqrt(d))/2a */
   IntPts[0].X = Ln->P1.X + t*V.X;      /* Pt = Ln->P1 + t V */
   IntPts[0].Y = Ln->P1.Y + t*V.Y;
   t = (-b - sqrt_d)/(a + a);
   IntPts[1].X = Ln->P1.X + t*V.X;
   IntPts[1].Y = Ln->P1.Y + t*V.Y;
   return 2;
   }

/* compute all intersections of two ellipses */
/* E1 and E2 are the two ellipses */
/* IntPts points to an array of twelve points
    (some duplicates may be returned) */
/* The number of intersections found is returned */
/* Both ellipses are assumed to have non-zero radii */
int Int2Elip(Point *IntPts,MyEllipse *E1,MyEllipse *E2)
   {
   TMat ElpCirMat1,ElpCirMat2,InvMat,TempMat;
   Conic Conic1,Conic2,Conic3,TempConic;
   float Roots[3],qRoots[2];
   static Circle TestCir = {{0.0,0.0},1.0};
   Line TestLine[2];
   Point TestPoint;
   float PolyCoef[4];		/* coefficients of the polynomial */
   float D;			/* discriminant: B^2 - AC */
   float Phi;			/* ellipse rotation */
   float m,n;			/* ellipse translation */
   float Scl;			/* scaling factor */
   int NumRoots,NumLines;
   int CircleInts;		/* intersections between line and circle */
   int IntCount = 0;		/* number of intersections found */
   int i,j,k;

   /* compute the transformations which turn E1 and E2 into circles */
   Elp2Cir(E1,&ElpCirMat1);
   Elp2Cir(E2,&ElpCirMat2);

   /* compute the inverse transformation of ElpCirMat1 */
   InvElp2Cir(E1,&InvMat);

   /* Compute the characteristic matrices */
   GenEllipseCoefs(E1,&Conic1);
   GenEllipseCoefs(E2,&Conic2);

   /* Find x such that Det(Conic1 + x Conic2) = 0 */
   PolyCoef[0] = -Conic1.C*Conic1.D*Conic1.D + 2.0*Conic1.B*Conic1.D*Conic1.E -
      Conic1.A*Conic1.E*Conic1.E - Conic1.B*Conic1.B*Conic1.F +
      Conic1.A*Conic1.C*Conic1.F;
   PolyCoef[1] = -(Conic2.C*Conic1.D*Conic1.D) -
      2.0*Conic1.C*Conic1.D*Conic2.D + 2.0*Conic2.B*Conic1.D*Conic1.E +
      2.0*Conic1.B*Conic2.D*Conic1.E - Conic2.A*Conic1.E*Conic1.E +
      2.0*Conic1.B*Conic1.D*Conic2.E - 2.0*Conic1.A*Conic1.E*Conic2.E -
      2.0*Conic1.B*Conic2.B*Conic1.F + Conic2.A*Conic1.C*Conic1.F +
      Conic1.A*Conic2.C*Conic1.F - Conic1.B*Conic1.B*Conic2.F +
      Conic1.A*Conic1.C*Conic2.F;
   PolyCoef[2] = -2.0*Conic2.C*Conic1.D*Conic2.D - Conic1.C*Conic2.D*Conic2.D +
      2.0*Conic2.B*Conic2.D*Conic1.E + 2.0*Conic2.B*Conic1.D*Conic2.E +
      2.0*Conic1.B*Conic2.D*Conic2.E - 2.0*Conic2.A*Conic1.E*Conic2.E -
      Conic1.A*Conic2.E*Conic2.E - Conic2.B*Conic2.B*Conic1.F +
      Conic2.A*Conic2.C*Conic1.F - 2.0*Conic1.B*Conic2.B*Conic2.F +
      Conic2.A*Conic1.C*Conic2.F + Conic1.A*Conic2.C*Conic2.F;
   PolyCoef[3] = -Conic2.C*Conic2.D*Conic2.D + 2.0*Conic2.B*Conic2.D*Conic2.E -
      Conic2.A*Conic2.E*Conic2.E - Conic2.B*Conic2.B*Conic2.F +
      Conic2.A*Conic2.C*Conic2.F;
   NumRoots = SolveCubic(PolyCoef,Roots);

   if (NumRoots == 0)
      return 0;

   /* we try all the roots, even though it's redundant, so that we
      avoid some pathological situations */
   for (i=0;i<NumRoots;i++)
      {
      NumLines = 0;

      /* Conic3 = Conic1 + mu Conic2 */
      Conic3.A = Conic1.A + Roots[i]*Conic2.A;
      Conic3.B = Conic1.B + Roots[i]*Conic2.B;
      Conic3.C = Conic1.C + Roots[i]*Conic2.C;
      Conic3.D = Conic1.D + Roots[i]*Conic2.D;
      Conic3.E = Conic1.E + Roots[i]*Conic2.E;
      Conic3.F = Conic1.F + Roots[i]*Conic2.F;

      D = Conic3.B*Conic3.B - Conic3.A*Conic3.C;
      if (IsZero(Conic3.A) && IsZero(Conic3.B) && IsZero(Conic3.C))
	 {
	 /* Case 1 - Single line */
	 NumLines = 1;
	 /* compute endpoints of the line, avoiding division by zero */
	 if (fabs(Conic3.D) > fabs(Conic3.E))
	    {
	    TestLine[0].P1.Y = 0.0;
	    TestLine[0].P1.X = -Conic3.F/(Conic3.D + Conic3.D);
	    TestLine[0].P2.Y = 1.0;
	    TestLine[0].P2.X = -(Conic3.E + Conic3.E + Conic3.F)/
	       (Conic3.D + Conic3.D);
	    }
	 else
	    {
	    TestLine[0].P1.X = 0.0;
	    TestLine[0].P1.Y = -Conic3.F/(Conic3.E + Conic3.E);
	    TestLine[0].P2.X = 1.0;
	    TestLine[0].P2.X = -(Conic3.D + Conic3.D + Conic3.F)/
	       (Conic3.E + Conic3.E);
	    }
	 }
      else
	 {
	 /* use the espresion for Phi that takes atan of the
	    smallest argument */
	 Phi = (fabs(Conic3.B + Conic3.B) < fabs(Conic3.A-Conic3.C)?
	    atan((Conic3.B + Conic3.B)/(Conic3.A - Conic3.C)):
	    M_PI_2 - atan((Conic3.A - Conic3.C)/(Conic3.B + Conic3.B)))/2.0;
	 if (IsZero(D))
	    {
	    /* Case 2 - Parallel lines */
	    TempConic = Conic3;
	    TempMat = IdentMat;
	    RotateMat(&TempMat,-Phi);
	    TransformConic(&TempConic,&TempMat);
	    if (IsZero(TempConic.C))   /* vertical */
	       {
	       PolyCoef[0] = TempConic.F;
	       PolyCoef[1] = TempConic.D;
	       PolyCoef[2] = TempConic.A;
	       if ((NumLines=SolveQuadric(PolyCoef,qRoots))!=0)
		  {
		  TestLine[0].P1.X = qRoots[0];
		  TestLine[0].P1.Y = -1.0;
		  TestLine[0].P2.X = qRoots[0];
		  TestLine[0].P2.Y = 1.0;
		  if (NumLines==2)
		     {
		     TestLine[1].P1.X = qRoots[1];
		     TestLine[1].P1.Y = -1.0;
		     TestLine[1].P2.X = qRoots[1];
		     TestLine[1].P2.Y = 1.0;
		     }
		  }
	       }
	    else		    /* horizontal */
	       {
	       PolyCoef[0] = TempConic.F;
	       PolyCoef[1] = TempConic.E;
	       PolyCoef[2] = TempConic.C;
	       if ((NumLines=SolveQuadric(PolyCoef,qRoots))!=0)
		  {
		  TestLine[0].P1.X = -1.0;
		  TestLine[0].P1.Y = qRoots[0];
		  TestLine[0].P2.X = 1.0;
		  TestLine[0].P2.Y = qRoots[0];
		  if (NumLines==2)
		     {
		     TestLine[1].P1.X = -1.0;
		     TestLine[1].P1.Y = qRoots[1];
		     TestLine[1].P2.X = 1.0;
		     TestLine[1].P2.Y = qRoots[1];
		     }
		  }
	       }
	    TempMat = IdentMat;
	    RotateMat(&TempMat,Phi);
	    TransformPoint(&TestLine[0].P1,&TempMat);
	    TransformPoint(&TestLine[0].P2,&TempMat);
	    if (NumLines==2)
	       {
	       TransformPoint(&TestLine[1].P1,&TempMat);
	       TransformPoint(&TestLine[1].P2,&TempMat);
	       }
	    }
	 else
	    {
	    /* Case 3 - Crossing lines */
	    NumLines = 2;

	    /* translate the system so that the intersection of the lines
	       is at the origin */
	    TempConic = Conic3;
	    m = (Conic3.C*Conic3.D - Conic3.B*Conic3.E)/D;
	    n = (Conic3.A*Conic3.E - Conic3.B*Conic3.D)/D;
	    TempMat = IdentMat;
	    TranslateMat(&TempMat,-m,-n);
	    RotateMat(&TempMat,-Phi);
	    TransformConic(&TempConic,&TempMat);

	    /* Compute the line endpoints */
	    TestLine[0].P1.X = sqrt(fabs(1.0/TempConic.A));
	    TestLine[0].P1.Y = sqrt(fabs(1.0/TempConic.C));
	    Scl = max(TestLine[0].P1.X,TestLine[0].P1.Y);  /* adjust range */
	    TestLine[0].P1.X /= Scl;
	    TestLine[0].P1.Y /= Scl;
	    TestLine[0].P2.X = - TestLine[0].P1.X;
	    TestLine[0].P2.Y = - TestLine[0].P1.Y;
	    TestLine[1].P1.X = TestLine[0].P1.X;
	    TestLine[1].P1.Y = - TestLine[0].P1.Y;
	    TestLine[1].P2.X = - TestLine[1].P1.X;
	    TestLine[1].P2.Y = - TestLine[1].P1.Y;

	    /* translate the lines back */
	    TempMat = IdentMat;
	    RotateMat(&TempMat,Phi);
	    TranslateMat(&TempMat,m,n);
	    TransformPoint(&TestLine[0].P1,&TempMat);
	    TransformPoint(&TestLine[0].P2,&TempMat);
	    TransformPoint(&TestLine[1].P1,&TempMat);
	    TransformPoint(&TestLine[1].P2,&TempMat);
	    }
	 }

      /* find the ellipse line intersections */
      for (j = 0;j < NumLines;j++)
	 {
	 /* transform the line endpts into the circle space of the ellipse */
	 TransformPoint(&TestLine[j].P1,&ElpCirMat1);
	 TransformPoint(&TestLine[j].P2,&ElpCirMat1);

	 /* compute the number of intersections of the transformed line
	    and test circle */
	 CircleInts = IntCirLine(&IntPts[IntCount],&TestCir,&TestLine[j]);
	 if (CircleInts>0)
	    {
	    /* transform the intersection points back into ellipse space */
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
      if (TestPoint.X < 2.0 && TestPoint.Y < 2.0 &&
	 IsZero(1.0 - sqrt(TestPoint.X*TestPoint.X +
	 TestPoint.Y*TestPoint.Y)))
	 IntPts[IntCount++]=IntPts[i];
      }
   return IntCount;
   }

/* Test routines */

/* MyEllipse with center at (1,2), major radius 2, minor radius 1,
   and rotation 0. */
MyEllipse TestEllipse={{2.0,1.0},2.0,1.0,0.0};

/* Transform matrix for shear of 45 degrees from vertical */
TMat TestTransform={1.0,0.0,1.0,1.0,0.0,0.0};

/* Display an ellipse. This version lists the structure values. */
void DisplayEllipse(MyEllipse *EllipseP)
   {
   printf("\tCenter at (%6.3f,%6.3f)\n"
      "\tMajor radius: %6.3f\n"
      "\tMinor radius: %6.3f\n"
      "\tRotation: %6.3f degrees\n\n",
      EllipseP->Center.X,EllipseP->Center.Y,
      EllipseP->MaxRad,EllipseP->MinRad,
      180.0*EllipseP->Phi/M_PI);
   }


/* test ellipses for intersection */
MyEllipse Elp1={{5.0,4.0},1.0,0.5,M_PI/3.0};
MyEllipse Elp2={{4.0,3.0},2.0,1.0,0.0};
MyEllipse Elp3={{1.0,1.0},2.0,1.0,M_PI_2};
MyEllipse Elp4={{1.0,1.0},2.0,0.5,0.0};
MyEllipse Sal1={{0.0,-1.0},8.0,2.0,0.0};
MyEllipse Sal2={{0.0,-2.0},sqrt(12),sqrt(2),M_PI_2};

//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{

		//==============================
	SMF::Debug::setDebugAll();
	SMF::Debug::setFilename("debug_math.txt");
	SMF::Debug::setDebugMethod(SMF::Debug::File);

	
	SMF::MATH::CSIMD::initHeap();
   Point IntPts[12];
   int IntCount;
   int i;
   for (i = 0;i < 11;i++) { IntPts[i].X=0; IntPts[i].Y=0;}
   /* find ellipse intersections */
    printf("Intersections of ellipses Sal1 & Sal2:\n\n");
   IntCount=Int2Elip(IntPts,&Sal1,&Sal2);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].X,IntPts[i].Y);
  printf("Intersections of ellipses 1 & 2:\n\n");
   IntCount=Int2Elip(IntPts,&Elp1,&Elp2);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].X,IntPts[i].Y);
   printf("Intersections of ellipses 3 & 4:\n");
   IntCount=Int2Elip(IntPts,&Elp3,&Elp4);
   for (i = 0;i < IntCount;i++)
      printf("   %f, %f\n",IntPts[i].X,IntPts[i].Y);

   /* transform ellipses */
   printf("\n\nBefore transformation:\n");
   DisplayEllipse(&TestEllipse);

   TransformEllipse(&TestEllipse,&TestTransform);

   printf("After transformation:\n");
   DisplayEllipse(&TestEllipse);

	SMF::MATH::CSIMD::shutdown();

    return 0;

}



