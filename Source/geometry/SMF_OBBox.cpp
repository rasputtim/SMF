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


#include "geometry/SMF_OBBox.h"
#include "math/SMF_MathDefs.h"
#include "math/SMF_JointTransform.h"
#include "geometry/all.h"

namespace SMF{
using namespace MATH;
namespace GEO {

COBBox box_zero( CVec3D::zero, CVec3D::zero, mat3_identity );

/*
            4---{4}---5
 +         /|        /|
 Z      {7} {8}   {5} |
 -     /    |    /    {9}
      7--{6}----6     |
      |     |   |     |
    {11}    0---|-{0}-1
      |    /    |    /       -
      | {3}  {10} {1}       Y
      |/        |/         +
      3---{2}---2

	    - X +

  plane bits:
  0 = min x
  1 = max x
  2 = min y
  3 = max y
  4 = min z
  5 = max z

*/

/*
static int boxVertPlanes[8] = {
	( (1<<0) | (1<<2) | (1<<4) ),
	( (1<<1) | (1<<2) | (1<<4) ),
	( (1<<1) | (1<<3) | (1<<4) ),
	( (1<<0) | (1<<3) | (1<<4) ),
	( (1<<0) | (1<<2) | (1<<5) ),
	( (1<<1) | (1<<2) | (1<<5) ),
	( (1<<1) | (1<<3) | (1<<5) ),
	( (1<<0) | (1<<3) | (1<<5) )
};

static int boxVertEdges[8][3] = {
	// bottom
	{ 3, 0, 8 },
	{ 0, 1, 9 },
	{ 1, 2, 10 },
	{ 2, 3, 11 },
	// top
	{ 7, 4, 8 },
	{ 4, 5, 9 },
	{ 5, 6, 10 },
	{ 6, 7, 11 }
};

static int boxEdgePlanes[12][2] = {
	// bottom
	{ 4, 2 },
	{ 4, 1 },
	{ 4, 3 },
	{ 4, 0 },
	// top
	{ 5, 2 },
	{ 5, 1 },
	{ 5, 3 },
	{ 5, 0 },
	// sides
	{ 0, 2 },
	{ 2, 1 },
	{ 1, 3 },
	{ 3, 0 }
};

static int boxEdgeVerts[12][2] = {
	// bottom
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 0 },
	// top
	{ 4, 5 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 4 },
	// sides
	{ 0, 4 },
	{ 1, 5 },
	{ 2, 6 },
	{ 3, 7 }
};
*/

static int boxPlaneBitsSilVerts[64][7] = {
	{ 0, 0, 0, 0, 0, 0, 0 }, // 000000 = 0
	{ 4, 7, 4, 0, 3, 0, 0 }, // 000001 = 1
	{ 4, 5, 6, 2, 1, 0, 0 }, // 000010 = 2
	{ 0, 0, 0, 0, 0, 0, 0 }, // 000011 = 3
	{ 4, 4, 5, 1, 0, 0, 0 }, // 000100 = 4
	{ 6, 3, 7, 4, 5, 1, 0 }, // 000101 = 5
	{ 6, 4, 5, 6, 2, 1, 0 }, // 000110 = 6
	{ 0, 0, 0, 0, 0, 0, 0 }, // 000111 = 7
	{ 4, 6, 7, 3, 2, 0, 0 }, // 001000 = 8
	{ 6, 6, 7, 4, 0, 3, 2 }, // 001001 = 9
	{ 6, 5, 6, 7, 3, 2, 1 }, // 001010 = 10
	{ 0, 0, 0, 0, 0, 0, 0 }, // 001011 = 11
	{ 0, 0, 0, 0, 0, 0, 0 }, // 001100 = 12
	{ 0, 0, 0, 0, 0, 0, 0 }, // 001101 = 13
	{ 0, 0, 0, 0, 0, 0, 0 }, // 001110 = 14
	{ 0, 0, 0, 0, 0, 0, 0 }, // 001111 = 15
	{ 4, 0, 1, 2, 3, 0, 0 }, // 010000 = 16
	{ 6, 0, 1, 2, 3, 7, 4 }, // 010001 = 17
	{ 6, 3, 2, 6, 5, 1, 0 }, // 010010 = 18
	{ 0, 0, 0, 0, 0, 0, 0 }, // 010011 = 19
	{ 6, 1, 2, 3, 0, 4, 5 }, // 010100 = 20
	{ 6, 1, 2, 3, 7, 4, 5 }, // 010101 = 21
	{ 6, 2, 3, 0, 4, 5, 6 }, // 010110 = 22
	{ 0, 0, 0, 0, 0, 0, 0 }, // 010111 = 23
	{ 6, 0, 1, 2, 6, 7, 3 }, // 011000 = 24
	{ 6, 0, 1, 2, 6, 7, 4 }, // 011001 = 25
	{ 6, 0, 1, 5, 6, 7, 3 }, // 011010 = 26
	{ 0, 0, 0, 0, 0, 0, 0 }, // 011011 = 27
	{ 0, 0, 0, 0, 0, 0, 0 }, // 011100 = 28
	{ 0, 0, 0, 0, 0, 0, 0 }, // 011101 = 29
	{ 0, 0, 0, 0, 0, 0, 0 }, // 011110 = 30
	{ 0, 0, 0, 0, 0, 0, 0 }, // 011111 = 31
	{ 4, 7, 6, 5, 4, 0, 0 }, // 100000 = 32
	{ 6, 7, 6, 5, 4, 0, 3 }, // 100001 = 33
	{ 6, 5, 4, 7, 6, 2, 1 }, // 100010 = 34
	{ 0, 0, 0, 0, 0, 0, 0 }, // 100011 = 35
	{ 6, 4, 7, 6, 5, 1, 0 }, // 100100 = 36
	{ 6, 3, 7, 6, 5, 1, 0 }, // 100101 = 37
	{ 6, 4, 7, 6, 2, 1, 0 }, // 100110 = 38
	{ 0, 0, 0, 0, 0, 0, 0 }, // 100111 = 39
	{ 6, 6, 5, 4, 7, 3, 2 }, // 101000 = 40
	{ 6, 6, 5, 4, 0, 3, 2 }, // 101001 = 41
	{ 6, 5, 4, 7, 3, 2, 1 }, // 101010 = 42
	{ 0, 0, 0, 0, 0, 0, 0 }, // 101011 = 43
	{ 0, 0, 0, 0, 0, 0, 0 }, // 101100 = 44
	{ 0, 0, 0, 0, 0, 0, 0 }, // 101101 = 45
	{ 0, 0, 0, 0, 0, 0, 0 }, // 101110 = 46
	{ 0, 0, 0, 0, 0, 0, 0 }, // 101111 = 47
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110000 = 48
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110001 = 49
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110010 = 50
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110011 = 51
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110100 = 52
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110101 = 53
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110110 = 54
	{ 0, 0, 0, 0, 0, 0, 0 }, // 110111 = 55
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111000 = 56
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111001 = 57
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111010 = 58
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111011 = 59
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111100 = 60
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111101 = 61
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111110 = 62
	{ 0, 0, 0, 0, 0, 0, 0 }, // 111111 = 63
};


/*
============
COBBox::addPoint
============
*/
bool COBBox::addPoint( const CVec3D &v ) {
	CMat3D axis2;
	CAABBox bounds1, bounds2;

	if ( extents[0] < 0.0f ) {
		extents.toZero();
		center = v;
		axis.toIdentity();
		return true;
	}

	bounds1[0][0] = bounds1[1][0] = center * axis[0];
	bounds1[0][1] = bounds1[1][1] = center * axis[1];
	bounds1[0][2] = bounds1[1][2] = center * axis[2];
	bounds1[0] -= extents;
	bounds1[1] += extents;
	if ( !bounds1.addPoint( CVec3D( v * axis[0], v * axis[1], v * axis[2] ) ) ) {
		// point is contained in the box
		return false;
	}

	axis2[0] = v - center;
	axis2[0].toNormal();
	axis2[1] = axis[ Min3Index( axis2[0] * axis[0], axis2[0] * axis[1], axis2[0] * axis[2] ) ];
	axis2[1] = axis2[1] - ( axis2[1] * axis2[0] ) * axis2[0];
	axis2[1].toNormal();
	axis2[2].cross( axis2[0], axis2[1] );

	projectToAxis( axis2, bounds2 );
	bounds2.addPoint( CVec3D( v * axis2[0], v * axis2[1], v * axis2[2] ) );

	// create new box based on the smallest bounds
	if ( bounds1.getVolume() < bounds2.getVolume() ) {
		center = ( bounds1[0] + bounds1[1] ) * 0.5f;
		extents = bounds1[1] - center;
		center *= axis;
	}
	else {
		center = ( bounds2[0] + bounds2[1] ) * 0.5f;
		extents = bounds2[1] - center;
		center *= axis2;
		axis = axis2;
	}
	return true;
}

/*
============
COBBox::AddBox
============
*/
bool COBBox::AddBox( const COBBox &a ) {
	int i, besti;
	float v, bestv;
	CVec3D dir;
	CMat3D ax[4];
	CAABBox bounds[4], b;

	if ( a.extents[0] < 0.0f ) {
		return false;
	}

	if ( extents[0] < 0.0f ) {
		center = a.center;
		extents = a.extents;
		axis = a.axis;
		return true;
	}

	// test axis of this box
	ax[0] = axis;
	bounds[0][0][0] = bounds[0][1][0] = center * ax[0][0];
	bounds[0][0][1] = bounds[0][1][1] = center * ax[0][1];
	bounds[0][0][2] = bounds[0][1][2] = center * ax[0][2];
	bounds[0][0] -= extents;
	bounds[0][1] += extents;
	a.projectToAxis( ax[0], b );
	if ( !bounds[0].addBounds( b ) ) {
		// the other box is contained in this box
		return false;
	}

	// test axis of other box
	ax[1] = a.axis;
	bounds[1][0][0] = bounds[1][1][0] = a.center * ax[1][0];
	bounds[1][0][1] = bounds[1][1][1] = a.center * ax[1][1];
	bounds[1][0][2] = bounds[1][1][2] = a.center * ax[1][2];
	bounds[1][0] -= a.extents;
	bounds[1][1] += a.extents;
	projectToAxis( ax[1], b );
	if ( !bounds[1].addBounds( b ) ) {
		// this box is contained in the other box
		center = a.center;
		extents = a.extents;
		axis = a.axis;
		return true;
	}

	// test axes aligned with the vector between the box centers and one of the box axis
	dir = a.center - center;
	dir.toNormal();
	for ( i = 2; i < 4; i++ ) {
		ax[i][0] = dir;
		ax[i][1] = ax[i-2][ Min3Index( dir * ax[i-2][0], dir * ax[i-2][1], dir * ax[i-2][2] ) ];
		ax[i][1] = ax[i][1] - ( ax[i][1] * dir ) * dir;
		ax[i][1].toNormal();
		ax[i][2].cross( dir, ax[i][1] );

		projectToAxis( ax[i], bounds[i] );
		a.projectToAxis( ax[i], b );
		bounds[i].addBounds( b );
	}

	// get the bounds with the smallest volume
	bestv = CMath::INFINITY_FLOAT;
	besti = 0;
	for ( i = 0; i < 4; i++ ) {
		v = bounds[i].getVolume();
		if ( v < bestv ) {
			bestv = v;
			besti = i;
		}
	}

	// create a box from the smallest bounds axis pair
	center = ( bounds[besti][0] + bounds[besti][1] ) * 0.5f;
	extents = bounds[besti][1] - center;
	center *= ax[besti];
	axis = ax[besti];

	return false;
}

/*
================
COBBox::planeDistance
================
*/
float COBBox::planeDistance( const CPlane &plane ) const {
	float d1, d2;

	d1 = plane.getDistance( center );
	d2 = CMath::fabs( extents[0] * plane.getNormal()[0] ) +
			CMath::fabs( extents[1] * plane.getNormal()[1] ) +
				CMath::fabs( extents[2] * plane.getNormal()[2] );

	if ( d1 - d2 > 0.0f ) {
		return d1 - d2;
	}
	if ( d1 + d2 < 0.0f ) {
		return d1 + d2;
	}
	return 0.0f;
}

/*
================
COBBox::planeSide
================
*/
int COBBox::planeSide( const CPlane &plane, const float epsilon ) const {
	float d1, d2;

	d1 = plane.getDistance( center );
	d2 = CMath::fabs( extents[0] * plane.getNormal()[0] ) +
			CMath::fabs( extents[1] * plane.getNormal()[1] ) +
				CMath::fabs( extents[2] * plane.getNormal()[2] );

	if ( d1 - d2 > epsilon ) {
		return PLANESIDE_FRONT;
	}
	if ( d1 + d2 < -epsilon ) {
		return PLANESIDE_BACK;
	}
	return PLANESIDE_CROSS;
}

/*
============
COBBox::intersectsBox
============
*/
bool COBBox::intersectsBox( const COBBox &a ) const {
    CVec3D dir;			// vector between centers
    float c[3][3];		// matrix c = axis.transpose() * a.axis
    float ac[3][3];		// absolute values of c
    float axisdir[3];	// axis[i] * dir
    float d, e0, e1;	// distance between centers and projected extents

	dir = a.center - center;
    
    // axis C0 + t * A0
    c[0][0] = axis[0] * a.axis[0];
    c[0][1] = axis[0] * a.axis[1];
    c[0][2] = axis[0] * a.axis[2];
    axisdir[0] = axis[0] * dir;
    ac[0][0] = CMath::fabs( c[0][0] );
    ac[0][1] = CMath::fabs( c[0][1] );
    ac[0][2] = CMath::fabs( c[0][2] );

    d = CMath::fabs( axisdir[0] );
	e0 = extents[0];
    e1 = a.extents[0] * ac[0][0] + a.extents[1] * ac[0][1] + a.extents[2] * ac[0][2];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A1
    c[1][0] = axis[1] * a.axis[0];
    c[1][1] = axis[1] * a.axis[1];
    c[1][2] = axis[1] * a.axis[2];
    axisdir[1] = axis[1] * dir;
    ac[1][0] = CMath::fabs( c[1][0] );
    ac[1][1] = CMath::fabs( c[1][1] );
    ac[1][2] = CMath::fabs( c[1][2] );

    d = CMath::fabs( axisdir[1] );
	e0 = extents[1];
    e1 = a.extents[0] * ac[1][0] + a.extents[1] * ac[1][1] + a.extents[2] * ac[1][2];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A2
    c[2][0] = axis[2] * a.axis[0];
    c[2][1] = axis[2] * a.axis[1];
    c[2][2] = axis[2] * a.axis[2];
    axisdir[2] = axis[2] * dir;
    ac[2][0] = CMath::fabs( c[2][0] );
    ac[2][1] = CMath::fabs( c[2][1] );
    ac[2][2] = CMath::fabs( c[2][2] );

    d = CMath::fabs( axisdir[2] );
	e0 = extents[2];
    e1 = a.extents[0] * ac[2][0] + a.extents[1] * ac[2][1] + a.extents[2] * ac[2][2];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * B0
    d = CMath::fabs( a.axis[0] * dir );
    e0 = extents[0] * ac[0][0] + extents[1] * ac[1][0] + extents[2] * ac[2][0];
	e1 = a.extents[0];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * B1
    d = CMath::fabs( a.axis[1] * dir );
    e0 = extents[0] * ac[0][1] + extents[1] * ac[1][1] + extents[2] * ac[2][1];
	e1 = a.extents[1];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * B2
    d = CMath::fabs( a.axis[2] * dir );
    e0 = extents[0] * ac[0][2] + extents[1] * ac[1][2] + extents[2] * ac[2][2];
	e1 = a.extents[2];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A0xB0
    d = CMath::fabs( axisdir[2] * c[1][0] - axisdir[1] * c[2][0] );
    e0 = extents[1] * ac[2][0] + extents[2] * ac[1][0];
    e1 = a.extents[1] * ac[0][2] + a.extents[2] * ac[0][1];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A0xB1
    d = CMath::fabs( axisdir[2] * c[1][1] - axisdir[1] * c[2][1] );
    e0 = extents[1] * ac[2][1] + extents[2] * ac[1][1];
    e1 = a.extents[0] * ac[0][2] + a.extents[2] * ac[0][0];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A0xB2
    d = CMath::fabs( axisdir[2] * c[1][2] - axisdir[1] * c[2][2] );
    e0 = extents[1] * ac[2][2] + extents[2] * ac[1][2];
    e1 = a.extents[0] * ac[0][1] + a.extents[1] * ac[0][0];
    if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A1xB0
    d = CMath::fabs( axisdir[0] * c[2][0] - axisdir[2] * c[0][0] );
    e0 = extents[0] * ac[2][0] + extents[2] * ac[0][0];
    e1 = a.extents[1] * ac[1][2] + a.extents[2] * ac[1][1];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A1xB1
    d = CMath::fabs( axisdir[0] * c[2][1] - axisdir[2] * c[0][1] );
    e0 = extents[0] * ac[2][1] + extents[2] * ac[0][1];
    e1 = a.extents[0] * ac[1][2] + a.extents[2] * ac[1][0];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A1xB2
    d = CMath::fabs( axisdir[0] * c[2][2] - axisdir[2] * c[0][2] );
    e0 = extents[0] * ac[2][2] + extents[2] * ac[0][2];
    e1 = a.extents[0] * ac[1][1] + a.extents[1] * ac[1][0];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A2xB0
    d = CMath::fabs( axisdir[1] * c[0][0] - axisdir[0] * c[1][0] );
    e0 = extents[0] * ac[1][0] + extents[1] * ac[0][0];
    e1 = a.extents[1] * ac[2][2] + a.extents[2] * ac[2][1];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A2xB1
    d = CMath::fabs( axisdir[1] * c[0][1] - axisdir[0] * c[1][1] );
    e0 = extents[0] * ac[1][1] + extents[1] * ac[0][1];
    e1 = a.extents[0] * ac[2][2] + a.extents[2] * ac[2][0];
	if ( d > e0 + e1 ) {
        return false;
	}

    // axis C0 + t * A2xB2
    d = CMath::fabs( axisdir[1] * c[0][2] - axisdir[0] * c[1][2] );
    e0 = extents[0] * ac[1][2] + extents[1] * ac[0][2];
    e1 = a.extents[0] * ac[2][1] + a.extents[1] * ac[2][0];
	if ( d > e0 + e1 ) {
        return false;
	}
    return true;
}

/*
============
COBBox::lineIntersection

  Returns true if the line intersects the box between the start and end point.
============
*/
bool COBBox::lineIntersection( const CVec3D &start, const CVec3D &end ) const {
    float ld[3];
    CVec3D lineDir = 0.5f * ( end - start );
    CVec3D lineCenter = start + lineDir;
    CVec3D dir = lineCenter - center;

    ld[0] = CMath::fabs( lineDir * axis[0] );
	if ( CMath::fabs( dir * axis[0] ) > extents[0] + ld[0] ) {
        return false;
	}

    ld[1] = CMath::fabs( lineDir * axis[1] );
	if ( CMath::fabs( dir * axis[1] ) > extents[1] + ld[1] ) {
        return false;
	}

    ld[2] = CMath::fabs( lineDir * axis[2] );
	if ( CMath::fabs( dir * axis[2] ) > extents[2] + ld[2] ) {
        return false;
	}

    CVec3D cross = lineDir.cross( dir );

	if ( CMath::fabs( cross * axis[0] ) > extents[1] * ld[2] + extents[2] * ld[1] ) {
        return false;
	}

	if ( CMath::fabs( cross * axis[1] ) > extents[0] * ld[2] + extents[2] * ld[0] ) {
        return false;
	}

	if ( CMath::fabs( cross * axis[2] ) > extents[0] * ld[1] + extents[1] * ld[0] ) {
        return false;
	}

    return true;
}

/*
============
BoxPlaneClip
============
*/
static bool BoxPlaneClip( const float denom, const float numer, float &scale0, float &scale1 ) {
	if ( denom > 0.0f ) {
		if ( numer > denom * scale1 ) {
			return false;
		}
		if ( numer > denom * scale0 ) {
			scale0 = numer / denom;
		}
		return true;
	}
	else if ( denom < 0.0f ) {
		if ( numer > denom * scale0 ) {
			return false;
		}
		if ( numer > denom * scale1 ) {
			scale1 = numer / denom;
		}
		return true;
	}
	else {
		return ( numer <= 0.0f );
	}
}

/*
============
COBBox::rayIntersection

  Returns true if the ray intersects the box.
  The ray can intersect the box in both directions from the start point.
  If start is inside the box then scale1 < 0 and scale2 > 0.
============
*/
bool COBBox::rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale1, float &scale2 ) const {
	CVec3D localStart, localDir;

	localStart = ( start - center ) * axis.transpose();
	localDir = dir * axis.transpose();

	scale1 = -CMath::INFINITY_FLOAT;
	scale2 = CMath::INFINITY_FLOAT;
    return	BoxPlaneClip(  localDir.x, -localStart.x - extents[0], scale1, scale2 ) &&
			BoxPlaneClip( -localDir.x,  localStart.x - extents[0], scale1, scale2 ) &&
			BoxPlaneClip(  localDir.y, -localStart.y - extents[1], scale1, scale2 ) &&
			BoxPlaneClip( -localDir.y,  localStart.y - extents[1], scale1, scale2 ) &&
			BoxPlaneClip(  localDir.z, -localStart.z - extents[2], scale1, scale2 ) &&
			BoxPlaneClip( -localDir.z,  localStart.z - extents[2], scale1, scale2 );
}

/*
============
COBBox::fromPoints

  Tight box for a collection of points.
============
*/
void COBBox::fromPoints( const CVec3D *points, const int numPoints ) {
	int i;
	float invNumPoints, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ;
	CVec3D dir;
	CAABBox bounds;
	CMatXD eigenVectors;
	CVecXD eigenValues;

	// compute mean of points
	center = points[0];
	for ( i = 1; i < numPoints; i++ ) {
		center += points[i];
	}
	invNumPoints = 1.0f / numPoints;
	center *= invNumPoints;

	// compute covariances of points
	sumXX = 0.0f; sumXY = 0.0f; sumXZ = 0.0f;
	sumYY = 0.0f; sumYZ = 0.0f; sumZZ = 0.0f;
	for ( i = 0; i < numPoints; i++ ) {
		dir = points[i] - center;
		sumXX += dir.x * dir.x;
		sumXY += dir.x * dir.y;
		sumXZ += dir.x * dir.z;
		sumYY += dir.y * dir.y;
		sumYZ += dir.y * dir.z;
		sumZZ += dir.z * dir.z;
	}
	sumXX *= invNumPoints;
	sumXY *= invNumPoints;
	sumXZ *= invNumPoints;
	sumYY *= invNumPoints;
	sumYZ *= invNumPoints;
	sumZZ *= invNumPoints;

	// compute eigenvectors for covariance matrix
	eigenValues.setData( 3, VECX_ALLOCA( 3 ) );
	eigenVectors.setData( 3, 3, MATX_ALLOCA( 9 ) ); //(3 * 3)

	eigenVectors[0][0] = sumXX;
	eigenVectors[0][1] = sumXY;
	eigenVectors[0][2] = sumXZ;
	eigenVectors[1][0] = sumXY;
	eigenVectors[1][1] = sumYY;
	eigenVectors[1][2] = sumYZ;
	eigenVectors[2][0] = sumXZ;
	eigenVectors[2][1] = sumYZ;
	eigenVectors[2][2] = sumZZ;
	eigenVectors.eigen_solveSymmetric( eigenValues );
	eigenVectors.eigen_SortIncreasing( eigenValues );

	axis[0][0] = eigenVectors[0][0];
	axis[0][1] = eigenVectors[0][1];
	axis[0][2] = eigenVectors[0][2];
	axis[1][0] = eigenVectors[1][0];
	axis[1][1] = eigenVectors[1][1];
	axis[1][2] = eigenVectors[1][2];
	axis[2][0] = eigenVectors[2][0];
	axis[2][1] = eigenVectors[2][1];
	axis[2][2] = eigenVectors[2][2];

	extents[0] = eigenValues[0];
	extents[1] = eigenValues[0];
	extents[2] = eigenValues[0];

	// refine by calculating the bounds of the points projected onto the axis and adjusting the center and extents
	bounds.clear();
    for ( i = 0; i < numPoints; i++ ) {
		bounds.addPoint( CVec3D( points[i] * axis[0], points[i] * axis[1], points[i] * axis[2] ) );
    }
	center = ( bounds[0] + bounds[1] ) * 0.5f;
	extents = bounds[1] - center;
	center *= axis;
}

/*
============
COBBox::fromPointTranslation

  Most tight box for the translational movement of the given point.
============
*/
void COBBox::fromPointTranslation( const CVec3D &point, const CVec3D &translation ) {
	// FIXME: implement
}

/*
============
COBBox::fromBoxTranslation

  Most tight box for the translational movement of the given box.
============
*/
void COBBox::fromBoxTranslation( const COBBox &box, const CVec3D &translation ) {
	// FIXME: implement
}
#if 0
/*
============
COBBox::fromPointRotation

  Most tight bounds for the rotational movement of the given point.
============
*/
void COBBox::fromPointRotation( const CVec3D &point, const CRotation &rotation ) {
	// FIXME: implement
}

/*
============
COBBox::FromBoxRotation

  Most tight box for the rotational movement of the given box.
============
*/
void COBBox::FromBoxRotation( const COBBox &box, const CRotation &rotation ) {
	// FIXME: implement
}
#endif 
/*
============
COBBox::toPoints
============
*/
void COBBox::toPoints( CVec3D points[8] ) const {
	CMat3D ax;
	CVec3D temp[4];

	ax[0] = extents[0] * axis[0];
	ax[1] = extents[1] * axis[1];
	ax[2] = extents[2] * axis[2];
	temp[0] = center - ax[0];
	temp[1] = center + ax[0];
	temp[2] = ax[1] - ax[2];
	temp[3] = ax[1] + ax[2];
	points[0] = temp[0] - temp[3];
	points[1] = temp[1] - temp[3];
	points[2] = temp[1] + temp[2];
	points[3] = temp[0] + temp[2];
	points[4] = temp[0] - temp[2];
	points[5] = temp[1] - temp[2];
	points[6] = temp[1] + temp[3];
	points[7] = temp[0] + temp[3];
}

/*
============
COBBox::getProjectionSilhouetteVerts
============
*/
int COBBox::getProjectionSilhouetteVerts( const CVec3D &projectionOrigin, CVec3D silVerts[6] ) const {
	float f;
	int i, planeBits, *index;
	CVec3D points[8], dir1, dir2;

	toPoints( points );

	dir1 = points[0] - projectionOrigin;
	dir2 = points[6] - projectionOrigin;
	f = dir1 * axis[0];
	planeBits = IEEE_FLT_SIGNBITNOTSET( f );
	f = dir2 * axis[0];
	planeBits |= IEEE_FLT_SIGNBITSET( f ) << 1;
	f = dir1 * axis[1];
	planeBits |= IEEE_FLT_SIGNBITNOTSET( f ) << 2;
	f = dir2 * axis[1];
	planeBits |= IEEE_FLT_SIGNBITSET( f ) << 3;
	f = dir1 * axis[2];
	planeBits |= IEEE_FLT_SIGNBITNOTSET( f ) << 4;
	f = dir2 * axis[2];
	planeBits |= IEEE_FLT_SIGNBITSET( f ) << 5;

	index = boxPlaneBitsSilVerts[planeBits];
	for ( i = 0; i < index[0]; i++ ) {
		silVerts[i] = points[index[i+1]];
	}

	return index[0];
}

/*
============
COBBox::getParallelProjectionSilhouetteVerts
============
*/
int COBBox::getParallelProjectionSilhouetteVerts( const CVec3D &projectionDir, CVec3D silVerts[6] ) const {
	float f;
	int i, planeBits, *index;
	CVec3D points[8];

	toPoints( points );

	planeBits = 0;
	f = projectionDir * axis[0];
	if ( IEEE_FLT_ISNOTZERO( f ) ) {
		planeBits = 1 << IEEE_FLT_SIGNBITSET( f );
	}
	f = projectionDir * axis[1];
	if ( IEEE_FLT_ISNOTZERO( f ) ) {
		planeBits |= 4 << IEEE_FLT_SIGNBITSET( f );
	}
	f = projectionDir * axis[2];
	if ( IEEE_FLT_ISNOTZERO( f ) ) {
		planeBits |= 16 << IEEE_FLT_SIGNBITSET( f );
	}

	index = boxPlaneBitsSilVerts[planeBits];
	for ( i = 0; i < index[0]; i++ ) {
		silVerts[i] = points[index[i+1]];
	}

	return index[0];
}

/// The implementation of this function is from Christer Ericson's Real-Time Collision Detection, p.133.
CVec3D COBBox::closestPoint(const CVec3D &targetPoint) const
{
	CVec3D d = targetPoint - center;
	CVec3D closestPoint = center; // Start at the center point of the COBBox.
	for(int i = 0; i < 3; ++i) // project the target onto the COBBox axes and walk towards that point.
		closestPoint += CLAMP(d*axis[i], -extents[i], extents[i]) * axis[i];

	return closestPoint;
}

float COBBox::distance(const CVec3D &point) const
{
	///\todo This code can be optimized a bit. See Christer Ericson's Real-Time Collision Detection,
	/// p.134.
	CVec3D closestPoint_ = closestPoint(point);
	return point.distance(closestPoint_);
}

float COBBox::distance(const CSphere &sphere) const
{
	return MAX(0.f, distance(sphere.getOrigin()) - sphere.getRadius());
}



bool COBBox::intersects(const CAABBox &aabb) const
{
	return intersects(COBBox(aabb));
}


bool COBBox::intersects(const COBBox &b, float epsilon) const
{
	SMF_ASSERT(center.isFinite());
	SMF_ASSERT(b.center.isFinite());
	SMF_ASSERT(CVec3D::areOrthonormal(axis[0], axis[1], axis[2]));
	SMF_ASSERT(CVec3D::areOrthonormal(b.axis[0], b.axis[1], b.axis[2]));

	// Generate a rotation matrix that transforms from world space to this COBBox's coordinate space.
	CMat3D R;
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			R[i][j] = (axis[i]* b.axis[j]);

	CVec3D t = b.center - center;
	// Express the translation vector in a's coordinate frame.
	t = CVec3D((t* axis[0]), (t*axis[1]), (t* axis[2]));

	CMat3D AbsR;
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			AbsR[i][j] = CMath::fabs(R[i][j]) + epsilon;

	// Test the three major axes of this COBBox.
	for(int i = 0; i < 3; ++i)
	{
		float ra = extents[i];
		float rb = ((b.extents)* AbsR[i]);
		if (CMath::fabs(t[i]) > ra + rb)
			return false;
	}

	// Test the three major axes of the COBBox b.
	for(int i = 0; i < 3; ++i)
	{
		float ra = extents[0] * AbsR[0][i] + extents[1] * AbsR[1][i] + extents[2] * AbsR[2][i];
		float rb = b.extents[i];
		if (CMath::fabs(t.x * R[0][i] + t.y * R[1][i] + t.z * R[2][i]) > ra + rb)
			return false;
	}

	// Test the 9 different cross-axes.

	// A.x <cross> B.x
	float ra = extents.y * AbsR[2][0] + extents.z * AbsR[1][0];
	float rb = b.extents.y * AbsR[0][2] + b.extents.z * AbsR[0][1];
	if (CMath::fabs(t.z * R[1][0] - t.y * R[2][0]) > ra + rb)
		return false;

	// A.x < cross> B.y
	ra = extents.y * AbsR[2][1] + extents.z * AbsR[1][1];
	rb = b.extents.x * AbsR[0][2] + b.extents.z * AbsR[0][0];
	if (CMath::fabs(t.z * R[1][1] - t.y * R[2][1]) > ra + rb)
		return false;

	// A.x <cross> B.z
	ra = extents.y * AbsR[2][2] + extents.z * AbsR[1][2];
	rb = b.extents.x * AbsR[0][1] + b.extents.y * AbsR[0][0];
	if (CMath::fabs(t.z * R[1][2] - t.y * R[2][2]) > ra + rb)
		return false;

	// A.y <cross> B.x
	ra = extents.x * AbsR[2][0] + extents.z * AbsR[0][0];
	rb = b.extents.y * AbsR[1][2] + b.extents.z * AbsR[1][1];
	if (CMath::fabs(t.x * R[2][0] - t.z * R[0][0]) > ra + rb)
		return false;

	// A.y <cross> B.y
	ra = extents.x * AbsR[2][1] + extents.z * AbsR[0][1];
	rb = b.extents.x * AbsR[1][2] + b.extents.z * AbsR[1][0];
	if (CMath::fabs(t.x * R[2][1] - t.z * R[0][1]) > ra + rb)
		return false;

	// A.y <cross> B.z
	ra = extents.x * AbsR[2][2] + extents.z * AbsR[0][2];
	rb = b.extents.x * AbsR[1][1] + b.extents.y * AbsR[1][0];
	if (CMath::fabs(t.x * R[2][2] - t.z * R[0][2]) > ra + rb)
		return false;

	// A.z <cross> B.x
	ra = extents.x * AbsR[1][0] + extents.y * AbsR[0][0];
	rb = b.extents.y * AbsR[2][2] + b.extents.z * AbsR[2][1];
	if (CMath::fabs(t.y * R[0][0] - t.x * R[1][0]) > ra + rb)
		return false;

	// A.z <cross> B.y
	ra = extents.x * AbsR[1][1] + extents.y * AbsR[0][1];
	rb = b.extents.x * AbsR[2][2] + b.extents.z * AbsR[2][0];
	if (CMath::fabs(t.y * R[0][1] - t.x * R[1][1]) > ra + rb)
		return false;

	// A.z <cross> B.z
	ra = extents.x * AbsR[1][2] + extents.y * AbsR[0][2];
	rb = b.extents.x * AbsR[2][1] + b.extents.y * AbsR[2][0];
	if (CMath::fabs(t.y * R[0][2] - t.x * R[1][2]) > ra + rb)
		return false;

	// No separating axis exists, so the two COBBox don't intersect.
	return true;
}

/// The implementation of COBBox-Plane intersection test follows Christer Ericson's Real-Time Collision Detection, p. 163. [groupSyntax]
bool COBBox::intersects(const CPlane &p) const
{
	// Compute the projection interval radius of this COBBox onto L(t) = this->center + x * p.getNormal();
	float t = extents[0] * CMath::fabs((p.getNormal()* axis[0])) +
			  extents[1] * CMath::fabs((p.getNormal()* axis[1])) +
			  extents[2] * CMath::fabs((p.getNormal()* axis[2]));
	// Compute the distance of this COBBox center from the plane.
	float s = (p.getNormal()* center) - p.d;
	return CMath::fabs(s) <= t;
}

bool COBBox::intersects(const CRay &ray) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CRay r = worldToLocal() * ray;
	return aabb.intersects(r);
}

bool COBBox::intersects(const CRay &ray, float &dNear, float &dFar) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CRay r = worldToLocal() * ray;
	return aabb.intersects(r, dNear, dFar);
}

bool COBBox::intersects(const CLine &line) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CLine l = worldToLocal() * line;
	return aabb.intersects(l);
}

bool COBBox::intersects(const CLine &line, float &dNear, float &dFar) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CLine l = worldToLocal() * line;
	return aabb.intersects(l, dNear, dFar);
}

bool COBBox::intersects(const CLineSegment &lineSegment) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CLineSegment l = worldToLocal() * lineSegment;
	return aabb.intersects(l);
}

bool COBBox::intersects(const CLineSegment &lineSegment, float &dNear, float &dFar) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CLineSegment l = worldToLocal() * lineSegment;
	return aabb.intersects(l, dNear, dFar);
}

/// The implementation of the COBBox-Sphere intersection test follows Christer Ericson's Real-Time Collision Detection, p. 166. [groupSyntax]
bool COBBox::intersects(const CSphere &sphere, CVec3D *closestPointOnOBB) const
{
	// find the point on this CAABBox closest to the sphere center.
	CVec3D pt = closestPoint(sphere.getOrigin());

	// If that point is inside sphere, the CAABBox and sphere intersect.
	if (closestPointOnOBB)
		*closestPointOnOBB = pt;

	return pt.distanceSq(sphere.getOrigin()) <= sphere.getRadius() * sphere.getRadius();
}



bool COBBox::intersects(const CTriangle &triangle) const
{
	CAABBox aabb(CVec3D(0,0,0), CVec3D(size()));
	CTriangle t = worldToLocal() * triangle;
	return t.intersects(aabb);
}

bool COBBox::intersects(const CPolygon &polygon) const
{
	return toPolyhedron().intersects(polygon);
}



bool COBBox::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

CVec3D COBBox::size() const
{
	return extents * 2.f;
}

CVec3D COBBox::halfSize() const
{
	return extents;
}

CVec3D COBBox::diagonal() const
{
	return 2.f * halfDiagonal();
}

CVec3D COBBox::halfDiagonal() const
{
	return axis[0] * extents[0] + axis[1] * extents[1] + axis[2] * extents[2];
}

CMatJoint3x4 COBBox::worldToLocal() const
{
	CMatJoint3x4 m = localToWorld();
	m.inverseOrthonormal();
	return m;
}

CMatJoint3x4 COBBox::localToWorld() const
{
	// To produce a normalized local->world matrix, do the following.
	/*
	CMatJoint3x4 m;
	CVec3D x = axis[0] * extents.x;
	CVec3D y = axis[1] * extents.y;
	CVec3D z = axis[2] * extents.z;
	m.setCol(0, 2.f * x);
	m.setCol(1, 2.f * y);
	m.setCol(2, 2.f * z);
	m.setCol(3, center - x - y - z);
	return m;
	*/

	SMF_ASSERT(axis[0].isNormalized());
	SMF_ASSERT(axis[1].isNormalized());
	SMF_ASSERT(axis[2].isNormalized());
	CMatJoint3x4 m;
	m.setCol(0, axis[0]);
	m.setCol(1, axis[1]);
	m.setCol(2, axis[2]);
	m.setCol(3, center - axis[0] * extents.x - axis[1] * extents.y - axis[2] * extents.z);
	SMF_ASSERT(m.isOrthonormal());
	return m;
}

CAABBox COBBox::minimalEnclosingAABB() const
{
	CAABBox aabb;
	aabb.setFrom(*this);
	return aabb;
}

#if 0

CAABBox COBBox::maximalContainedAABB() const
{
#ifdef _MSC_VER
#pragma WARNING(COBBox::maximalContainedAABB not implemented!)
#else
#warning COBBox::maximalContainedAABB not implemented!
#endif
	SMF_ASSERT(false && "COBBox::maximalContainedAABB not implemented!"); /// \todo Implement.
	return CAABBox();
}
#endif

CSphere COBBox::minimalEnclosingSphere() const
{
	CSphere s(center, halfDiagonal().getLenght());
	
	return s;
}

CSphere COBBox::maximalContainedSphere() const
{
	CSphere s(center, extents.minElement());
	return s;
}

CLineSegment COBBox::edge(int edgeIndex) const
{
	SMF_ASSERT(0 <= edgeIndex && edgeIndex <= 11);
	switch(edgeIndex)
	{
		default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
		case 0: return CLineSegment(cornerPoint(0), cornerPoint(1));
		case 1: return CLineSegment(cornerPoint(0), cornerPoint(2));
		case 2: return CLineSegment(cornerPoint(0), cornerPoint(4));
		case 3: return CLineSegment(cornerPoint(1), cornerPoint(3));
		case 4: return CLineSegment(cornerPoint(1), cornerPoint(5));
		case 5: return CLineSegment(cornerPoint(2), cornerPoint(3));
		case 6: return CLineSegment(cornerPoint(2), cornerPoint(6));
		case 7: return CLineSegment(cornerPoint(3), cornerPoint(7));
		case 8: return CLineSegment(cornerPoint(4), cornerPoint(5));
		case 9: return CLineSegment(cornerPoint(4), cornerPoint(6));
		case 10: return CLineSegment(cornerPoint(5), cornerPoint(7));
		case 11: return CLineSegment(cornerPoint(6), cornerPoint(7));
	}
}

CVec3D COBBox::cornerPoint(int cornerIndex) const
{	
	SMF_ASSERT(0 <= cornerIndex && cornerIndex <= 7);
	switch(cornerIndex)
	{
		default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
		case 0: return center - extents.x * axis[0] - extents.y * axis[1] - extents.z * axis[2];
		case 1: return center - extents.x * axis[0] - extents.y * axis[1] + extents.z * axis[2];
		case 2: return center - extents.x * axis[0] + extents.y * axis[1] - extents.z * axis[2];
		case 3: return center - extents.x * axis[0] + extents.y * axis[1] + extents.z * axis[2];
		case 4: return center + extents.x * axis[0] - extents.y * axis[1] - extents.z * axis[2];
		case 5: return center + extents.x * axis[0] - extents.y * axis[1] + extents.z * axis[2];
		case 6: return center + extents.x * axis[0] + extents.y * axis[1] - extents.z * axis[2];
		case 7: return center + extents.x * axis[0] + extents.y * axis[1] + extents.z * axis[2];
	}
}

CVec3D COBBox::extremePoint(const CVec3D &direction) const
{
	CVec3D pt = center;
	pt += axis[0] * ((direction* axis[0]) >= 0.f ? extents.x : -extents.x);
	pt += axis[1] * ((direction* axis[1]) >= 0.f ? extents.y : -extents.y);
	pt += axis[2] * ((direction* axis[2]) >= 0.f ? extents.z : -extents.z);
	return pt;
}

CPolyhedron COBBox::toPolyhedron() const
{
	// Note to maintainer: This function is an exact copy of AABB:toPolyhedron() and Frustum::toPolyhedron().

	CPolyhedron p;
	// Populate the corners of this COBBox.
	// The will be in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++.
	for(int i = 0; i < 8; ++i)
		p.v.push_back(cornerPoint(i));

	// Generate the 6 faces of this COBBox.
	const int faces[6][4] =
	{
		{ 0, 1, 3, 2 }, // X-
		{ 4, 6, 7, 5 }, // X+
		{ 0, 4, 5, 1 }, // Y-
		{ 7, 6, 2, 3 }, // Y+
		{ 0, 2, 6, 4 }, // Z-
		{ 1, 5, 7, 3 }, // Z+
	};

	for(int f = 0; f < 6; ++f)
	{
		CPolyhedron::Face face;
		for(int v = 0; v < 4; ++v)
			face.v.push_back(faces[f][v]);
		p.f.push_back(face);
	}

	return p;
}

/// See Christer Ericson's book Real-Time Collision Detection, page 83.
void COBBox::extremePointsAlongDirection(const CVec3D &dir, const CVec3D *pointArray, int numPoints, int &idxSmallest, int &idxLargest)
{
	SMF_ASSERT(pointArray || numPoints == 0);

	idxSmallest = idxLargest = 0;

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!pointArray)
		return;
#endif

	float smallestD = CMath::INFINITY_FLOAT;
	float largestD = -CMath::INFINITY_FLOAT;
	for(int i = 0; i < numPoints; ++i)
	{
		float d = (pointArray[i]* dir);
		if (d < smallestD)
		{
			smallestD = d;
			idxSmallest = i;
		}
		if (d > largestD)
		{
			largestD = d;
			idxLargest = i;
		}
	}
}

bool COBBox::contains(const CVec3D &point) const
{
	CVec3D pt = point - center;
	return CMath::fabs((pt* axis[0])) <= extents[0] &&
	       CMath::fabs((pt* axis[1])) <= extents[1] &&
	       CMath::fabs((pt* axis[2])) <= extents[2];
}

bool COBBox::contains(const CLineSegment &lineSegment) const
{
	return contains(lineSegment.begin) && contains(lineSegment.end);
}

bool COBBox::contains(const CAABBox &aabb) const
{
	// Since both AABB and COBBox are convex objects, this COBBox contains the AABB
	// if and only if it contains all its corner points.
	for(int i = 0; i < 8; ++i)
		if (!contains(aabb.cornerPoint(i)))
			return false;

	return true;
}

bool COBBox::contains(const COBBox &obb) const
{
	for(int i = 0; i < 8; ++i)
		if (!contains(obb.cornerPoint(i)))
			return false;

	return true;
}

bool COBBox::contains(const CTriangle &triangle) const
{
	return contains(triangle.a) && contains(triangle.b) && contains(triangle.c);
}

bool COBBox::contains(const CPolygon &polygon) const
{
	for(int i = 0; i < polygon.numVertices(); ++i)
		if (!contains(polygon.vertex(i)))
			return false;
	return true;
}



bool COBBox::contains(const CPolyhedron &polyhedron) const
{
	SMF_ASSERT(polyhedron.IsClosed());
	for(int i = 0; i < polyhedron.numVertices(); ++i)
		if (!contains(polyhedron.vertex(i)))
			return false;

	return true;
}

CVec3D COBBox::faceCenterPoint(int faceIndex) const
{
	//assume(0 <= faceIndex && faceIndex <= 5);

	switch(faceIndex)
	{
	default: // For release builds where assume() is disabled, return always the first option if out-of-bounds.
	case 0: return center - extents.x * axis[0];
	case 1: return center + extents.x * axis[0];
	case 2: return center - extents.y * axis[1];
	case 3: return center + extents.y * axis[1];
	case 4: return center - extents.z * axis[2];
	case 5: return center + extents.z * axis[2];
	}
}

CPlane COBBox::facePlane(int faceIndex) const
{
	//assume(0 <= faceIndex && faceIndex <= 5);
	switch(faceIndex)
	{
	default: // For release builds where assume() is disabled, return always the first option if out-of-bounds.
	case 0: return CPlane(faceCenterPoint(0), -axis[0]);
	case 1: return CPlane(faceCenterPoint(1), axis[0]);
	case 2: return CPlane(faceCenterPoint(2), -axis[1]);
	case 3: return CPlane(faceCenterPoint(3), axis[1]);
	case 4: return CPlane(faceCenterPoint(4), -axis[2]);
	case 5: return CPlane(faceCenterPoint(5), axis[2]);
	}
}

} //end GEO
}  //end SMF