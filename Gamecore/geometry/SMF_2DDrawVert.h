/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#ifndef _SMF__DRAWVERT_2D_2D_H__
#define _SMF__DRAWVERT_2D_2D_H__
#include "../SMF_Config.h"
#include "../math/SMF_Vector.h"
#include "../math/SMF_Simd.h"
#include "../math/SMF_JointTransform.h"

namespace SMF{
using namespace MATH;
namespace GEO {
/*
===============================================================================

	Draw vertex.
	A vertex is the smallest graphical entity that you can manipulate. In short, it is a graphical point: it contains of course a 2D position (x, y), but also a color
===============================================================================
*/


/**
 * \class CVertex2D
 *
 * \ingroup SMF_Geometric
 * \if pt_br
 * \brief Implementa um vertex 3D, con sua tangente, bitangente, normal e coordenada e textura
 * \elseif us_en
 * \brief Draw vertex 3D
 * \note Tangents and Normals
	\p With a 3-D surface you can have a tangent line or a tangent plane. Either of these touches the surface at a single point rather than cutting through it or missing it entirely. Think of a straight piece of wire (tangent line) or a flat piece of cardboard (tangent plane) resting on top of a soccer ball.
	\p The normal to a surface is a only a line (you can't have a normal plane). The normal is a line which is perpendicular to the tangent plane at the same point. By the way, a line is perpendicular to a plane if it is perpendicular to all lines in the plane through the point of contact. 
* \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 * 
  */
class SMF_API CVertex2D {
public:
	/**
	\brief localização do vertex no espaço
	**/
	MATH::CVec2D			xy;  //8 bytes
	sf_u8					color[4];     // 4bytes
	sf_u8					color2[4];		// 4 bytes -- weights for skinning


	float				operator[]( const int index ) const;
	float &				operator[]( const int index );

	void				clear();

	const CVec2D		getNormal() const;
	const CVec2D		getNormalRaw() const;		// not re-normalized for renderbump
	//smoothly transition
	void				lerp( const CVertex2D &a, const CVertex2D &b, const float f );
	void				lerpAll( const CVertex2D &a, const CVertex2D &b, const float f );

	void				setColor( sf_u32 color );
	void				setNativeOrderColor( sf_u32 color );
	sf_u32				getColor() const;

	void				setColor2( sf_u32 color );
	void				setNativeOrderColor2( sf_u32 color );
	void				clearColor2();
	sf_u32				getColor2() const;

};

#define DRAWVERT_2D_SIZE				16
#define DRAWVERT_2D_SIZE_STR			"16"

#define DRAWVERT_2D_XYZ_OFFSET			(0*4)
#define DRAWVERT_2D_XYZ_OFFSET_STR		"(0*4)"
#define DRAWVERT_2D_COLOR_OFFSET		(1*8)
#define DRAWVERT_2D_COLOR2_OFFSET		(1*8)

//assert_offsetof( CVertex2D, xy,		DRAWVERT_2D_XYZ_OFFSET );
//assert_offsetof( CVertex2D, normal,	DRAWVERT_2D_NORMAL_OFFSET );
//SMF_ASSERT( (int)&((CVertex2D *)0)->tangents[0] == DRAWVERT_2D_TANGENT0_OFFSET );
//SMF_ASSERT( (int)&((CVertex2D *)0)->tangents[1] == DRAWVERT_2D_TANGENT1_OFFSET );
//assert_offsetof( CVertex2D, tangent,	DRAWVERT_2D_TANGENT_OFFSET );

/*
========================
VertexFloatToByte

Assumes input is in the range [-1, 1]
========================
*/
SMF_INLINE void VertexFloatToByte( const float & x, const float & y, const float & z, sf_u8 * bval ) {
	assert_4_byte_aligned( bval );	// for __stvebx



/*
========================
CVertex2D::operator[]
========================
*/
SMF_INLINE float CVertex2D::operator[]( const int index ) const {
	SMF_ASSERT( index >= 0 && index < 5 );
	return ((float *)(&xy))[index];
}

/*
========================
CVertex2D::operator[]
========================
*/
SMF_INLINE float	&CVertex2D::operator[]( const int index ) {
	SMF_ASSERT( index >= 0 && index < 5 );
	return ((float *)(&xy))[index];
}

/*
========================
CVertex2D::clear
========================
*/
SMF_INLINE void CVertex2D::clear() {
	*reinterpret_cast<sf_u32 *>(&this->xy.x) = 0;
	*reinterpret_cast<sf_u32 *>(&this->xy.y) = 0;
	*reinterpret_cast<sf_u32 *>(this->color) = 0;
	*reinterpret_cast<sf_u32 *>(this->color2) = 0;
}

/*
========================
CVertex2D::getNormal
========================
*/
SMF_INLINE const CVec2D CVertex2D::getNormal() const {
	CVec2D n(normal[0],normal[1]);
	n.toNormal();	// after the normal has been compressed & uncompressed, it may not be normalized anymore
	return n;
}

/*
========================
CVertex2D::getNormalRaw
========================
*/
SMF_INLINE const CVec2D CVertex2D::getNormalRaw() const {
	CVec2D n(normal[0],normal[1]);
	// don't re-normalize just like we do in the vertex programs
	return n;
}






/*
========================
CVertex2D::lerp
========================
*/
SMF_INLINE void CVertex2D::lerp( const CVertex2D &a, const CVertex2D &b, const float f ) {
	xy = a.xy + f * ( b.xy - a.xy );
	//setTexCoord( MATH::lerp( a.getTexCoord(), b.getTexCoord(), f ) );
}

/*
========================
CVertex2D::lerpAll
========================
*/
SMF_INLINE void CVertex2D::lerpAll( const CVertex2D &a, const CVertex2D &b, const float f ) {
	xy = MATH::lerp( a.xy, b.xy, f );
	
	color[0] = (sf_u8)( a.color[0] + f * ( b.color[0] - a.color[0] ) );
	color[1] = (sf_u8)( a.color[1] + f * ( b.color[1] - a.color[1] ) );
	color[2] = (sf_u8)( a.color[2] + f * ( b.color[2] - a.color[2] ) );
	color[3] = (sf_u8)( a.color[3] + f * ( b.color[3] - a.color[3] ) );

	color2[0] = (sf_u8)( a.color2[0] + f * ( b.color2[0] - a.color2[0] ) );
	color2[1] = (sf_u8)( a.color2[1] + f * ( b.color2[1] - a.color2[1] ) );
	color2[2] = (sf_u8)( a.color2[2] + f * ( b.color2[2] - a.color2[2] ) );
	color2[3] = (sf_u8)( a.color2[3] + f * ( b.color2[3] - a.color2[3] ) );
}

/*
========================
CVertex2D::setNativeOrderColor
========================
*/
SMF_INLINE void CVertex2D::setNativeOrderColor( sf_u32 color ) {
	*reinterpret_cast<sf_u32 *>(this->color) = color;
}

/*
========================
CVertex2D::setColor
========================
*/
SMF_INLINE void CVertex2D::setColor( sf_u32 color ) {
	*reinterpret_cast<sf_u32 *>(this->color) = color;
}

/*
========================
CVertex2D::setColor
========================
*/
SMF_INLINE sf_u32 CVertex2D::getColor() const {
	return *reinterpret_cast<const sf_u32 *>(this->color);
}


/*
========================
CVertex2D::setColor
========================
*/
SMF_INLINE void CVertex2D::setColor2( sf_u32 color2 ) {
	*reinterpret_cast<sf_u32 *>(this->color2) = color2;
}

/*
========================
CVertex2D::clearColor2
========================
*/
SMF_INLINE void CVertex2D::clearColor2() {
	*reinterpret_cast<sf_u32 *>(this->color2) = 0x80808080;
}

/*
========================
CVertex2D::getColor2
========================
*/
SMF_INLINE sf_u32 CVertex2D::getColor2() const {
	return *reinterpret_cast<const sf_u32 *>(this->color2);
}


} //end GEO
} //end SMF
#endif /* _SMF__DRAWVERT_2D_2D_H__ */
