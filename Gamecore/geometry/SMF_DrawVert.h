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

#ifndef _SMF__DRAWVERT_H__
#define _SMF__DRAWVERT_H__
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

// The hardware converts a sf_u8 to a float by division with 255 and in the
// vertex programs we convert the floating-point value in the range [0, 1]
// to the range [-1, 1] by multiplying with 2 and subtracting 1.
#define VERTEX_BYTE_TO_FLOAT( x )		( (x) * ( 2.0f / 255.0f ) - 1.0f )
#define VERTEX_FLOAT_TO_BYTE( x )		CMath::ftob( ( (x) + 1.0f ) * ( 255.0f / 2.0f ) + 0.5f )

// The hardware converts a sf_u8 to a float by division with 255 and in the
// fragment programs we convert the floating-point value in the range [0, 1]
// to the range [-1, 1] by multiplying with 2 and subtracting 1.
// This is the conventional OpenGL mapping which specifies an exact
// representation for -1 and +1 but not 0. The DirectX 10 mapping is
// in the comments which specifies a non-linear mapping with an exact
// representation of -1, 0 and +1 but -1 is represented twice.
#define NORMALMAP_BYTE_TO_FLOAT( x )	VERTEX_BYTE_TO_FLOAT( x )	//( (x) - 128.0f ) * ( 1.0f / 127.0f )
#define NORMALMAP_FLOAT_TO_BYTE( x )	VERTEX_FLOAT_TO_BYTE( x )	//idMath::ftob( 128.0f + 127.0f * (x) + 0.5f )


/**
 * \class CVertex
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
class SMF_API CVertex {
public:
	/**
	\brief localização do vertex no espaço
	**/
	MATH::CVec3D			xyz;  //12 bytes
	/**
	\brief coordenadas de textura
	**/
	MATH::CVec2D			st;	   //8bytes. texture coord
	/**
	\brief A normal do vertex
	**/
	MATH::CVec3D			normal; //12 bytes
	/**
	\brief A tangents[0] representa o plano tangente ao vertex
	      \n tangents[1] representa a bitangente
	**/
	MATH::CVec3D			tangents[2]; //24bytes
	sf_u8					color[4];     // 4bytes
	sf_u8					color2[4];		// 4 bytes -- weights for skinning


	float				operator[]( const int index ) const;
	float &				operator[]( const int index );

	void				clear();

	const CVec3D		getNormal() const;
	const CVec3D		getNormalRaw() const;		// not re-normalized for renderbump

	// must be normalized already!
	void				setNormal( float x, float y, float z );
	void				setNormal( const CVec3D & n );

	const CVec3D		getTangent() const;
	const CVec3D		getTangentRaw() const;		// not re-normalized for renderbump

	// must be normalized already!
	void				setTangent( float x, float y, float z );
	void				setTangent( const CVec3D & t );

	 const CVec3D 		getBiTangent() const;
	 const CVec3D 		getBiTangentRaw() const;	// not re-normalized for renderbump

	 void				    setBiTangent( float x, float y, float z );
	 SMF_INLINE void		setBiTangent( const CVec3D & t );

	//s float				getBiTangentSign() const;
	//s sf_u8				getBiTangentSignBit() const;

	void				setTexCoordNative( const halfFloat_t s, const halfFloat_t t );
	void				setTexCoordNative( const float s, const float t );
	void				setTexCoord( const CVec2D & st );
	void				setTexCoord( float s, float t );
	void				setTexCoordS( float s );
	void				setTexCoordT( float t );
	const CVec2D		getTexCoord() const;
	const halfFloat_t	getTexCoordNativeS_half() const;
	const halfFloat_t	getTexCoordNativeT_half() const;
	const float			getTexCoordNativeS_float() const;
	const float			getTexCoordNativeT_float() const;

	// either 1.0f or -1.0f
	SMF_INLINE void		setBiTangentSign( float sign );
	SMF_INLINE void		setBiTangentSignBit( sf_u8 bit );

	void				lerp( const CVertex &a, const CVertex &b, const float f );
	void				lerpAll( const CVertex &a, const CVertex &b, const float f );

	void				setColor( sf_u32 color );
	void				setNativeOrderColor( sf_u32 color );
	sf_u32				getColor() const;

	void				setColor2( sf_u32 color );
	void				setNativeOrderColor2( sf_u32 color );
	void				clearColor2();
	sf_u32				getColor2() const;

	static CVertex	getSkinnedDrawVert( const CVertex & vert, const CMatJoint3x4 * joints );
	static CVec3D		getSkinnedDrawVertPosition( const CVertex & vert, const CMatJoint3x4 * joints );
};

#define JOINTQUAT_SIZE				(7*4)
#define JOINTMAT_SIZE				(4*3*4)
#define JOINTWEIGHT_SIZE			(4*4)


#define DRAWVERT_SIZE				64
#define DRAWVERT_SIZE_STR			"64"

#define DRAWVERT_XYZ_OFFSET			(0*4)
#define DRAWVERT_XYZ_OFFSET_STR			"(0*4)"
#define DRAWVERT_ST_OFFSET			(3*4)
#define DRAWVERT_NORMAL_OFFSET		(5*4)
#define DRAWVERT_NORMAL_OFFSET_STR		"(5*4)"

#define DRAWVERT_TANGENT0_OFFSET	(8*4)
#define DRAWVERT_TANGENT0_OFFSET_STR	"(8*4)"

#define DRAWVERT_TANGENT1_OFFSET	(11*4)
#define DRAWVERT_TANGENT1_OFFSET_STR	"(11*4)"

#define DRAWVERT_COLOR_OFFSET		(14*4)

assert_offsetof( CVertex, xyz,		DRAWVERT_XYZ_OFFSET );
assert_offsetof( CVertex, normal,	DRAWVERT_NORMAL_OFFSET );
//SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
//SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );
//assert_offsetof( CVertex, tangent,	DRAWVERT_TANGENT_OFFSET );

/*
========================
VertexFloatToByte

Assumes input is in the range [-1, 1]
========================
*/
SMF_INLINE void VertexFloatToByte( const float & x, const float & y, const float & z, sf_u8 * bval ) {
	assert_4_byte_aligned( bval );	// for __stvebx

#ifdef ID_WIN_X86_SSE2_INTRIN

	const sf_m128 vector_float_one			= { 1.0f, 1.0f, 1.0f, 1.0f };
	const sf_m128 vector_float_half			= { 0.5f, 0.5f, 0.5f, 0.5f };
	const sf_m128 vector_float_255_over_2	= { 255.0f / 2.0f, 255.0f / 2.0f, 255.0f / 2.0f, 255.0f / 2.0f };

	const sf_m128 xyz = _mm_unpacklo_ps( _mm_unpacklo_ps( _mm_load_ss( &x ), _mm_load_ss( &z ) ), _mm_load_ss( &y ) );
	const sf_m128 xyzScaled = _mm_madd_ps( _mm_add_ps( xyz, vector_float_one ), vector_float_255_over_2, vector_float_half );
	const __m128i xyzInt = _mm_cvtps_epi32( xyzScaled );
	const __m128i xyzShort = _mm_packs_epi32( xyzInt, xyzInt );
	const __m128i xyzChar = _mm_packus_epi16( xyzShort, xyzShort );
	const __m128i xyz16 = _mm_unpacklo_epi8( xyzChar, _mm_setzero_si128() );

	bval[0] = (sf_u8)_mm_extract_epi16( xyz16, 0 );	// cannot use _mm_extract_epi8 because it is an SSE4 instruction
	bval[1] = (sf_u8)_mm_extract_epi16( xyz16, 1 );
	bval[2] = (sf_u8)_mm_extract_epi16( xyz16, 2 );

#else

	bval[0] = VERTEX_FLOAT_TO_BYTE( x );
	bval[1] = VERTEX_FLOAT_TO_BYTE( y );
	bval[2] = VERTEX_FLOAT_TO_BYTE( z );

#endif
}

/*
========================
CVertex::operator[]
========================
*/
SMF_INLINE float CVertex::operator[]( const int index ) const {
	SMF_ASSERT( index >= 0 && index < 5 );
	return ((float *)(&xyz))[index];
}

/*
========================
CVertex::operator[]
========================
*/
SMF_INLINE float	&CVertex::operator[]( const int index ) {
	SMF_ASSERT( index >= 0 && index < 5 );
	return ((float *)(&xyz))[index];
}

/*
========================
CVertex::clear
========================
*/
SMF_INLINE void CVertex::clear() {
	*reinterpret_cast<sf_u32 *>(&this->xyz.x) = 0;
	*reinterpret_cast<sf_u32 *>(&this->xyz.y) = 0;
	*reinterpret_cast<sf_u32 *>(&this->xyz.z) = 0;
	*reinterpret_cast<sf_u32 *>(&this->st.x) = 0;
	*reinterpret_cast<sf_u32 *>(&this->st.y) = 0;

	//*reinterpret_cast<sf_u32 *>(this->normal) = 0x00FF8080;	// x=0, y=0, z=1
	normal=CVec3D::unitZ;
	//*reinterpret_cast<sf_u32 *>(this->tangent) = 0xFF8080FF;	// x=1, y=0, z=0
	tangents[0]=CVec3D::unitX;
	tangents[1]=CVec3D::unitX;
	*reinterpret_cast<sf_u32 *>(this->color) = 0;
	*reinterpret_cast<sf_u32 *>(this->color2) = 0;
}

/*
========================
CVertex::getNormal
========================
*/
SMF_INLINE const CVec3D CVertex::getNormal() const {
//	CVec3D n(	VERTEX_BYTE_TO_FLOAT( normal[0] ),
//				VERTEX_BYTE_TO_FLOAT( normal[1] ),
//				VERTEX_BYTE_TO_FLOAT( normal[2] ) );
	CVec3D n(normal[0],normal[1],normal[2]);
	n.toNormal();	// after the normal has been compressed & uncompressed, it may not be normalized anymore
	return n;
}

/*
========================
CVertex::getNormalRaw
========================
*/
SMF_INLINE const CVec3D CVertex::getNormalRaw() const {
//	CVec3D n(	VERTEX_BYTE_TO_FLOAT( normal[0] ),
//				VERTEX_BYTE_TO_FLOAT( normal[1] ),
//				VERTEX_BYTE_TO_FLOAT( normal[2] ) );
	CVec3D n(normal[0],normal[1],normal[2]);
	// don't re-normalize just like we do in the vertex programs
	return n;
}

/*
========================
CVertex::setNormal
must be normalized already!
========================
*/
SMF_INLINE void CVertex::setNormal( const CVec3D & n ) {
	//VertexFloatToByte( n.x, n.y, n.z, normal );
	normal=n;
}

/*
========================
CVertex::setNormal
========================
*/
SMF_INLINE void CVertex::setNormal( float x, float y, float z ) {
	//VertexFloatToByte( x, y, z, normal );
	normal.x=x;
	normal.x=y;
	normal.x=z;
}

/*
========================
&CVertex::getTangent
========================
*/
SMF_INLINE const CVec3D CVertex::getTangent() const {
//	CVec3D t(	VERTEX_BYTE_TO_FLOAT( tangent[0] ),
//				VERTEX_BYTE_TO_FLOAT( tangent[1] ),
//				VERTEX_BYTE_TO_FLOAT( tangent[2] ) );
	CVec3D t(tangents[0].x,tangents[0].y,tangents[0].z);
	t.toNormal();
	return t;
}

/*
========================
&CVertex::getTangentRaw
========================
*/
SMF_INLINE const CVec3D CVertex::getTangentRaw() const {
//	CVec3D t(	VERTEX_BYTE_TO_FLOAT( tangent[0] ),
//				VERTEX_BYTE_TO_FLOAT( tangent[1] ),
//				VERTEX_BYTE_TO_FLOAT( tangent[2] ) );
	CVec3D t(tangents[0].x,tangents[0].y,tangents[0].z);
	// don't re-normalize just like we do in the vertex programs
	return t;
}

/*
========================
CVertex::setTangent
========================
*/
SMF_INLINE void CVertex::setTangent( float x, float y, float z ) {
	//VertexFloatToByte( x, y, z, tangent );
	tangents[0].x=x;
	tangents[0].y=y;
	tangents[0].z=z;
}

/*
========================
CVertex::setTangent
========================
*/
SMF_INLINE void CVertex::setTangent( const CVec3D & t ) {
	//VertexFloatToByte( t.x, t.y, t.z, tangent );
	tangents[0]=t;
}


/*
========================
CVertex::getBiTangent
========================
*/
SMF_INLINE const CVec3D CVertex::getBiTangent() const {
#if 0	// derive from the normal, tangent, and bitangent direction flag
	CVec3D bitangent;
	bitangent.cross( getNormal(), getTangent() );
	bitangent *= getBiTangentSign();
	return bitangent;
#endif
		CVec3D t(tangents[1].x,tangents[1].y,tangents[1].z);
	t.toNormal();
	return t;
}

/*
========================
CVertex::getBiTangentRaw
========================
*/
SMF_INLINE const CVec3D CVertex::getBiTangentRaw() const {
#if 0
	// derive from the normal, tangent, and bitangent direction flag
	// don't re-normalize just like we do in the vertex programs
	CVec3D bitangent;
	bitangent.cross( getNormalRaw(), getTangentRaw() );
	bitangent *= getBiTangentSign();
	return bitangent;
#endif
		CVec3D t(tangents[1].x,tangents[1].y,tangents[1].z);
	// don't re-normalize just like we do in the vertex programs
	return t;

}

/*
========================
CVertex::setBiTangent
========================
*/
SMF_INLINE void CVertex::setBiTangent( float x, float y, float z ) {
	setBiTangent( CVec3D( x, y, z ) );
}

/*
========================
CVertex::setBiTangent
The way I calculated them is this if the normal is [x,y,z]:

Tangent: [-y,x,z] x [x,y,z] (thats cross product)

Bitangent: Tangent x Normal (again cross product)
========================
*/
SMF_INLINE void CVertex::setBiTangent( const CVec3D &t ) {
	CVec3D bitangent;
	bitangent.cross( getNormal(), getTangent() );
	//setBiTangentSign( bitangent * t );
	tangents[1]=bitangent;
}
#if 0
/*
========================
CVertex::getBiTangentSign
========================
*/

SMF_INLINE float CVertex::getBiTangentSign() const {
	return ( tangent[3] < 128 ) ? -1.0f : 1.0f;
}

/*
========================
CVertex::getBiTangentSignBit
========================
*/
SMF_INLINE sf_u8 CVertex::getBiTangentSignBit() const {
	return ( tangent[3] < 128 ) ? 1 : 0;
}

/*
========================
CVertex::setBiTangentSign
========================
*/
SMF_INLINE void CVertex::setBiTangentSign( float sign ) {
	tangent[3] = ( sign < 0.0f ) ? 0 : 255;
}

/*
========================
CVertex::setBiTangentSignBit
========================
*/
SMF_INLINE void CVertex::setBiTangentSignBit( sf_u8 sign ) {
	tangent[3] = sign ? 0 : 255;
}
#endif
/*
========================
CVertex::lerp
========================
*/
SMF_INLINE void CVertex::lerp( const CVertex &a, const CVertex &b, const float f ) {
	xyz = a.xyz + f * ( b.xyz - a.xyz );
	setTexCoord( MATH::lerp( a.getTexCoord(), b.getTexCoord(), f ) );
}

/*
========================
CVertex::lerpAll
========================
*/
SMF_INLINE void CVertex::lerpAll( const CVertex &a, const CVertex &b, const float f ) {
	xyz = MATH::lerp( a.xyz, b.xyz, f );
	setTexCoord( MATH::lerp( a.getTexCoord(), b.getTexCoord(), f ) );

	CVec3D normal = MATH::lerp( a.getNormal(), b.getNormal(), f );
	CVec3D tangent = MATH::lerp( a.getTangent(), b.getTangent(), f );
	//CVec3D bitangent = MATH::lerp( a.getBiTangent(), b.getBiTangent(), f );
	normal.toNormal();
	tangents[0].toNormal();
	tangents[1].toNormal();
	//bitangent.toNormal();
	setNormal( normal );
	setTangent( tangent );
	//setBiTangent( bitangent );

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
CVertex::setNativeOrderColor
========================
*/
SMF_INLINE void CVertex::setNativeOrderColor( sf_u32 color ) {
	*reinterpret_cast<sf_u32 *>(this->color) = color;
}

/*
========================
CVertex::setColor
========================
*/
SMF_INLINE void CVertex::setColor( sf_u32 color ) {
	*reinterpret_cast<sf_u32 *>(this->color) = color;
}

/*
========================
CVertex::setColor
========================
*/
SMF_INLINE sf_u32 CVertex::getColor() const {
	return *reinterpret_cast<const sf_u32 *>(this->color);
}

/*
========================
CVertex::setTexCoordNative
========================
*/
SMF_INLINE void CVertex::setTexCoordNative( const halfFloat_t s, const halfFloat_t t ) {
	st[0] = s;
	st[1] = t;
}

/*
========================
CVertex::setTexCoord
========================
*/
SMF_INLINE void CVertex::setTexCoord( const CVec2D & st ) {
	setTexCoordS( st.x );
	setTexCoordT( st.y );
}

/*
========================
CVertex::setTexCoord
========================
*/
SMF_INLINE void CVertex::setTexCoord( float s, float t ) {
	setTexCoordS( s );
	setTexCoordT( t );
}

/*
========================
CVertex::setTexCoordS
========================
*/
SMF_INLINE void CVertex::setTexCoordS( float s ) {
	//st[0] = F32toF16( s );
	st[0]=s;
}

/*
========================
CVertex::setTexCoordT
========================
*/
SMF_INLINE void CVertex::setTexCoordT( float t ) {
	//st[1] = F32toF16( t );
	st[1]=t;
}

/*
========================
CVertex::getTexCoord
========================
*/
SMF_INLINE const CVec2D	CVertex::getTexCoord() const {
	//return CVec2D( F16toF32( st[0] ), F16toF32( st[1] ) );
	return st;
}

/*
========================
CVertex::getTexCoordNativeS
========================
*/
SMF_INLINE const halfFloat_t CVertex::getTexCoordNativeS_half() const {
	return F32toF16(st[0]);
}

/*
========================
CVertex::getTexCoordNativeT
========================
*/
SMF_INLINE const halfFloat_t CVertex::getTexCoordNativeT_half() const {
	return F32toF16(st[1]);
}
SMF_INLINE const float CVertex::getTexCoordNativeS_float() const {
	return st[0];
}


SMF_INLINE const float CVertex::getTexCoordNativeT_float() const {
	return st[1];
}

/*
========================
CVertex::setNativeOrderColor2
========================
*/
SMF_INLINE void CVertex::setNativeOrderColor2( sf_u32 color2 ) {
	*reinterpret_cast<sf_u32 *>(this->color2) = color2;
}

/*
========================
CVertex::setColor
========================
*/
SMF_INLINE void CVertex::setColor2( sf_u32 color2 ) {
	*reinterpret_cast<sf_u32 *>(this->color2) = color2;
}

/*
========================
CVertex::clearColor2
========================
*/
SMF_INLINE void CVertex::clearColor2() {
	*reinterpret_cast<sf_u32 *>(this->color2) = 0x80808080;
}

/*
========================
CVertex::getColor2
========================
*/
SMF_INLINE sf_u32 CVertex::getColor2() const {
	return *reinterpret_cast<const sf_u32 *>(this->color2);
}

/*
========================
WriteDrawVerts16

Use 16-sf_u8 in-order SIMD writes because the destVerts may live in write-combined memory
========================
*/
SMF_INLINE void WriteDrawVerts16( CVertex * destVerts, const CVertex * localVerts, int numVerts ) {
	//assert_sizeof( CVertex, 32 );
	assert_16_byte_aligned( destVerts );
	assert_16_byte_aligned( localVerts );

#ifdef ID_WIN_X86_SSE2_INTRIN

	for ( int i = 0; i < numVerts; i++ ) {
		__m128i v0 = _mm_load_si128( (const __m128i *)( (sf_u8 *)( localVerts + i ) +  0 ) );
		__m128i v1 = _mm_load_si128( (const __m128i *)( (sf_u8 *)( localVerts + i ) + 16 ) );
		_mm_stream_si128( (__m128i *)( (sf_u8 *)( destVerts + i ) +  0 ), v0 );
		_mm_stream_si128( (__m128i *)( (sf_u8 *)( destVerts + i ) + 16 ), v1 );
	}

#else

	memcpy( destVerts, localVerts, numVerts * sizeof( CVertex ) );

#endif
}

/*
=====================
CVertex::getSkinnedDrawVert
=====================
*/
SMF_INLINE CVertex CVertex::getSkinnedDrawVert( const CVertex & vert, const CMatJoint3x4 * joints ) {
	if ( joints == NULL ) {
		return vert;
	}
#if 0
	const CMatJoint3x4 & j0 = joints[vert.color[0]];
	const CMatJoint3x4 & j1 = joints[vert.color[1]];
	const CMatJoint3x4 & j2 = joints[vert.color[2]];
	const CMatJoint3x4 & j3 = joints[vert.color[3]];

	const float w0 = vert.color2[0] * ( 1.0f / 255.0f );
	const float w1 = vert.color2[1] * ( 1.0f / 255.0f );
	const float w2 = vert.color2[2] * ( 1.0f / 255.0f );
	const float w3 = vert.color2[3] * ( 1.0f / 255.0f );

	CMatJoint3x4 accum;
	CMatJoint3x4::mul( accum, j0, w0 );
	CMatJoint3x4::Mad( accum, j1, w1 );
	CMatJoint3x4::Mad( accum, j2, w2 );
	CMatJoint3x4::Mad( accum, j3, w3 );

	CVertex outVert;
	outVert.xyz = accum * CVec4D( vert.xyz.x, vert.xyz.y, vert.xyz.z, 1.0f );
	outVert.setTexCoordNative( vert.getTexCoordNativeS(), vert.getTexCoordNativeT() );
	outVert.setNormal( accum * vert.getNormal() );
	outVert.setTangent( accum * vert.getTangent() );
	outVert.tangent[3] = vert.tangent[3];
	for ( int i = 0; i < 4; i++ ) {
		outVert.color[i] = vert.color[i];
		outVert.color2[i] = vert.color2[i];
	}
	return outVert;
#endif
	CVertex newvertex;
	return newvertex;
}

/*
=====================
CVertex::getSkinnedDrawVertPosition
=====================
*/
SMF_INLINE CVec3D CVertex::getSkinnedDrawVertPosition( const CVertex & vert, const CMatJoint3x4 * joints ) {
	if ( joints == NULL ) {
		return vert.xyz;
	}
#if 0
	const CMatJoint3x4 & j0 = joints[vert.color[0]];
	const CMatJoint3x4 & j1 = joints[vert.color[1]];
	const CMatJoint3x4 & j2 = joints[vert.color[2]];
	const CMatJoint3x4 & j3 = joints[vert.color[3]];

	const float w0 = vert.color2[0] * ( 1.0f / 255.0f );
	const float w1 = vert.color2[1] * ( 1.0f / 255.0f );
	const float w2 = vert.color2[2] * ( 1.0f / 255.0f );
	const float w3 = vert.color2[3] * ( 1.0f / 255.0f );
#endif
	CMatJoint3x4 accum;
	//s CMatJoint3x4::mul( accum, j0, w0 );
	//s CMatJoint3x4::Mad( accum, j1, w1 );
	//s CMatJoint3x4::Mad( accum, j2, w2 );
	//s CMatJoint3x4::Mad( accum, j3, w3 );

	return accum * CVec4D( vert.xyz.x, vert.xyz.y, vert.xyz.z, 1.0f );
}

/*
===============================================================================
Shadow vertex
===============================================================================
*/
class CShadowVert {
public:
	CVec4D			xyzw;

	void			clear();
	static int		 createShadowCache( CShadowVert * vertexCache, const CVertex *verts, const int numVerts );
};

#define SHADOWVERT_XYZW_OFFSET		(0)

assert_offsetof( CShadowVert, xyzw, SHADOWVERT_XYZW_OFFSET );

SMF_INLINE void CShadowVert::clear() {
	xyzw.toZero();
}

/*
===============================================================================
Skinned Shadow vertex
===============================================================================
*/
class CShadowVertSkinned {
public:
	CVec4D			xyzw;
	sf_u8			color[4];
	sf_u8			color2[4];
	sf_u8			pad[8];		// pad to multiple of 32-sf_u8 for glDrawElementsBaseVertex

	void			clear();
	static int		 createShadowCache( CShadowVertSkinned * vertexCache, const CVertex *verts, const int numVerts );
};

#define SHADOWVERTSKINNED_XYZW_OFFSET		(0)
#define SHADOWVERTSKINNED_COLOR_OFFSET		(16)
#define SHADOWVERTSKINNED_COLOR2_OFFSET		(20)

assert_offsetof( CShadowVertSkinned, xyzw, SHADOWVERTSKINNED_XYZW_OFFSET );
assert_offsetof( CShadowVertSkinned, color, SHADOWVERTSKINNED_COLOR_OFFSET );
assert_offsetof( CShadowVertSkinned, color2, SHADOWVERTSKINNED_COLOR2_OFFSET );

SMF_INLINE void CShadowVertSkinned::clear() {
	xyzw.toZero();
}

} //end GEO
} //end SMF
#endif /* !__DRAWVERT_H__ */
