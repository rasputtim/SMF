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

#ifndef _SMF__MATH_SIMD_GENERIC_H__
#define _SMF__MATH_SIMD_GENERIC_H__
#include "../SMF_Config.h"
#include "SMF_Simd.h"

namespace SMF {
namespace GEO{
class CAABBox;
}
using namespace GEO;
namespace MATH{

/*
===============================================================================

	Generic implementation of CSIMDProcessor

===============================================================================
*/
/**
 * \class CSIMD_Generic
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação genérica de um processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 *
 **/
class SMF_API CSIMD_Generic : public CSIMDProcessor {
public:
	CSIMD_Generic():CSIMDProcessor(CPUID_GENERIC){};
	CSIMD_Generic(cpuid_t id):CSIMDProcessor(id){};
	~CSIMD_Generic(){};
	virtual const char * VPCALL getName() const;
	/**
	\brief   dst[i] = constant + src[i];
	**/
	virtual void VPCALL add( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL add( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL sub( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL sub( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mul( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mul( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL div( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL div( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mulAdd( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mulAdd( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mulSub( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mulSub( float *dst,			const float *src0,		const float *src1,		const int count );

	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CVec3D *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CPlane *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CVertex *src,	const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVec3D *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CPlane *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVertex *src,	const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D *src0,		const CVec3D *src1,		const int count );
	virtual void VPCALL dot( float &dot,			const float *src1,		const float *src2,		const int count );

	virtual void VPCALL cmpGT( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGE( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLT( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLE( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );

	virtual void VPCALL minMax( float &min,			float &max,				const float *src,		const int count );
	virtual	void VPCALL minMax( CVec2D &min,		CVec2D &max,			const CVec2D *src,		const int count );
	virtual void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVec3D *src,		const int count );
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int count );
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int *indexes,		const int count );

	virtual void VPCALL clamp( float *dst,			const float *src,		const float min,		const float max,		const int count );
	virtual void VPCALL clampMin( float *dst,		const float *src,		const float min,		const int count );
	virtual void VPCALL clampMax( float *dst,		const float *src,		const float max,		const int count );

	virtual void VPCALL memCopy( void *dst,			const void *src,		const int count );
	virtual void VPCALL memSet( void *dst,			const int val,			const int count );

	virtual void VPCALL zero16( float *dst,			const int count );
	virtual void VPCALL negate16( float *dst,		const int count );
	virtual void VPCALL copy16( float *dst,			const float *src,		const int count );
	virtual void VPCALL add16( float *dst,			const float *src1,		const float *src2,		const int count );
	virtual void VPCALL sub16( float *dst,			const float *src1,		const float *src2,		const int count );
	virtual void VPCALL mul16( float *dst,			const float *src1,		const float constant,	const int count );
	virtual void VPCALL addAssign16( float *dst,	const float *src,		const int count );
	virtual void VPCALL subAssign16( float *dst,	const float *src,		const int count );
	virtual void VPCALL mulAssign16( float *dst,	const float constant,	const int count );

	virtual bool VPCALL matX_inverse_4x4(float* mat);
	virtual void VPCALL matX_MultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 );
	virtual void VPCALL matX_TransposeMultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 );
	virtual void VPCALL matX_LowerTriangularsolve( const CMatXD &L, float *x, const float *b, const int n, int skip = 0 );
	virtual void VPCALL matX_LowerTriangularsolveTranspose( const CMatXD &L, float *x, const float *b, const int n );
	virtual bool VPCALL matX_LDLTFactor( CMatXD &mat, CVecXD &invDiag, const int n );

	virtual void VPCALL blendJoints( CJointQuaternion *joints, const CJointQuaternion *blendJoints, const float lerp, const int *index, const int numJoints );
	virtual void VPCALL convertJointQuatsToJointMats( CMatJoint3x4 *jointMats, const CJointQuaternion *jointQuats, const int numJoints );
	virtual void VPCALL convertJointMatsToJointQuats( CJointQuaternion *jointQuats, const CMatJoint3x4 *jointMats, const int numJoints );
	virtual void VPCALL transformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint );
	virtual void VPCALL untransformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint );
	virtual void VPCALL transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights );
	virtual void VPCALL tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL decalPointCull( sf_u8 *cullBits, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL overlayPointCull( sf_u8 *cullBits, CVec2D *texCoords, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL deriveTriPlanes( CPlane *planes, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual void VPCALL  deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	//virtual void VPCALL deriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts );
	virtual void VPCALL  normalizeTangents( CVertex *verts, const int numVerts );
	virtual void VPCALL  createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual void VPCALL  createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual int  VPCALL  createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts );
	virtual int  VPCALL  createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts );

	virtual void VPCALL  upSamplePCMTo44kHz( float *dest, const short *pcm, const int numSamples, const int kHz, const int numChannels );
	virtual void VPCALL  upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels );
	virtual void VPCALL  mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] );
	virtual void VPCALL  mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] );
	virtual void VPCALL  mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] );
	virtual void VPCALL  mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] );
	virtual void VPCALL  mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples );

	//=================================
	//métodos Vector3D

	virtual void VPCALL vector3D_Sum(CVec3D* pOut, const CVec3D* pIn);
	virtual void VPCALL vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2);
	virtual void VPCALL vector3D_Diff(CVec3D* pLeft, CVec3D* pRight);
	virtual void VPCALL vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_Scale(CVec3D* pOut, float scalar);
	virtual void VPCALL vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar);
	virtual float VPCALL vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2);
	virtual float VPCALL vector3D_LengthSq(const CVec3D* pVec);
	virtual float VPCALL vector3D_Length(const CVec3D* pVec);
	virtual void VPCALL vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_Normalize(CVec3D* pVec);
	virtual void VPCALL vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec);
	virtual float VPCALL vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2);

	virtual void VPCALL vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn);
	virtual void VPCALL vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2);
	virtual void VPCALL vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight);
	virtual void VPCALL vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_AlignedScale(CVec3D* pOut, float scalar);
	virtual void VPCALL vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar);
	virtual float VPCALL vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2);
	virtual float VPCALL vector3D_AlignedLengthSq(const CVec3D* pVec);
	virtual float VPCALL vector3D_AlignedLength(const CVec3D* pVec);
	virtual void VPCALL vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
	virtual void VPCALL vector3D_AlignedNormalize(CVec3D* pVec);
	virtual void VPCALL vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec);
	virtual float VPCALL vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2);

	virtual float VPCALL vector4D_Dot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2);
	virtual void VPCALL vector4D_Sum(CVec4D* pOut, const CVec4D* pIn);
	virtual void VPCALL vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2);
	virtual void VPCALL vector4D_Diff(CVec4D* pLeft, CVec4D* pRight);
	virtual void VPCALL vector4D_DiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight);
	virtual void VPCALL vector4D_Scale(CVec4D* pOut, float scalar);
	virtual void VPCALL vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar);
	virtual float VPCALL vector4D_LengthSq(const CVec4D* pVec);
	virtual float VPCALL vector4D_Length(const CVec4D* pVec);
	virtual void VPCALL vector4D_Normalize(CVec4D* pVec);
	virtual void VPCALL vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec);
	virtual float VPCALL vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2);

	virtual float VPCALL vector4D_AlignedDot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2);
	virtual void VPCALL vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn);
	virtual void VPCALL vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2);
	virtual void VPCALL vector4D_AlignedDiff(CVec4D* pLeft, CVec4D* pRight);
	virtual void VPCALL vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight);
	virtual void VPCALL vector4D_AlignedScale(CVec4D* pOut, float scalar);
	virtual void VPCALL vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar);
	virtual float VPCALL vector4D_AlignedLengthSq(const CVec4D* pVec);
	virtual float VPCALL vector4D_AlignedLength(const CVec4D* pVec);
	virtual void VPCALL vector4D_AlignedNormalize(CVec4D* pVec);
	virtual void VPCALL vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec);
	virtual float VPCALL vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2);
//===========trigonometry=====================================
	virtual float VPCALL  invSqrt( float x );
	virtual void  VPCALL  InvSqrt4( float x[4] );
	virtual float VPCALL  sinZeroHalfPI( float x );
	virtual void  VPCALL  sin4ZeroHalfPI( float entrada[4], float saida[4] ) ;
	virtual float VPCALL  sin( float a );
	virtual void  VPCALL  sin4( float entrada[4], float saida[4] );
	virtual float VPCALL  cosZeroHalfPI( float x );
	virtual void  VPCALL  cos4ZeroHalfPI( float entrada[4], float saida[4] );
	virtual float VPCALL  cos( float a );
	virtual void  VPCALL  cos4( float entrada[4], float saida[4] );
	virtual void  VPCALL  sincos( float a, float &sin, float &cos );
	virtual void  VPCALL  sincos4( float a[4], float sin[4], float cos[4] );
	virtual float VPCALL  aTanPositive( float y, float x );
	virtual void  VPCALL  aTan4Positive( float y[4], float x[4], float resultado[4] );
	virtual float VPCALL  atan( float y, float x );
	virtual void  VPCALL  aTan4( float y[4], float x[4], float resultado[4] );

//=============================MATRIZ================================================
//============CMat4D===================================================
/**
\brief soma pMat e pIn e guarda o resultado em pMat
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matriz que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Sum(CMat4D* pMat, const CMat4D* pIn);
/**
\brief soma pIn1 e pIn2 e armazena o resultado em pMat
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\param pMat Matriz que guardará o resultado da soma
\param pIn1 Matriz que será somada a pIn2
\param pIn2 Matriz que será somada a pIn1
**/
virtual void  VPCALL mat4D_SumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
**/
virtual void  VPCALL mat4D_Diff(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\brief Subtracts pRight from pLeft and stores the result in pMat.
**/
virtual void  VPCALL mat4D_DiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Scale(CMat4D* pMat, float scalar);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar);
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight);
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\brief Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\note Transposes the matrix pIn stores the result in pIn.
**/
virtual void  VPCALL mat4D_Transpose(CMat4D* pIn);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
**/
virtual void  VPCALL mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn).

**/
virtual void  VPCALL mat4D_VectorMultiply(CVec3D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat).

**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pIn)
\param pOut4D Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
**/
virtual void  VPCALL mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat);
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
**/
virtual void  VPCALL mat4D_ToRotate(CMat4D* pMat, float yaw, float pitch, float roll);
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação

**/
virtual void  VPCALL mat4D_ToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll);

/**
\brief soma pMat e pIn e guarda o resultado em pMat
\note pMat,pIn devem estar alinhados em 16bytes
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matria que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
\note pMat, pIn must be 16 bytes aligned
**/
virtual void  VPCALL mat4D_AlignedSum(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\note  pMat,pIn1,pIn2 devem estar alinhados em 16bytes.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn1,pIn2 must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedSumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn devem estar alinhados em 16bytes.
\note pMat,pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiff(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\note Subtracts pRight from pLeft and stores the result in pMat.
\note pMat,pLeft,pRight devem estar alinhados em 16bytes.
\note pMat,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScale(CMat4D* pMat, float scalar);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar);
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pLeft,pRight devem estar alinhados em 16bytes.
\note pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight);
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pOut,pLeft,pRight devem estar alinhados em 16bytes.
\note pOut,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\brief Transposes the matrix pIn stores the result in pIn.
\note pIn deve estar alinhada em 16bytes.
\note pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTranspose(CMat4D* pIn);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut4D * pIn)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut Vetor multiplicador e que armazenará o resultado
\param pIn Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with SIMDx86Matrix_Vector4Multiply().
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec3D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with SIMDx86Matrix_Vector4MultiplyOf().
\note pOut, pIn, pMat devem estar alinhados em 16bytes.
\note pOut, pIn, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pMat)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut4D Vetor multiplicador e que armazenará o resultado
\param pMat Matriz multiplicadora
\note pOut4D, pMat devem estar alinhados em 16bytes.
\note pOut4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
\note pOut4D, pIn4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat);
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedToRotate(CMat4D* pMat, float yaw, float pitch, float roll);
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll);

virtual void  VPCALL quat_to_mat4x4(sf_m128 q, sf_m128 t, sf_m128 *m);

/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE 4.1, we can use the dpps (dot product) instruction, _mm_dp_ps intrinsic.
/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE3, we can repeatedly use haddps to accumulate the result.
virtual sf_m128  VPCALL mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);
virtual  void   VPCALL mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);

virtual sf_m128  VPCALL colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);


/// Compute the product M*v, where M is a 3x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
virtual sf_m128  VPCALL mat3x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);

virtual CVec3D  VPCALL mat3x4_mul_vec(const sf_m128 *matrix, sf_m128 vector);

/**
\brief multiplica duas matrizes CMat4D  m1*m2
\param [out] out resultado da multiplicação
\param m1 matriz a ser mutiplicada
\param m2 matria a ser multiplicada
**/
virtual void  VPCALL mat4x4_mul_dpps(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual void  VPCALL mat4x4_mul_dpps_2(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual void  VPCALL mat4x4_mul_dpps_3(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);


virtual void  VPCALL mat3x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual float  VPCALL mat4x4_inverse(const CMat4D *mat, CMat4D *out);

/// Inverts a 3x4 affine transformation matrix (in row-major format) that only consists of rotation (+possibly mirroring) and translation.
virtual void  VPCALL mat3x4_inverse_orthonormal(sf_m128 *mat, sf_m128 *out);
virtual sf_m128  VPCALL  newtonRhapsonRecipStep(sf_m128 recip, sf_m128 estimate);

virtual sf_m128  VPCALL  newtonRhapsonRecip(sf_m128 recip);

/// Computes the determinant of a 4x4 matrix.
virtual float  VPCALL mat4x4_determinant(const CMat4D *row);

/// Computes the determinant of a 3x4 matrix stored in row-major format. (Treated as a square matrix with last row [0,0,0,1])
virtual float  VPCALL mat3x4_determinant(const sf_m128 *row);

virtual void  VPCALL mat3x4_transpose(const sf_m128 *src, sf_m128 *dst);





//=========================Quaternion========================================
/**
\brief Normaiza o Quaternion e armazena o resultado em pQuat
\param pQuat Quaternion a ser normalizado
\brief Normalizes pQuat and stores it in pQuat.
**/
virtual void  VPCALL  quaternion_Normalize(CQuaternion* pQuat);

/**
\brief Normaiza o Quaternion e armazena o resultado em pOut
\param pQuat Quaternion a ser normalizado
\param pOut Quaternion que recebera o resultado
\brief Normalizes pQuat and stores it in pOut.
**/

virtual void VPCALL  quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat);
/**
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\param pLeft Quaternion que será multiplicado e receberá o resultado
\param pRight Quaternion que será multiplicado
\brief Multiplies pLeft by pRight and stores the result in pLeft. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
**/
virtual void VPCALL  quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight);
/**
\if pt_br
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\elseif us_en
\brief Multiplies pLeft by pRight and stores the result in pOut. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
\endif
\param pLeft Quaternion que será multiplicado
\param pRight Quaternion que será multiplicado
\param pOut Quaternion que receberá o resultado
**/
virtual void VPCALL  quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight);

/// Converts a quaternion to a row-major matrix.
/// From http://renderfeather.googlecode.com/hg-history/034a1900d6e8b6c92440382658d2b01fc732c5de/Doc/optimized%20Matrix%20quaternion%20conversion.pdf
virtual void VPCALL quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m);
virtual sf_m128 VPCALL quat_transform_vec4(sf_m128 quat, sf_m128 vec);
virtual sf_m128 VPCALL quat_mul_quat(sf_m128 q1, sf_m128 q2);
virtual sf_m128 VPCALL quat_div_quat(sf_m128 q1, sf_m128 q2);

//================CPlane=======================================
/**
\brief constrói um plano através de três pontos. O Plano não é normalizado
\param pOut O Plano que será criado
\param pA primeiro ponto
\param pB segundo ponto
\param pC terceiro ponto
\brief This function constructs a plane from three points. The plane is not normalized.
**/
virtual void   VPCALL  plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC);
/**
\brief retorna a distância sem sinal(módulo), entre um ponto e um plano
\param pPlane Plano a seu utilizado no cálculo
\param pPoint ponto a ser utilizado no cálculo
\brief This function returns the unsigned distance between a point and a plane.
**/
virtual float  VPCALL  plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint);
/**
\brief calcula o produto escalar entre o plano e um ponto com w=1.0f. Util para classificar
       um ponto em relação ao plano: um valor menor que zero significa que o ponto está atraz do plano
	   e um valor maoir que zero significa que o ponto está na frente do plano
\param pPlane Plano que será utilizado no cálculo
\param pVec Ponto que será utilizado no cálculo
\brief This function returns the dot product between a plane and a point with w = 1.0.
This is useful for classifying a point in relation to a plane: a value of less than zero
implies the point it behind the plane, a value of zero implies the point is on the plane,
and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot(const CPlane* pPlane, const CVec3D* pVec);
/**
\brief retorna o produto escalar entre o plano e um ponto 4D
\param pPlane plano que será utilizado no cálculo do produto
\param pVec4 ponto de 4 dimensões que será utilizado no calculo do produto
\brief This function returns the dot product between a plane and a 4D point.
This is useful for classifying a point in relation to a plane: a value of less than zero implies the point it behind the plane, a value of zero implies the point is on the plane, and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4);
/**
\brief retorna o produoto escalar entre o plano e outra normal. representa o cosseno do angulo entre os dois, se ambos estiverem normalizados
\param pPlane Plano que será utilizado para calcular o produto
\param pVec normal do segundo plano
\return O resultado do produto
\brief This functions returns the dot product between the plane's normal and another normal.
This represents the cosine of the angle between the two, if both are normalized.

**/
virtual float  VPCALL  plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec);
/**
\brief calcula o produto escalar entre dois planos normalizados.
\param pPlane1 Plano que será multiplicado
\param pPlane2 Plano que será multiplicado
\return o resultado do produto
\brief This functions returns the dot product between the two planes' normals.
This represents the cosine of the dihedral angle between the two, if both are normalized.
**/
virtual float  VPCALL  plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2);
/**
\brief Normaliza o plano pOut e armazena o resultado nele mesmo
\param pOut plano que será normalizado e armazenará o resultado
\brief This function normalizes the plane such that the magnitude of its normal is 1 and modifies 'd' component appropriately.

**/
virtual void   VPCALL  plane_Normalize(CPlane* pOut);
/**
\brief Normaliza o plano e armazena o resultado em pOut
\param pOut plano que armazenará o resultado
\param pIn Plano que será normalizado
\brief This function stores the normalized plane pIn into pOut, while pIn is unaffected.
**/
virtual void  VPCALL  plane_NormalizeOf(CPlane* pOut, CPlane* pIn);

virtual	bool  VPCALL  intersectLineAABB(const CAABBox &box, const CVec4D &rayPos, const CVec4D &rayDir, float tNear, float tFar);

};

} //end MATH
} //end SMF
#endif /* !__MATH_SIMD_GENERIC_H__ */
