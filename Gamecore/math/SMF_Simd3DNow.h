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

#ifndef _SMF__MATH_SIMD_3DNOW_H__
#define _SMF__MATH_SIMD_3DNOW_H__
#include "../SMF_Config.h"
#include "SMF_SimdMMX.h"
namespace SMF {
namespace MATH{
#if !defined (__GNUC__)
/*
===============================================================================

	3DNow! implementation of CSIMDProcessor
	http://en.wikipedia.org/wiki/3DNow!
===============================================================================
*/
/**
 * \class CSIMD_3DNow
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação do processador 3DNow!
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \note http://en.wikipedia.org/wiki/3DNow!
 **/
class SMF_API CSIMD_3DNow : public CSIMD_MMX {
public:
	CSIMD_3DNow():CSIMD_MMX(CPUID_3DNOW){};
	CSIMD_3DNow(cpuid_t id):CSIMD_MMX(id){};
	~CSIMD_3DNow(){};
#if defined (WIN32) || defined (__GNUC__)
public:
	virtual const char * VPCALL getName() const;

	virtual void VPCALL memCopy( void *dst,			const void *src,		const int count );
//=================================
//métodos GNU C
#if 0
virtual void VPCALL vector3D_Sum(CVec3D* pOut, const CVec3D* pIn);
virtual void VPCALL vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2);
virtual void VPCALL vector3D_Diff(CVec3D* pLeft, CVec3D* pRight);
virtual void VPCALL vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_Scale(CVec3D* pOut, float scalar);
virtual void VPCALL vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar);
virtual float VPCALL vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2);
virtual float VPCALL vector3D_Dot4(const float* pSrc4D1, const float* pSrc4D2);
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
virtual float VPCALL vector3D_AlignedDot4(const float* pSrc4D1, const float* pSrc4D2);
virtual float VPCALL vector3D_AlignedLengthSq(const CVec3D* pVec);
virtual float VPCALL vector3D_AlignedLength(const CVec3D* pVec);
virtual void VPCALL vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_AlignedNormalize(CVec3D* pVec);
virtual void VPCALL vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec);
virtual float VPCALL vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2);
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

#endif

#endif
};
#endif //__GNUC__
} //end MATH
}// end SMF
#endif /* !__MATH_SIMD_3DNOW_H__ */
