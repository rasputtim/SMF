//    Open Mugen is a redevelopment of Elecbyte's M.U.G.E.N wich will be 100% compatible to it
//    Copyright (C) 2004  Sahin Vardar
//
//    If you know bugs or have a wish on Open Muegn or (money/girls/a car) for me ;-)
//    Feel free and email me: sahin_v@hotmail.com  ICQ:317502935
//    Web: http://openmugen.sourceforge.net/
//    --------------------------------------------------------------------------
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.



#include <io.h>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>

__attribute__((  aligned(16)  ))  float _SIMDx86_float_SSE_NO_W_MASK[4]={0xffffffff,0xffffffff,0xffffffff,0x0};

float  Vector4D_Dot(const float* pSrc4D1, const float* pSrc4D2)
{


	float dummy;
	asm(
	"movups (%1), %%xmm0\n"
	"movups (%2), %%xmm1\n"
	"mulps %%xmm0, %%xmm1\n"
	"movaps %%xmm1, %%xmm2\n"
	"shufps $0x1B, %%xmm2, %%xmm2\n"	/* xmm2 = x | y | z | w */
	"addps %%xmm1, %%xmm2\n"			/* xmm2 = w+x | y+z | y+z | w+x */
	"movss %%xmm2, %%xmm3\n"				/* xmm3 = ??? | ??? | ??? | w+x */
	"shufps $0x01, %%xmm2, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | y+z */
	"addss %%xmm3, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | dot4 */
	"movss %%xmm2, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;



}

void Vector3D_AlignedSumOf(float* pOut, const float* pIn1, const float* pIn2)
{
/* SSE/SSE2/SSE3 Implementation */
void *SSE_NO_W_MASK=&_SIMDx86_float_SSE_NO_W_MASK;

	asm(
	"movaps (%0), %%xmm0\n"
	"movaps (%1),%%xmm1\n"
	"movaps (%3),%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"addps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%2)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "r" (SSE_NO_W_MASK)
	);


}
//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{
__attribute__((  aligned(16)  )) float var1[4];
__attribute__((  aligned(16)  )) float var2[4];
__attribute__((  aligned(16)  )) float var3[4];
var1[0]=2.0f;
var1[1]=3.0f;
var1[2]=4.0f;
var1[3]=5.0f;
var3[0]=1.0f;
var3[1]=1.0f;
var3[2]=1.0f;
var3[3]=1.0f;
int count=4;

printf("Before Sum -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
Vector3D_AlignedSumOf(var3,var1,var2);
printf("After Sum -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
float x=Vector4D_Dot(var3,var1);
printf("After Dot -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
printf("Dot is -> %f \n", x);

return 0;
}



