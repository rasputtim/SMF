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

void Copy16( float *dst, const float *src, const int count ) {
//	Debug::debug(Debug::math, __FUNCTION__) << "size of float: " << sizeof(float) << endl;

	__asm__ (
		"mov		ecx, %0\n"
		"mov		edx, %1\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		doneCopy16\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopCopy16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"movaps	[edx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopCopy16\n"
	"doneCopy16:\n"
	:
	    : "m" (src), "m" (dst),"m"(count)
	        :);
}
__attribute__((  aligned(16)  ))  float _SIMDx86_float_SSE_NO_W_MASK[4]={0xffffffff,0xffffffff,0xffffffff,0x0};

float  Vector4D_Dot(const float* pSrc1, const float* pSrc2)
{

	/* SSE/SSE2 Implementation */
	float dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"mulps xmm1, xmm0\n"
			"movaps xmm2,xmm1\n"
			"shufps xmm2, xmm2,  0x1B\n" /*0x1B= 00011011 / xmm2 = x | y | z | w */

			"addps  xmm2, xmm1\n"			/* xmm2 = w+x | y+z | y+z | w+x */
			"movss  xmm3, xmm2\n"			/* xmm3 = ??? | ??? | ??? | w+x */
			"shufps xmm1, xmm1, 0x01\n"  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			"addss  xmm2, xmm3\n"			/* xmm2 = ??? | ??? | ??? | dot4 */


			//coloca 0-31 de xmm2 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm2\n"

	:
	: "r" (pSrc1), "r" (pSrc2),"m" (min)
	);

	return dummy;



}

void Vector3D_AlignedSumOf(float* pOut, const float* pIn1, const float* pIn2)
{
		/* SSE/SSE2/SSE3 Implementation */
		asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"addps xmm0, xmm1\n"                   // Add vector elements
                //"movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				//"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				//"movss	[ecx+ 8], xmm0\n"
                /* Store*/
                "movlpd [%2] , xmm0\n"
                "shufps  xmm0, xmm0, 0x02\n"  //shuffle colocar nos bits 0-31 o 3o byte
                "movss	[%2+ 8], xmm0\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
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
printf("Before Copy-> %f %f %f %f\n", var2[0], var2[1], var2[2], var2[3]);
Copy16(var2,var1,count);

printf("After Copy -> %f %f %f %f\n", var2[0], var2[1], var2[2], var2[3]);
Vector3D_AlignedSumOf(var3,var1,var2);
printf("After Sum -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
float x=Vector4D_Dot(var3,var1);
printf("After Dot -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
printf("Dot is -> %f \n", x);

return 0;
}



