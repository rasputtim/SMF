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
#include "xmmintrin.h"
#include <malloc.h>

//#define ALIGNED                         __attribute__((  aligned(16)  ))
#define ALIGN16( x )					__attribute__((  aligned(16)  )) x
#define ALIGN4_INIT4( X, I0, I1, I2, I3 )	ALIGN16( static X[4] ) = { I0, I1, I2, I3 }


ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_W_MASK, 0xFFFFFFFF,  0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 );

//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{
//printf("Sin -> %f\n", sin);
ALIGN16( CVector3D vetor1(1.0f,2.0f,3.0f));
ALIGNED float result;
result= Vector3D_Length(vetor1);
    return 0;

}



