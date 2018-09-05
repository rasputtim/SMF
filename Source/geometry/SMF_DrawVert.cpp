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

#include "geometry/SMF_DrawVert.h"
//#pragma hdrstop
namespace SMF {

/*
=============
CVertex::Normalize
=============
*/
#if 0
void CVertex::Normalize( void ) {
	normal.toNormal();
	tangents[1].cross( normal, tangents[0] );
	tangents[1].toNormal();
	tangents[0].cross( tangents[1], normal );
	tangents[0].toNormal();
}
#endif

/*
============
CShadowVert:: createShadowCache
============
*/
int CShadowVert:: createShadowCache( CShadowVert * vertexCache, const CVertex *verts, const int numVerts ) {
	for ( int i = 0; i < numVerts; i++ ) {
		vertexCache[i*2+0].xyzw[0] = verts[i].xyz[0];
		vertexCache[i*2+0].xyzw[1] = verts[i].xyz[1];
		vertexCache[i*2+0].xyzw[2] = verts[i].xyz[2];
		vertexCache[i*2+0].xyzw[3] = 1.0f;

		vertexCache[i*2+1].xyzw[0] = verts[i].xyz[0];
		vertexCache[i*2+1].xyzw[1] = verts[i].xyz[1];
		vertexCache[i*2+1].xyzw[2] = verts[i].xyz[2];
		vertexCache[i*2+1].xyzw[3] = 0.0f;
	}
	return numVerts * 2;
}

/*
===================
CShadowVertSkinned:: createShadowCache
===================
*/
int CShadowVertSkinned:: createShadowCache( CShadowVertSkinned * vertexCache, const CVertex *verts, const int numVerts ) {
	for ( int i = 0; i < numVerts; i++ ) {
		vertexCache[0].xyzw[0] = verts[i].xyz[0];
		vertexCache[0].xyzw[1] = verts[i].xyz[1];
		vertexCache[0].xyzw[2] = verts[i].xyz[2];
		vertexCache[0].xyzw[3] = 1.0f;
		*(unsigned int *)vertexCache[0].color = *(unsigned int *)verts[i].color;
		*(unsigned int *)vertexCache[0].color2 = *(unsigned int *)verts[i].color2;

		vertexCache[1].xyzw[0] = verts[i].xyz[0];
		vertexCache[1].xyzw[1] = verts[i].xyz[1];
		vertexCache[1].xyzw[2] = verts[i].xyz[2];
		vertexCache[1].xyzw[3] = 0.0f;
		*(unsigned int *)vertexCache[1].color = *(unsigned int *)verts[i].color;
		*(unsigned int *)vertexCache[1].color2 = *(unsigned int *)verts[i].color2;

		vertexCache += 2;
	}
	return numVerts * 2;
}
} //end SMF
