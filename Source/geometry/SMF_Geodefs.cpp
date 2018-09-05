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

#include "SMF_Config.h"
#include "geometry/SMF_GeoDefs.h"

namespace SMF {
namespace GEO{
static const int RightHandSide        = -1;
static const int LeftHandSide         = +1;
static const int Clockwise            = -1;
static const int CounterClockwise     = +1;
static const int CW = CounterClockwise; //same as Clockwise
static const int CCW = Clockwise; //same as CounterClockwise
static const int CollinearOrientation =  0;
static const int AboveOrientation     = +1;
static const int BelowOrientation     = -1;
static const int CoplanarOrientation  =  0;

} //end GEO
}  //end SMF
