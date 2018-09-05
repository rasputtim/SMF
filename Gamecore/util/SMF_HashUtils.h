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

#ifndef SMF_HASHUTILS_H_
#define SMF_HASHUTILS_H_

#include "../SMF_Config.h"

namespace SMF{
namespace Util{


int getHash( const char *string );

int getHash( const char *string, int length ) ;

int getIHash( const char *string );

int getIHash( const char *string, int length );

} //end Util
} //end SMF



#endif //!SMF_HASHUTILS_H