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

#include "exceptions/SMF_Error.h"
#include <sstream>


namespace SMF {
CGeneralError::CGeneralError():
Exception::CBase("", -1),
reason("unspecified"),
where("?"),
line(0){
}

CGeneralError::CGeneralError(const string & reason, const string & where, int line):
Exception::CBase("", line),
reason(reason),
where(where),
line(line){
	
}
        
Exception::CBase * CGeneralError::copy() const {
    return new CGeneralError(reason, where, line);
}
        
const string CGeneralError::getFullReason() const {
    ::std::ostringstream out;
    out << where << ":" << line << " " << reason;
    return out.str();
}
	
CGeneralError::~CGeneralError() throw() {
}

} //end SMF