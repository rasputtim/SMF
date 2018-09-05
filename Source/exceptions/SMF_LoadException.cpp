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
#include <string>
#include "exceptions/SMF_LoadException.h"

using namespace std;

namespace SMF {



CLoadException::CLoadException(const string & file, int line, const string & reason):
Exception::CBase(file, line),
reason(reason){
}
    
CLoadException::CLoadException(const string & file, int line, const Exception::CBase & nested, const string & reason):
Exception::CBase(file, line, nested),
reason(reason){
}

CLoadException::~CLoadException() throw (){
}

CLoadException::CLoadException(const CLoadException & copy):
Exception::CBase(copy),
reason(copy.reason){
}

Exception::CBase * CLoadException::copy() const {
    return new CLoadException(*this);
}
    
const string CLoadException::getReason() const {
    return reason;
}

} //end SMF
