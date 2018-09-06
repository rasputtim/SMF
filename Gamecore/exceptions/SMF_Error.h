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

#ifndef _SMF__error_12345_h
#define _SMF__error_12345_h

#include "../SMF_Config.h"
#include "SMF_Exception.h"
#include <string>

using namespace std;
namespace SMF {
	namespace Exception {
		/**
		 * \class CGeneralError
		 *
		 * \ingroup SMF_Exceptions
		 *
		 * \brief Excessão de Erro Genérico
		 *
		 * \author (last to touch it) $Autor: Rasputtim $
		 *
		 * \version 1.0 $Revision: 1.0 $
		 *
		 * Contact: Rasputtim@hotmail.com
		 *
		 * Created on: 04 de Janeiro de 2012
		  */

		class  SMF_API CGeneralError : public CBase {
		public:
			CGeneralError();
			CGeneralError(const string & reason, const string & where = "?", int line = 0);

			virtual ~CGeneralError() throw();

			virtual CBase * copy() const;

			const string getFullReason() const;

			virtual void throwSelf() const {
				throw *this;
			}

			inline const string getReason() const {
				return reason;
			}

		protected:
			string reason;
			string where;
			int line;
		};
	} //end exception
} //end SMF
#endif