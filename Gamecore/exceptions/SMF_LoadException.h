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
#ifndef _SMF__load_exception_h
#define _SMF__load_exception_h

/* Generic load exception class. Thrown whenever a structure is being created
 * and an error occurs.
 */


#include <string>
#include "../SMF_Config.h"
#include "SMF_Exception.h"

using namespace std;
namespace SMF {
	namespace Exception {
		/**
		 * \class CLoadException
		 *
		 * \ingroup SMF_Exceptions
		 *
		 * \brief Excessão de carregamento
		 *
		 * Ocorre sempre que uma estrutura está sendo criada e ocorre um erro, ex: carregamento de arquivos, configurações, etc
		 *
		 *
		 * \author (last to touch it) $Autor: Rasputtim $
		 *
		 * \version 1.0 $Revision: 1.0 $
		 *
		 * Contact: Rasputtim@hotmail.com
		 *
		 * Created on: 04 de Janeiro de 2012
		  */

		class  SMF_API CLoadException : public Exception::CBase {
		public:

			CLoadException(const string reason);
			CLoadException(const string & file, int line, const string & reason);
			CLoadException(const string & file, int line, const Exception::CBase & nested, const string & reason);
			CLoadException(const CLoadException & copy);
			virtual void throwSelf() const {
				throw *this;
			}



			virtual ~CLoadException() throw();
			virtual const string getReason() const;  //temporariamente aqui até implementar catch do exception base
		protected:
			string reason;
			virtual Exception::CBase * copy() const;



		};
	} //end exception
} //end SMF
#endif
