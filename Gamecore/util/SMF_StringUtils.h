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

#ifndef _SMF__STRING_UTILS_H__
#define _SMF__STRING_UTILS_H__

#include "SMF_ArgsCmdLine.h"
#include "../exceptions/all.h"
#include <cstring>
#include <stdio.h>

/*
===============================================================================

	Character string

===============================================================================
*/
#if 0
// these library functions should not be used for cross platform compatibility
#define strcmp			CMyString::compare		// use_idStr_Cmp
#define strncmp			use_idStr_Cmpn

#if defined( StrCmpN )
#undef StrCmpN
#endif
#define StrCmpN			use_idStr_Cmpn

#if defined( strcmpi )
#undef strcmpi
#endif
#define strcmpi			use_idStr_Icmp

#if defined( StrCmpI )
#undef StrCmpI
#endif
#define StrCmpI			use_idStr_Icmp

#if defined( StrCmpNI )
#undef StrCmpNI
#endif
#define StrCmpNI		use_idStr_Icmpn

#define stricmp			CMyString::compareInsen		// use_idStr_Icmp
#define _stricmp		use_idStr_Icmp
#define strcasecmp		use_idStr_Icmp
#define strnicmp		use_idStr_Icmpn
#define _strnicmp		use_idStr_Icmpn
#define _memicmp		use_idStr_Icmpn
#define snprintf		use_idStr_snPrintf
#define _snprintf		use_idStr_snPrintf
#define vsnprintf		use_idStr_vsnPrintf
#define _vsnprintf		use_idStr_vsnPrintf
#endif
namespace SMF {
namespace MATH{
class CVec4D;
}
namespace Util{
#ifndef FILE_HASH_SIZE
#define FILE_HASH_SIZE		1024
#endif

// color escape character
const int C_COLOR_ESCAPE			= '^';
const int C_COLOR_DEFAULT			= '0';
const int C_COLOR_RED				= '1';
const int C_COLOR_GREEN				= '2';
const int C_COLOR_YELLOW			= '3';
const int C_COLOR_BLUE				= '4';
const int C_COLOR_CYAN				= '5';
const int C_COLOR_MAGENTA			= '6';
const int C_COLOR_WHITE				= '7';
const int C_COLOR_GRAY				= '8';
const int C_COLOR_BLACK				= '9';

// color escape string
#define S_COLOR_DEFAULT				"^0"
#define S_COLOR_RED					"^1"
#define S_COLOR_GREEN				"^2"
#define S_COLOR_YELLOW				"^3"
#define S_COLOR_BLUE				"^4"
#define S_COLOR_CYAN				"^5"
#define S_COLOR_MAGENTA				"^6"
#define S_COLOR_WHITE				"^7"
#define S_COLOR_GRAY				"^8"
#define S_COLOR_BLACK				"^9"

// make CMyString a multiple of 16 bytes long
// don't make too large to keep memory requirements to a minimum
const int STR_ALLOC_BASE			= 20;
const int STR_ALLOC_GRAN			= 32;

typedef enum {
	MEASURE_SIZE = 0,
	MEASURE_BANDWIDTH
} Measure_t;

/**
 * \class CMyString
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief string c++
 * \elseif us_en
 * \brief c++ string
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMyString {

public:
						CMyString();
						CMyString( const CMyString &text );
						CMyString( const CMyString &text, int start, int end );
						CMyString( const char *text );
						CMyString( const char *text, int start, int end );
						explicit CMyString( const bool b );
						explicit CMyString( const char c );
						explicit CMyString( const int i );
						explicit CMyString( const unsigned u );
						explicit CMyString( const float f );
						~CMyString();

	size_t				size() const;
	const char *		c_str() const;
	operator			const char *() const;
	operator			const char *();

	char				operator[]( int index ) const;
	char &				operator[]( int index );

	void				operator=( const CMyString &text );
	void				operator=( const char *text );
	//void				operator=( const char c);
	void				operator=( char c);

	friend CMyString		operator+( const CMyString &a, const CMyString &b );
	friend CMyString		operator+( const CMyString &a, const char *b );
	friend CMyString		operator+( const char *a, const CMyString &b );

	friend CMyString		operator+( const CMyString &a, const float b );
	friend CMyString		operator+( const CMyString &a, const int b );
	friend CMyString		operator+( const CMyString &a, const unsigned b );
	friend CMyString		operator+( const CMyString &a, const bool b );
	friend CMyString		operator+( const CMyString &a, const char b );

	CMyString &				operator<<( const CMyString &a );
	CMyString &				operator<<( const char *a );
	CMyString &				operator<<( const float a );
	CMyString &				operator<<( const char a );
	CMyString &				operator<<( const int a );
	CMyString &				operator<<( const unsigned a );
	CMyString &				operator<<( const bool a );



	CMyString &				operator+=( const CMyString &a );
	CMyString &				operator+=( const char *a );
	CMyString &				operator+=( const float a );
	CMyString &				operator+=( const char a );
	CMyString &				operator+=( const int a );
	CMyString &				operator+=( const unsigned a );
	CMyString &				operator+=( const bool a );

	/// case sensitive compare
	friend bool			operator==( const CMyString &a, const CMyString &b );
	friend bool			operator==( const CMyString &a, const char *b );
	friend bool			operator==( const char *a, const CMyString &b );

	/// case sensitive compare
	friend bool			operator!=( const CMyString &a, const CMyString &b );
	friend bool			operator!=( const CMyString &a, const char *b );
	friend bool			operator!=( const char *a, const CMyString &b );

	/// case sensitive compare
	int					compare( const char *text ) const;
	int					compare_n( const char *text, int n ) const;
	int					comparePrefix( const char *text ) const;

	/// case insensitive compare
	int					compareInsen( const char *text ) const;
	int					compareInsen_n( const char *text, int n ) const;
	int					compareInsenPrefix( const char *text ) const;

	/// case insensitive compare ignoring color
	int					compareInsenNoColor( const char *text ) const;

	/// compares paths and makes sure folders come first
	int					comparePathInsen( const char *text ) const;
	int					comparePathInsen_n( const char *text, int n ) const;
	int					comparePathInsenPrefix( const char *text ) const;

	int					getLenght() const;
	int					allocated() const;
	void				empty();
	/// verify if the string is empty or not
	bool				isEmpty() const;
	void				clear();
	void				append( const char a );
	void				append( const CMyString &text );
	void				append( const char *text );
	void				append( const char *text, int len );
	void				insert( const char a, int index );
	void				insert( const char *text, int index );
	void				toLower();
	void				toUpper();
	bool				isNumeric() const;
	bool				isColor() const;
	bool				hasLower() const;
	bool				hasUpper() const;
	int					lengthWithoutColors() const;
	CMyString &			removeColors();
	void				capLength( int );
	void				fill( const char ch, int newlen );

	int					find( const char c, int start = 0, int end = -1 ) const;
	int					find( const char *text, bool casesensitive = true, int start = 0, int end = -1 ) const;
	bool				filter( const char *filter, bool casesensitive ) const;
	/// return the index to the last occurance of 'c', returns -1 if not found
	int					last( const char c ) const;
	/// store the leftmost 'len' characters in the result
	const char *		left( int len, CMyString &result ) const;
	/// store the rightmost 'len' characters in the result
	const char *		right( int len, CMyString &result ) const;
	/// store 'len' characters starting at 'start' in result
	const char *		mid( int start, int len, CMyString &result ) const;
	/// return the leftmost 'len' characters
	CMyString			left( int len ) const;
	/// return the rightmost 'len' characters
	CMyString			right( int len ) const;
	/// return 'len' characters starting at 'start'
	CMyString			mid( int start, int len ) const;
	/// strip char from front as many times as the char occurs
	void				stripLeading( const char c );
	/// strip string from front as many times as the string occurs
	void				stripLeading( const char *string );
	/// strip string from front just once if it occurs
	bool				stripLeadingOnce( const char *string );
	/// strip char from end as many times as the char occurs
	void				stripTrailing( const char c );
	// strip string from end as many times as the string occurs
	void				stripTrailing( const char *string );
	/// strip string from end just once if it occurs
	bool				stripTrailingOnce( const char *string );
	/// strip char from front and end as many times as the char occurs
	void				strip( const char c );
	/// strip string from front and end as many times as the string occurs
	void				strip( const char *string );
	/// strip trailing white space characters
	void				stripTrailingWhitespace();
	/// strip quotes around string
	CMyString &			stripQuotes();
	void				replace( const char *old, const char *nw );

	// file name methods
	/// hash key for the filename (skips extension)
	int					fileNameHash() const;
	/// convert slashes
	CMyString &			backSlashesToSlashes();
	/// set the given file extension
	CMyString &			setFileExtension( const char *extension );
	/// remove any file extension
	CMyString &			stripFileExtension();
	/// remove any file extension looking from front (useful if there are multiple .'s)
	CMyString &			stripAbsoluteFileExtension();
	/// if there's no file extension use the default
	CMyString &			defaultFileExtension( const char *extension );
	/// if there's no path use the default
	CMyString &			defaultPath( const char *basepath );
	/// append a partial path
	void				appendPath( const char *text );
	/// remove the filename from a path
	CMyString &			stripFilename();
	/// remove the path from the filename
	CMyString &			stripPath();
	/// copy the file path to another string
	void				extractFilePath( CMyString &dest ) const;
	// copy the filename to another string
	void				extractFileName( CMyString &dest ) const;
	// copy the filename minus the extension to another string
	void				extractFileBase( CMyString &dest ) const;
	// copy the file extension to another string
	void				extractFileExtension( CMyString &dest ) const;

	bool				checkExtension( const char *ext );

	// char * methods to replace library functions
	static int			length( const char *s );
	static char *		toLower( char *s );
	static char *		toUpper( char *s );
	static bool			isNumeric( const char *s );
	static bool			isColor( const char *s );
	static bool			hasLower( const char *s );
	static bool			hasUpper( const char *s );
	static int			lengthWithoutColors( const char *s );
	static char *		removeColors( char *s );
	static int			compare( const char *s1, const char *s2 );
	static int			compare_n( const char *s1, const char *s2, int n );
	static int			compareInsen( const char *s1, const char *s2 );
	static int			compareInsen_n( const char *s1, const char *s2, int n );
	static int			compareInsenNoColor( const char *s1, const char *s2 );
	/// compares paths and makes sure folders come first
	static int			comparePathInsen( const char *s1, const char *s2 );
	/// compares paths and makes sure folders come first
	static int			comparePathInsen_n( const char *s1, const char *s2, int n );
	static void			append( char *dest, int size, const char *src );
	static void			copyNz( char *dest, const char *src, int destsize );
	static int			snPrintf( char *dest, int size, const char *fmt, ... ) id_attribute((format(printf,3,4)));
	static int			vsnPrintf( char *dest, int size, const char *fmt, va_list argptr );
	static int			findChar( const char *str, const char c, int start = 0, int end = -1 );
	static int			findText( const char *str, const char *text, bool casesensitive = true, int start = 0, int end = -1 );
	static bool			filter( const char *filter, const char *name, bool casesensitive );
	static void			stripMediaName( const char *name, CMyString &mediaName );
	static bool			checkExtension( const char *name, const char *ext );
	static const char *	floatArraytoString( const float *array, const int length, const int precision );

	/// hash keys
	static int			hash( const char *string );
	static int			hash( const char *string, int length );
	/// case insensitive hash
	static int			hashInsens( const char *string );
	/// case insensitive hash
	static int			hashInsens( const char *string, int length );

	// character methods
	static char			toLower( char c );
	static char			toUpper( char c );
	static bool			charIsPrintable( int c );
	static bool			charIsLower( int c );
	static bool			charIsUpper( int c );
	static bool			charIsAlpha( int c );
	static bool			charIsNumeric( int c );
	static bool			charIsNewLine( char c );
	static bool			charIsTab( char c );
	static int			colorIndex( int c );
	static MATH::CVec4D &		ColorForIndex( int i );

	friend int			sprintf( CMyString &dest, const char *fmt, ... );
	friend int			vsprintf( CMyString &dest, const char *fmt, va_list ap );

	void				reAllocate( int amount, bool keepold );				// reallocate string data buffer
	void				freeData();									// free allocated string memory

						// format value in the given measurement with the best unit, returns the best unit
	int					bestUnit( const char *format, float value, Measure_t measure );
						// format value in the requested unit and measurement
	void				setUnit( const char *format, float value, int unit, Measure_t measure );

	static void			initMemory();
	static void			shutdownMemory();
	static void			purgeMemory();
	static void			showMemoryUsage_f( const CCMDLineArgs &args );

	int					dynamicMemoryUsed() const;
	static CMyString	formatNumber( int number );

protected:

	void				init();	//! initialize string using base buffer
	void				ensureAlloced( int amount, bool keepold = true );	//! ensure string data buffer is large anough

	int					m_iLen;   //! the string length
	char *				m_pData;  //! the C String
	int					m_iAllocated;  //! ammount of memory allocateb by m_pData
	char				m_cBaseBuffer[ STR_ALLOC_BASE ]; //! base buffer for string
};

char *					va( const char *fmt, ... ) id_attribute((format(printf,1,2)));

/**
\brief ensure string data buffer is large anough
**/
SMF_INLINE_FORCED void CMyString::ensureAlloced( int amount, bool keepold ) {
	if ( amount > m_iAllocated ) {
		reAllocate( amount, keepold );
	}
}
/**
\brief init the string buffer
**/
SMF_INLINE_FORCED void CMyString::init() {
	m_iLen = 0;
	m_iAllocated = STR_ALLOC_BASE;
	m_pData = m_cBaseBuffer;
	m_pData[ 0 ] = '\0';
#ifdef ID_DEBUG_UNINITIALIZED_MEMORY
	memset( baseBuffer, 0, sizeof( baseBuffer ) );
#endif
}

SMF_INLINE_FORCED CMyString::CMyString() {
	init();
}

SMF_INLINE_FORCED CMyString::CMyString( const CMyString &text ) {
	int l;

	init();
	l = text.getLenght();
	ensureAlloced( l + 1 );
	strcpy( m_pData, text.m_pData );
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::CMyString( const CMyString &text, int start, int end ) {
	int i;
	int l;

	init();
	if ( end > text.getLenght() ) {
		end = text.getLenght();
	}
	if ( start > text.getLenght() ) {
		start = text.getLenght();
	} else if ( start < 0 ) {
		start = 0;
	}

	l = end - start;
	if ( l < 0 ) {
		l = 0;
	}

	ensureAlloced( l + 1 );

	for ( i = 0; i < l; i++ ) {
		m_pData[ i ] = text[ start + i ];
	}

	m_pData[ l ] = '\0';
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::CMyString( const char *text ) {
	int l;

	init();
	if ( text ) {

#ifdef _WIN32
	l = std::strlen( text );
#else
	l = strlen( text );
#endif
		ensureAlloced( l + 1 );
		strcpy( m_pData, text );
		m_iLen = l;
	}
}

SMF_INLINE_FORCED CMyString::CMyString( const char *text, int start, int end ) {
	int i;
#ifdef _WIN32
	int l = std::strlen( text );
#else
	int l = strlen( text );
#endif
	init();
	if ( end > l ) {
		end = l;
	}
	if ( start > l ) {
		start = l;
	} else if ( start < 0 ) {
		start = 0;
	}

	l = end - start;
	if ( l < 0 ) {
		l = 0;
	}

	ensureAlloced( l + 1 );

	for ( i = 0; i < l; i++ ) {
		m_pData[ i ] = text[ start + i ];
	}

	m_pData[ l ] = '\0';
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::CMyString( const bool b ) {
	init();
	ensureAlloced( 2 );
	m_pData[ 0 ] = b ? '1' : '0';
	m_pData[ 1 ] = '\0';
	m_iLen = 1;
}

SMF_INLINE_FORCED CMyString::CMyString( const char c ) {
	init();
	ensureAlloced( 2 );
	m_pData[ 0 ] = c;
	m_pData[ 1 ] = '\0';
	m_iLen = 1;
}

SMF_INLINE_FORCED CMyString::CMyString( const int i ) {
	char text[ 64 ];
	int l;

	init();
#if defined(_MSC_VER)
	l = std::sprintf( text, "%d", i );
#else
	l = sprintf( text, "%d", i );
#endif

	ensureAlloced( l + 1 );
	strcpy( m_pData, text );
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::CMyString( const unsigned u ) {
	char text[ 64 ];
	int l;

	init();
#if defined(_MSC_VER)
	l = std::sprintf( text, "%u", u );
#else
	l = sprintf( text, "%u", u );
#endif
	ensureAlloced( l + 1 );
	strcpy( m_pData, text );
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::CMyString( const float f ) {
	char text[ 64 ];
	int l;

	init();
	l = CMyString::snPrintf( text, sizeof( text ), "%f", f );
	while( l > 0 && text[l-1] == '0' ) text[--l] = '\0';
	while( l > 0 && text[l-1] == '.' ) text[--l] = '\0';
	ensureAlloced( l + 1 );
	strcpy( m_pData, text );
	m_iLen = l;
}

SMF_INLINE_FORCED CMyString::~CMyString() {
	freeData();
}

SMF_INLINE_FORCED size_t CMyString::size() const {
	return sizeof( *this ) + allocated();
}

SMF_INLINE_FORCED const char *CMyString::c_str() const {
	return m_pData;
}

SMF_INLINE_FORCED CMyString::operator const char *() {
	return c_str();
}

SMF_INLINE_FORCED CMyString::operator const char *() const {
	return c_str();
}

SMF_INLINE_FORCED char CMyString::operator[]( int index ) const {
	SMF_ASSERT( ( index >= 0 ) && ( index <= m_iLen ) );
	return m_pData[ index ];
}

SMF_INLINE_FORCED char &CMyString::operator[]( int index ) {
	SMF_ASSERT( ( index >= 0 ) && ( index <= m_iLen ) );
	return m_pData[ index ];
}

SMF_INLINE_FORCED void CMyString::operator=( const CMyString &text ) {
	int l;

	l = text.getLenght();
	ensureAlloced( l + 1, false );
	memcpy( m_pData, text.m_pData, l );
	m_pData[l] = '\0';
	m_iLen = l;
}
SMF_INLINE_FORCED void	CMyString::operator=( char c){
	ensureAlloced( 2 );
	m_pData[ 0 ] = c;
	m_pData[ 1 ] = '\0';
	m_iLen = 1;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const CMyString &b ) {
	CMyString result( a );
	result.append( b );
	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const char *b ) {
	CMyString result( a );
	result.append( b );
	return result;
}

SMF_INLINE_FORCED CMyString operator+( const char *a, const CMyString &b ) {
	CMyString result( a );
	result.append( b );
	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const bool b ) {
	CMyString result( a );
	result.append( b ? "true" : "false" );
	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const char b ) {
	CMyString result( a );
	result.append( b );
	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const float b ) {
	char	text[ 64 ];
	CMyString	result( a );
#if defined(_MSC_VER)
	std::sprintf( text, "%f", b );
#else
	sprintf( text, "%f", b );
#endif
	result.append( text );

	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const int b ) {
	char	text[ 64 ];
	CMyString	result( a );
#if defined(_MSC_VER)
	std::sprintf( text, "%d", b );
#else
	sprintf( text, "%d", b );
#endif

	result.append( text );

	return result;
}

SMF_INLINE_FORCED CMyString operator+( const CMyString &a, const unsigned b ) {
	char	text[ 64 ];
	CMyString	result( a );
#if defined(_MSC_VER)
	std::sprintf( text, "%u", b );
#else
	sprintf( text, "%u", b );
#endif

	result.append( text );

	return result;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const float a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%f", a );
#else
	sprintf( text, "%f", a );
#endif

	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const int a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%d", a );
#else
	sprintf( text, "%d", a );
#endif

	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const unsigned a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%u", a );
#else
	sprintf( text, "%u", a );
#endif
	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const CMyString &a ) {
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const char *a ) {
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const char a ) {
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator+=( const bool a ) {
	append( a ? "true" : "false" );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const CMyString &a ){
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const float a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%f", a );
#else
	sprintf( text, "%f", a );
#endif

	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const int a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%d", a );
#else
	sprintf( text, "%d", a );
#endif

	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const unsigned a ) {
	char text[ 64 ];
#if defined(_MSC_VER)
	std::sprintf( text, "%u", a );
#else
	sprintf( text, "%u", a );
#endif
	append( text );

	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const char *a ) {
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const char a ) {
	append( a );
	return *this;
}

SMF_INLINE_FORCED CMyString &CMyString::operator<<( const bool a ) {
	append( a ? "true" : "false" );
	return *this;
}


SMF_INLINE_FORCED bool operator==( const CMyString &a, const CMyString &b ) {
	return ( !CMyString::compare( a.m_pData, b.m_pData ) );
}

SMF_INLINE_FORCED bool operator==( const CMyString &a, const char *b ) {
	SMF_ASSERT( b );
	return ( !CMyString::compare( a.m_pData, b ) );
}

SMF_INLINE_FORCED bool operator==( const char *a, const CMyString &b ) {
	SMF_ASSERT( a );
	return ( !CMyString::compare( a, b.m_pData ) );
}

SMF_INLINE_FORCED bool operator!=( const CMyString &a, const CMyString &b ) {
	return !( a == b );
}

SMF_INLINE_FORCED bool operator!=( const CMyString &a, const char *b ) {
	return !( a == b );
}

SMF_INLINE_FORCED bool operator!=( const char *a, const CMyString &b ) {
	return !( a == b );
}

SMF_INLINE_FORCED int CMyString::compare( const char *text ) const {
	SMF_ASSERT( text );
	return CMyString::compare( m_pData, text );
}

SMF_INLINE_FORCED int CMyString::compare_n( const char *text, int n ) const {
	SMF_ASSERT( text );
	return CMyString::compare_n( m_pData, text, n );
}

SMF_INLINE_FORCED int CMyString::comparePrefix( const char *text ) const {
	SMF_ASSERT( text );
#ifdef WIN32
	return CMyString::compare_n( m_pData, text, std::strlen( text ) );
#else
	return CMyString::compare_n( m_pData, text, strlen( text ) );
#endif

}

SMF_INLINE_FORCED int CMyString::compareInsen( const char *text ) const {
	SMF_ASSERT( text );
	return CMyString::compareInsen( m_pData, text );
}

SMF_INLINE_FORCED int CMyString::compareInsen_n( const char *text, int n ) const {
	SMF_ASSERT( text );
	return CMyString::compareInsen_n( m_pData, text, n );
}

SMF_INLINE_FORCED int CMyString::compareInsenPrefix( const char *text ) const {
	SMF_ASSERT( text );
#ifdef WIN32
	return CMyString::compareInsen_n( m_pData, text, std::strlen( text ) );
#else
	return CMyString::compareInsen_n( m_pData, text, strlen( text ) );
#endif

}

SMF_INLINE_FORCED int CMyString::compareInsenNoColor( const char *text ) const {
	SMF_ASSERT( text );
	return CMyString::compareInsenNoColor( m_pData, text );
}

SMF_INLINE_FORCED int CMyString::comparePathInsen( const char *text ) const {
	SMF_ASSERT( text );
	return CMyString::comparePathInsen( m_pData, text );
}

SMF_INLINE_FORCED int CMyString::comparePathInsen_n( const char *text, int n ) const {
	SMF_ASSERT( text );
	return CMyString::comparePathInsen_n( m_pData, text, n );
}

SMF_INLINE_FORCED int CMyString::comparePathInsenPrefix( const char *text ) const {
	SMF_ASSERT( text );
#ifdef WIN32
	return CMyString::comparePathInsen_n( m_pData, text, std::strlen( text ) );
#else
	return CMyString::comparePathInsen_n( m_pData, text, strlen( text ) );
#endif

}

SMF_INLINE_FORCED int CMyString::getLenght() const {
	return m_iLen;
}

SMF_INLINE_FORCED int CMyString::allocated() const {
	if ( m_pData != m_cBaseBuffer ) {
		return m_iAllocated;
	} else {
		return 0;
	}
}
//! clean the String
SMF_INLINE_FORCED void CMyString::empty() {
	ensureAlloced( 1 );
	m_pData[ 0 ] = '\0';
	m_iLen = 0;
}
//! verify if the string is empty or not
SMF_INLINE_FORCED bool CMyString::isEmpty() const {
	return ( CMyString::compare( m_pData, "" ) == 0 );
}

SMF_INLINE_FORCED void CMyString::clear() {
	freeData();
	init();
}

SMF_INLINE_FORCED void CMyString::append( const char a ) {
	ensureAlloced( m_iLen + 2 );
	m_pData[ m_iLen ] = a;
	m_iLen++;
	m_pData[ m_iLen ] = '\0';
}

SMF_INLINE_FORCED void CMyString::append( const CMyString &text ) {
	int newLen;
	int i;

	newLen = m_iLen + text.getLenght();
	ensureAlloced( newLen + 1 );
	for ( i = 0; i < text.m_iLen; i++ ) {
		m_pData[ m_iLen + i ] = text[ i ];
	}
	m_iLen = newLen;
	m_pData[ m_iLen ] = '\0';
}

SMF_INLINE_FORCED void CMyString::append( const char *text ) {
	int newLen;
	int i;

	if ( text ) {
#ifdef WIN32
	newLen = m_iLen + std::strlen( text );
#else
	newLen = m_iLen + strlen( text );
#endif

		ensureAlloced( newLen + 1 );
		for ( i = 0; text[ i ]; i++ ) {
			m_pData[ m_iLen + i ] = text[ i ];
		}
		m_iLen = newLen;
		m_pData[ m_iLen ] = '\0';
	}
}

SMF_INLINE_FORCED void CMyString::append( const char *text, int l ) {
	int newLen;
	int i;

	if ( text && l ) {
		newLen = m_iLen + l;
		ensureAlloced( newLen + 1 );
		for ( i = 0; text[ i ] && i < l; i++ ) {
			m_pData[ m_iLen + i ] = text[ i ];
		}
		m_iLen = newLen;
		m_pData[ m_iLen ] = '\0';
	}
}

SMF_INLINE_FORCED void CMyString::insert( const char a, int index ) {
	int i, l;

	if ( index < 0 ) {
		index = 0;
	} else if ( index > m_iLen ) {
		index = m_iLen;
	}

	l = 1;
	ensureAlloced( m_iLen + l + 1 );
	for ( i = m_iLen; i >= index; i-- ) {
		m_pData[i+l] = m_pData[i];
	}
	m_pData[index] = a;
	m_iLen++;
}

SMF_INLINE_FORCED void CMyString::insert( const char *text, int index ) {
	int i, l;

	if ( index < 0 ) {
		index = 0;
	} else if ( index > m_iLen ) {
		index = m_iLen;
	}
#ifdef WIN32
	l = std::strlen( text );
#else
	l = strlen( text );
#endif


	ensureAlloced( m_iLen + l + 1 );
	for ( i = m_iLen; i >= index; i-- ) {
		m_pData[i+l] = m_pData[i];
	}
	for ( i = 0; i < l; i++ ) {
		m_pData[index+i] = text[i];
	}
	m_iLen += l;
}

SMF_INLINE_FORCED void CMyString::toLower() {
	for (int i = 0; m_pData[i]; i++ ) {
		if ( charIsUpper( m_pData[i] ) ) {
			m_pData[i] += ( 'a' - 'A' );
		}
	}
}

SMF_INLINE_FORCED void CMyString::toUpper() {
	for (int i = 0; m_pData[i]; i++ ) {
		if ( charIsLower( m_pData[i] ) ) {
			m_pData[i] -= ( 'a' - 'A' );
		}
	}
}

SMF_INLINE_FORCED bool CMyString::isNumeric() const {
	return CMyString::isNumeric( m_pData );
}

SMF_INLINE_FORCED bool CMyString::isColor() const {
	return CMyString::isColor( m_pData );
}

SMF_INLINE_FORCED bool CMyString::hasLower() const {
	return CMyString::hasLower( m_pData );
}

SMF_INLINE_FORCED bool CMyString::hasUpper() const {
	return CMyString::hasUpper( m_pData );
}

SMF_INLINE_FORCED CMyString &CMyString::removeColors() {
	CMyString::removeColors( m_pData );
	m_iLen = length( m_pData );
	return *this;
}

SMF_INLINE_FORCED int CMyString::lengthWithoutColors() const {
	return CMyString::lengthWithoutColors( m_pData );
}

SMF_INLINE_FORCED void CMyString::capLength( int newlen ) {
	if ( m_iLen <= newlen ) {
		return;
	}
	m_pData[ newlen ] = 0;
	m_iLen = newlen;
}

SMF_INLINE_FORCED void CMyString::fill( const char ch, int newlen ) {
	ensureAlloced( newlen + 1 );
	m_iLen = newlen;
	memset( m_pData, ch, m_iLen );
	m_pData[ m_iLen ] = 0;
}

SMF_INLINE_FORCED int CMyString::find( const char c, int start, int end ) const {
	if ( end == -1 ) {
		end = m_iLen;
	}
	return CMyString::findChar( m_pData, c, start, end );
}

SMF_INLINE_FORCED int CMyString::find( const char *text, bool casesensitive, int start, int end ) const {
	if ( end == -1 ) {
		end = m_iLen;
	}
	return CMyString::findText( m_pData, text, casesensitive, start, end );
}

SMF_INLINE_FORCED bool CMyString::filter( const char *filter, bool casesensitive ) const {
	return CMyString::filter( filter, m_pData, casesensitive );
}

SMF_INLINE_FORCED const char *CMyString::left( int m_iLen, CMyString &result ) const {
	return mid( 0, m_iLen, result );
}

SMF_INLINE_FORCED const char *CMyString::right( int len, CMyString &result ) const {
	if ( len >= getLenght() ) {
		result = *this;
		return result;
	}
	return mid( getLenght() - len, len, result );
}

SMF_INLINE_FORCED CMyString CMyString::left( int len ) const {
	return mid( 0, len );
}

SMF_INLINE_FORCED CMyString CMyString::right( int len ) const {
	if ( len >= getLenght() ) {
		return *this;
	}
	return mid( getLenght() - len, len );
}

SMF_INLINE_FORCED void CMyString::strip( const char c ) {
	stripLeading( c );
	stripTrailing( c );
}

SMF_INLINE_FORCED void CMyString::strip( const char *string ) {
	stripLeading( string );
	stripTrailing( string );
}

SMF_INLINE_FORCED bool CMyString::checkExtension( const char *ext ) {
	return CMyString::checkExtension( m_pData, ext );
}

SMF_INLINE_FORCED int CMyString::length( const char *s ) {
	int i;
	for ( i = 0; s[i]; i++ ) {}
	return i;
}

SMF_INLINE_FORCED char *CMyString::toLower( char *s ) {
	for ( int i = 0; s[i]; i++ ) {
		if ( charIsUpper( s[i] ) ) {
			s[i] += ( 'a' - 'A' );
		}
	}
	return s;
}

SMF_INLINE_FORCED char *CMyString::toUpper( char *s ) {
	for ( int i = 0; s[i]; i++ ) {
		if ( charIsLower( s[i] ) ) {
			s[i] -= ( 'a' - 'A' );
		}
	}
	return s;
}

SMF_INLINE_FORCED int CMyString::hash( const char *string ) {
	int i, hash = 0;
	for ( i = 0; *string != '\0'; i++ ) {
		hash += ( *string++ ) * ( i + 119 );
	}
	return hash;
}

SMF_INLINE_FORCED int CMyString::hash( const char *string, int length ) {
	int i, hash = 0;
	for ( i = 0; i < length; i++ ) {
		hash += ( *string++ ) * ( i + 119 );
	}
	return hash;
}

SMF_INLINE_FORCED int CMyString::hashInsens( const char *string ) {
	int i, hash = 0;
	for( i = 0; *string != '\0'; i++ ) {
		hash += toLower( *string++ ) * ( i + 119 );
	}
	return hash;
}

SMF_INLINE_FORCED int CMyString::hashInsens( const char *string, int length ) {
	int i, hash = 0;
	for ( i = 0; i < length; i++ ) {
		hash += toLower( *string++ ) * ( i + 119 );
	}
	return hash;
}

SMF_INLINE_FORCED bool CMyString::isColor( const char *s ) {
	return ( s[0] == C_COLOR_ESCAPE && s[1] != '\0' && s[1] != ' ' );
}

SMF_INLINE_FORCED char CMyString::toLower( char c ) {
	if ( c <= 'Z' && c >= 'A' ) {
		return ( c + ( 'a' - 'A' ) );
	}
	return c;
}

SMF_INLINE_FORCED char CMyString::toUpper( char c ) {
	if ( c >= 'a' && c <= 'z' ) {
		return ( c - ( 'a' - 'A' ) );
	}
	return c;
}

SMF_INLINE_FORCED bool CMyString::charIsPrintable( int c ) {
	// test for regular ascii and western European high-ascii chars
	return ( c >= 0x20 && c <= 0x7E ) || ( c >= 0xA1 && c <= 0xFF );
}

SMF_INLINE_FORCED bool CMyString::charIsLower( int c ) {
	// test for regular ascii and western European high-ascii chars
	return ( c >= 'a' && c <= 'z' ) || ( c >= 0xE0 && c <= 0xFF );
}

SMF_INLINE_FORCED bool CMyString::charIsUpper( int c ) {
	// test for regular ascii and western European high-ascii chars
	return ( c <= 'Z' && c >= 'A' ) || ( c >= 0xC0 && c <= 0xDF );
}

SMF_INLINE_FORCED bool CMyString::charIsAlpha( int c ) {
	// test for regular ascii and western European high-ascii chars
	return ( ( c >= 'a' && c <= 'z' ) || ( c >= 'A' && c <= 'Z' ) ||
			 ( c >= 0xC0 && c <= 0xFF ) );
}

SMF_INLINE_FORCED bool CMyString::charIsNumeric( int c ) {
	return ( c <= '9' && c >= '0' );
}

SMF_INLINE_FORCED bool CMyString::charIsNewLine( char c ) {
	return ( c == '\n' || c == '\r' || c == '\v' );
}

SMF_INLINE_FORCED bool CMyString::charIsTab( char c ) {
	return ( c == '\t' );
}

SMF_INLINE_FORCED int CMyString::colorIndex( int c ) {
	return ( c & 15 );
}

SMF_INLINE_FORCED int CMyString::dynamicMemoryUsed() const {
	return ( m_pData == m_cBaseBuffer ) ? 0 : m_iAllocated;
}
} //end Util
} // end SMF

#endif /* !__SMF_STRING_UTILS_H__ */
