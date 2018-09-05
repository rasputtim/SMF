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

#include "sys/SMF_System.h"
#include "SMF_Config.h"
#ifdef _WIN32

#include <windows.h>
#include <fstream>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added
//#include "sys/win32/win_public.h"
#include  <io.h>
using namespace std;
#endif

#include <string>
#include  <stdio.h>
#include  <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include "sys/SMF_System.h"
#include "util/SMF_Debug.h"
#include "util/SMF_StringUtils.h"


namespace SMF{
using namespace Util;
namespace System{

#if defined  (WIN32) && defined(_MSC_VER)
/*
http://msdn.microsoft.com/en-us/library/windows/desktop/ms684877%28v=vs.85%29.aspx
*/
unsigned long memoryUsage(){
    HANDLE id = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS info;
    BOOL okay = GetProcessMemoryInfo(id, &info, sizeof(info));
    if (okay){
        return info.WorkingSetSize;
    }
    return 0;
}
#else
unsigned long memoryUsage(){
    return 0;
}

#endif
#ifndef _WIN32
//#ifndef WINDOWS
//#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
 #include <sys/time.h>  
 #include <stdio.h>  

uint64_t CTime::currentMicroseconds(){
    struct timeval hold;
    gettimeofday(&hold, NULL);
    return hold.tv_sec * 1000 * 1000 + hold.tv_usec;
}
#else

uint64_t CTime::currentMicroseconds(){
    LARGE_INTEGER ticksPerSecond;
    LARGE_INTEGER tick;
    QueryPerformanceFrequency(&ticksPerSecond);
    QueryPerformanceCounter(&tick);
    return (tick.QuadPart)/(ticksPerSecond.QuadPart/1000000);
}
#endif
uint64_t CTime::currentSeconds(){
    return currentMicroseconds() / (1000 * 1000);
}





/*
=================
timeStampToStr
=================
*/
const char *timeStampToStr( time_t timeStamp ) {
	static char timeString[MAX_STRING_CHARS];
	timeString[0] = '\0';

	tm*	time = localtime( &timeStamp );
	CMyString out;
	
	CMyString lang = "english";//s cvarSystem->GetCVarString( "sys_lang" );
	if ( lang.compareInsen( "english" ) == 0 ) {
		// english gets "month/day/year  hour:min" + "am" or "pm"
		out = va( "%02d", time->tm_mon + 1 );
		out += "/";
		out += va( "%02d", time->tm_mday );
		out += "/";
		out += va( "%d", time->tm_year + 1900 );
		out += "\t";
		if ( time->tm_hour > 12 ) {
			out += va( "%02d", time->tm_hour - 12 );
		} else if ( time->tm_hour == 0 ) {
				out += "12";
		} else {
			out += va( "%02d", time->tm_hour );
		}
		out += ":";
		out +=va( "%02d", time->tm_min );
		if ( time->tm_hour >= 12 ) {
			out += "pm";
		} else {
			out += "am";
		}
	} else {
		// europeans get "day/month/year  24hour:min"
		out = va( "%02d", time->tm_mday );
		out += "/";
		out += va( "%02d", time->tm_mon + 1 );
		out += "/";
		out += va( "%d", time->tm_year + 1900 );
		out += "\t";
		out += va( "%02d", time->tm_hour );
		out += ":";
		out += va( "%02d", time->tm_min );
	}
	SMF:: CMyString::copyNz( timeString, out, sizeof( timeString ) );

	return timeString;
}



#undef std::sprintf

} //end System
} //end SMF
	
	
