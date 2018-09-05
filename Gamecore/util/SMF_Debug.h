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
#ifndef _SMF__DEBUG_H_
#define _SMF__DEBUG_H_


#include <ostream>
#include <string>
#include <map>
#include <vector>
#include <string>



using namespace std;
namespace SMF {

class CLoadException;


namespace Debug
{
#ifdef ANDROID
class android_ostream: public std::ostream {
public:
    android_ostream(bool enabled = true);
    static android_ostream stream;
    /* make these private at some point */
public:
    bool enabled;
    ::std::ostringstream buffer;
};

typedef android_ostream stream_type;
android_ostream & operator<<(android_ostream & stream, const string & input);
android_ostream & operator<<(android_ostream & stream, const char * input);
android_ostream & operator<<(android_ostream & stream, const char);
android_ostream & operator<<(android_ostream & stream, const double);
android_ostream & operator<<(android_ostream & stream, const int);
android_ostream & operator<<(android_ostream & stream, const short int);
android_ostream & operator<<(android_ostream & stream, const short unsigned int);
android_ostream & operator<<(android_ostream & stream, const unsigned int);
android_ostream & operator<<(android_ostream & stream, const bool);
android_ostream & operator<<(android_ostream & stream, const long int);
android_ostream & operator<<(android_ostream & stream, const unsigned long int);
android_ostream & operator<<(android_ostream & stream, const void *);
android_ostream & operator<<(android_ostream & stream, std::ostream & (*f)(std::ostream &));
#elif defined(WII) && defined(DEBUG)
class wii_ostream: public std::ostream {
public:
    wii_ostream(bool enabled = true);
    static wii_ostream stream;
    /* make these private at some point */
public:
    bool enabled;
    ::std::ostringstream buffer;
};

typedef wii_ostream stream_type;
wii_ostream & operator<<(wii_ostream & stream, const string & input);
wii_ostream & operator<<(wii_ostream & stream, const char * input);
wii_ostream & operator<<(wii_ostream & stream, const char);
wii_ostream & operator<<(wii_ostream & stream, const double);
wii_ostream & operator<<(wii_ostream & stream, const int);
wii_ostream & operator<<(wii_ostream & stream, const short int);
wii_ostream & operator<<(wii_ostream & stream, const short unsigned int);
wii_ostream & operator<<(wii_ostream & stream, const unsigned int);
wii_ostream & operator<<(wii_ostream & stream, const bool);
wii_ostream & operator<<(wii_ostream & stream, const long int);
wii_ostream & operator<<(wii_ostream & stream, const unsigned long int);
wii_ostream & operator<<(wii_ostream & stream, const void *);
wii_ostream & operator<<(wii_ostream & stream, uint64_t); 
wii_ostream & operator<<(wii_ostream & stream, std::ostream & (*f)(std::ostream &));
#elif defined(NETWORK_DEBUG)
class network_ostream: public std::ostream {
public:
    network_ostream(const string & host, int port, bool enabled = true);
    static network_ostream stream;
    /* make these private at some point */
public:
    string host;
    int port;
    bool enabled;
    ::std::ostringstream buffer;
};

typedef network_ostream stream_type;
stream_type & operator<<(stream_type & stream, const string & input);
stream_type & operator<<(stream_type & stream, const char * input);
stream_type & operator<<(stream_type & stream, const char);
stream_type & operator<<(stream_type & stream, const double);
stream_type & operator<<(stream_type & stream, const int);
stream_type & operator<<(stream_type & stream, const short int);
stream_type & operator<<(stream_type & stream, const short unsigned int);
stream_type & operator<<(stream_type & stream, const unsigned int);
stream_type & operator<<(stream_type & stream, const bool);
stream_type & operator<<(stream_type & stream, const long int);
#ifndef PS3
stream_type & operator<<(stream_type & stream, const unsigned long int);
#endif
stream_type & operator<<(stream_type & stream, const void *);
stream_type & operator<<(stream_type & stream, uint64_t); 
stream_type & operator<<(stream_type & stream, std::ostream & (*f)(std::ostream &));

#else
typedef ::std::ostream stream_type;
#endif


void setFilename(string file);


stream_type & debug(int i, const string & context = "SMF");

enum Type {
	Default,
	File
};

void setDebugMethod(Type t);
/// Determina os módulos de software disponíveis para debugar
struct DebugModule{
DebugModule():bitmap(0),
font(0),
gui(0),
input(0),
compilers(0),
configuration(0),
parsers(0),
lists(0),
menu(0),
network(0),
resources(0),
objects(0),
exceptions(0),
environment(0),
filesystem(0),
sound(0),
console(0),
gameEngine(0),
sdlmanager(0),
xml(0),
gameObject(0),
gameStage(0),
gameStateMachine(0),
gameStory(0),
gameConfiguration(0),
heap(0),
math(0),
structures(0),
error(0)
{}

bool bitmap;
bool font;
bool gui;
bool input;
bool compilers;
bool configuration;
bool parsers;
bool lists;
bool menu;
bool network;
bool resources;
bool objects;
bool exceptions;
bool environment;
bool filesystem;
bool sound;
bool console;
bool gameEngine;
bool sdlmanager;
bool xml;
bool gameObject;
bool gameStage;
bool gameStateMachine;
bool gameStory;
bool gameConfiguration;
bool heap;
bool math;
bool structures;
bool error;
};

typedef enum {
MIN = 100,
error = MIN,
bitmap,
font,
gui,
input,
compilers,
configuration,
parsers,
lists,
menu,
network,
resources,
objects,
exceptions,
environment,
filesystem,
sound,
console,
gameEngine,
gameStage,
gameObject,
gameStateMachine,
gameStory,
gameConfiguration,
sdlmanager,
xml,
heap,
math,
structures,
MAX=structures
}Modules;

/** 
\brief método que inicializa o debug no módulo especificado
\param i módulo do debug que se deseja habilitar
*/
void setDebug( Modules i );
/** 
\brief método que inicializa o debug no módulo especificado
\param i módulo do debug que se deseja habilitar
*/
void setDebug( int i );

/** 
\brief método que retorna se um módulo do debug está habilitado ou não
* \param i módulo do debug que se deseja verificar
* \return true: se o módulo está ativo
* \return false: se o módulo está inativo
*/
bool getDebug( Modules i );

/** 
\brief método que inicializa o debug em todos os módulos
*/
void setDebugAll();
/** 
\brief método que desabilita o debug em todos os módulos
\param i módulo do debug a ser retirado 
*/
void resetDebug( int i );
/** 
\brief método que desabilita o debug em todos os módulos
\param i módulo do debug a ser retirado 
*/
void resetDebug( Modules i );


int getDebug();


}
} //end SMF
#endif
