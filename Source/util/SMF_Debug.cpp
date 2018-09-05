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

#include <iostream>
#include <fstream>
#include <sstream>
#include "util/SMF_Debug.h"
#include "exceptions/all.h"

using namespace std;
namespace SMF {


class nullstreambuf_t: public std::streambuf {
public:
	nullstreambuf_t():std::streambuf(){
	}
};

static nullstreambuf_t nullstreambuf;

class nullcout_t: public std::ostream {
public:
	nullcout_t():std::ostream(&nullstreambuf){
	}
};

static nullcout_t nullcout;



namespace Debug{
static int global_debug_level = 0;
static ofstream  debug_file;
static string debug_filename;
static string debug_parsers_filename;
static map<int,streambuf*> debug_buffers;
static DebugModule whatDebug;

void setDebugMethod(Type t){
	switch(t) {
	case File: {
				for( int x=MIN; x <= MAX; x++ ) {
				debug_buffers[x]= debug_file.rdbuf();
				}
				break;
	}
	case Default: {
				for( int x=MIN; x <= MAX; x++ ) {
					debug_buffers[x]= cout.rdbuf();
					
				}
				debug_buffers[0]=cout.rdbuf();
				break;
	}
	}
}

void setFilename(string filename){
	    debug_filename=filename;
		debug_file.open (filename.c_str());
		debug_buffers[0]= cout.rdbuf();		
}



std::ostream & debug(int i, const string & context){
	bool debug=false;

	switch (i) {
	case error: (whatDebug.error ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[error]); break;
	case bitmap: (whatDebug.bitmap ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[bitmap]); break;
	case font: (whatDebug.font ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[font]); break;
	case gui: (whatDebug.gui ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gui]);  break;
	case input: (whatDebug.input ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[input]); break;
	case compilers: (whatDebug.compilers ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[compilers]); break;
	case configuration: (whatDebug.configuration ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[configuration]); break;
	case parsers: (whatDebug.parsers ? debug = 1 : debug = 0); cout.rdbuf(debug_buffers[parsers]); break;
	case lists: (whatDebug.lists ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[lists]); break;
	case menu: (whatDebug.menu ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[menu]); break;
	case network: (whatDebug.network ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[network]); break;
	case resources: (whatDebug.resources ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[resources]); break;
	case objects: (whatDebug.objects ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[objects]); break;
	case exceptions: (whatDebug.exceptions ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[exceptions]); break;
	case environment: (whatDebug.environment ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[environment]); break;
	case filesystem: (whatDebug.filesystem ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[filesystem]); break;
	case sound: (whatDebug.sound ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[sound]); break;
	case console: (whatDebug.console ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[console]); break;
	case gameEngine: (whatDebug.gameEngine ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameEngine]); break;
	case gameStage: (whatDebug.gameStage ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameStage]); break;
	case gameObject: (whatDebug.gameObject ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameObject]); break;
	case gameStateMachine: (whatDebug.gameStateMachine ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameStateMachine]); break;
	case gameStory: (whatDebug.gameStory ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameStory]); break;
	case gameConfiguration: (whatDebug.gameConfiguration ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[gameConfiguration]); break;
	case sdlmanager: (whatDebug.sdlmanager ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[sdlmanager]); break;
	case xml: (whatDebug.xml ? debug = 1 : debug = 0);cout.rdbuf(debug_buffers[xml]); break;
	case heap: (whatDebug.heap ? debug = 1 : debug = 0); cout.rdbuf(debug_buffers[heap]); break;
	case math: (whatDebug.math ? debug = 1 : debug = 0); cout.rdbuf(debug_buffers[math]); break;
	case structures: (whatDebug.structures ? debug = 1 : debug = 0); cout.rdbuf(debug_buffers[structures]); break;
	default: debug = 0;cout.rdbuf(debug_buffers[xml]); break;
	}

  
	if (debug) { //0 means no debug
	 /*	switch (i) {
		//case parsers:
		//	 cout.rdbuf(debug_buffers[2]);

		default:
			 cout.rdbuf(debug_buffers[1]);
			
		} */
		//if(i!=parsers)  cout.rdbuf(debug_buffers[1]);
		cout << "[" << i << ":" << context << "] ";
		return cout;
	}
	return nullcout;
}
void setDebugAll() {
	cout << "Setting debug: " << endl;
	//Modules e;
	int x;
	int xmax = MAX;
	for(  x=MIN; x <= xmax; x++ ) {
		cout << "Setting debug: " << x << endl;
		//debug_buffers[x]= cout.rdbuf();
		setDebug( x );
	}



}

void setDebug( Debug::Modules i ){

	switch (i) {
	case error: (whatDebug.error=true);break;
	case bitmap: (whatDebug.bitmap=true);break;
	case font: (whatDebug.font=true);break;
	case gui: (whatDebug.gui=true); break;
	case input: (whatDebug.input=true);break;
	case compilers: (whatDebug.compilers=true);break;
	case configuration: (whatDebug.configuration=true);break;
	case parsers: (whatDebug.parsers=true);break;
	case lists: (whatDebug.lists=true);break;
	case menu: (whatDebug.menu=true);break;
	case network: (whatDebug.network=true);break;
	case resources: (whatDebug.resources=true);break;
	case objects: (whatDebug.objects=true);break;
	case exceptions: (whatDebug.exceptions=true);break;
	case environment: (whatDebug.environment=true);break;
	case filesystem: (whatDebug.filesystem=true);break;
	case sound: (whatDebug.sound=true);break;
	case console: (whatDebug.console=true);break;
	case gameEngine: (whatDebug.gameEngine=true);break;
	case gameStage: (whatDebug.gameStage=true);	break;
	case gameObject: (whatDebug.gameObject=true);break;
	case gameStateMachine: (whatDebug.gameStateMachine=true);break;
	case gameStory: (whatDebug.gameStory=true);break;
	case gameConfiguration: (whatDebug.gameConfiguration=true);break;
	case sdlmanager: (whatDebug.sdlmanager=true);break;
	case xml: (whatDebug.xml=true);break;
	case heap: (whatDebug.heap=true);break;
	case math: (whatDebug.math=true);break;
	case structures: (whatDebug.structures=true);break;
	}


}
void setDebug( int i ){

	switch (i) {
	case error: (whatDebug.error=true);break;
	case bitmap: (whatDebug.bitmap=true);break;
	case font: (whatDebug.font=true);break;
	case gui: (whatDebug.gui=true); break;
	case input: (whatDebug.input=true);break;
	case compilers: (whatDebug.compilers=true);break;
	case configuration: (whatDebug.configuration=true);break;
	case parsers: (whatDebug.parsers=true);break;
	case lists: (whatDebug.lists=true);break;
	case menu: (whatDebug.menu=true);break;
	case network: (whatDebug.network=true);break;
	case resources: (whatDebug.resources=true);break;
	case objects: (whatDebug.objects=true);break;
	case exceptions: (whatDebug.exceptions=true);break;
	case environment: (whatDebug.environment=true);break;
	case filesystem: (whatDebug.filesystem=true);break;
	case sound: (whatDebug.sound=true);break;
	case console: (whatDebug.console=true);break;
	case gameEngine: (whatDebug.gameEngine=true);break;
	case gameStage: (whatDebug.gameStage=true);	break;
	case gameObject: (whatDebug.gameObject=true);break;
	case gameStateMachine: (whatDebug.gameStateMachine=true);break;
	case gameStory: (whatDebug.gameStory=true);break;
	case gameConfiguration: (whatDebug.gameConfiguration=true);break;
	case sdlmanager: (whatDebug.sdlmanager=true);break;
	case xml: (whatDebug.xml=true);break;
	case heap: (whatDebug.heap=true);break;
	case math: (whatDebug.math=true);break;
	case structures: (whatDebug.structures=true);break;

	}


}

void resetDebug( Debug::Modules i ){

	switch (i) {
	case configuration: (whatDebug.configuration=false);break;
	case error: (whatDebug.error=false);break;
	case bitmap: (whatDebug.bitmap=false);break;
	case font: (whatDebug.font=false);break;
	case gui: (whatDebug.gui=false); break;
	case input: (whatDebug.input=false);break;
	case compilers: (whatDebug.compilers=false);break;
	case parsers: (whatDebug.parsers=false);break;
	case lists: (whatDebug.lists=false);break;
	case menu: (whatDebug.menu=false);break;
	case network: (whatDebug.network=false);break;
	case resources: (whatDebug.resources=false);break;
	case objects: (whatDebug.objects=false);break;
	case exceptions: (whatDebug.exceptions=false);break;
	case environment: (whatDebug.environment=false);break;
	case filesystem: (whatDebug.filesystem=false);break;
	case sound: (whatDebug.sound=false);break;
	case console: (whatDebug.console=false);break;
	case gameEngine: (whatDebug.gameEngine=false);break;
	case gameStage: (whatDebug.gameStage=false);break;
	case gameObject: (whatDebug.gameObject=false);break;
	case gameStateMachine: (whatDebug.gameStateMachine=false);break;
	case gameStory: (whatDebug.gameStory=false);break;
	case gameConfiguration: (whatDebug.gameConfiguration=false);break;
	case sdlmanager: (whatDebug.sdlmanager=false);break;
	case xml: (whatDebug.xml=false);break;
	case heap: (whatDebug.heap=false);break;
	case math: (whatDebug.math=false);break;
	case structures: (whatDebug.structures=false);break;

	}

	global_debug_level = i;
}
void resetDebug( int i ){

	switch (i) {
	case configuration: (whatDebug.configuration=false);break;
	case error: (whatDebug.error=false);break;
	case bitmap: (whatDebug.bitmap=false);break;
	case font: (whatDebug.font=false);break;
	case gui: (whatDebug.gui=false); break;
	case input: (whatDebug.input=false);break;
	case compilers: (whatDebug.compilers=false);break;
	case parsers: (whatDebug.parsers=false);break;
	case lists: (whatDebug.lists=false);break;
	case menu: (whatDebug.menu=false);break;
	case network: (whatDebug.network=false);break;
	case resources: (whatDebug.resources=false);break;
	case objects: (whatDebug.objects=false);break;
	case exceptions: (whatDebug.exceptions=false);break;
	case environment: (whatDebug.environment=false);break;
	case filesystem: (whatDebug.filesystem=false);break;
	case sound: (whatDebug.sound=false);break;
	case console: (whatDebug.console=false);break;
	case gameEngine: (whatDebug.gameEngine=false);break;
	case gameStage: (whatDebug.gameStage=false);break;
	case gameObject: (whatDebug.gameObject=false);break;
	case gameStateMachine: (whatDebug.gameStateMachine=false);break;
	case gameStory: (whatDebug.gameStory=false);break;
	case gameConfiguration: (whatDebug.gameConfiguration=false);break;
	case sdlmanager: (whatDebug.sdlmanager=false);break;
	case xml: (whatDebug.xml=false);break;
	case heap: (whatDebug.heap=false);break;
	case math: (whatDebug.math=false);break;
	case structures: (whatDebug.structures=false);break;

	}

	global_debug_level = i;
}


bool getDebug( Debug::Modules i ){

	switch (i) {
	case configuration: return whatDebug.configuration; break;
	case error: return whatDebug.error; break;
	case bitmap: return whatDebug.bitmap; break;
	case font: return whatDebug.font; break;
	case gui: return whatDebug.gui; break;
	case input: return whatDebug.input; break;
	case compilers: return whatDebug.compilers; break;
	case parsers: return whatDebug.parsers; break;
	case lists: return whatDebug.lists; break;
	case menu: return whatDebug.menu; break;
	case network: return whatDebug.network; break;
	case resources: return whatDebug.resources; break;
	case objects: return whatDebug.objects; break;
	case exceptions: return whatDebug.exceptions; break;
	case environment: return whatDebug.environment; break;
	case filesystem: return whatDebug.filesystem; break;
	case sound: return whatDebug.sound; break;
	case console: return whatDebug.console; break;
	case gameEngine: return whatDebug.gameEngine; break;
	case gameStage: return whatDebug.gameStage; break;
	case gameObject: return whatDebug.gameObject; break;
	case gameStateMachine: return whatDebug.gameStateMachine; break;
	case gameStory: return whatDebug.gameStory; break;
	case gameConfiguration: return whatDebug.gameConfiguration; break;
	case sdlmanager: return whatDebug.sdlmanager; break;
	case xml: return whatDebug.xml; break;
	case heap: return whatDebug.heap; break;
	case math: return whatDebug.math; break;
	case structures: return whatDebug.structures; break;

	}
}

} //end namespace debug
} //end SMF
