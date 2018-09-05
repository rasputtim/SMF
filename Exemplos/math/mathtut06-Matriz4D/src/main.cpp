
#include <io.h>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include "SMF.h"

using namespace std;
using std::set_terminate;
using std::set_unexpected;
using std::set_new_handler;

#undef printf




//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{

		//==============================
	SMF::Debug::setDebugAll();
	SMF::Debug::setFilename("debug_math.txt");
	SMF::Debug::setDebugMethod(SMF::Debug::File);


	SMF::MATH::CSIMD::initHeap();
	//===========================================
	SMF::MATH::CSIMD::TestMat4D_f( );


	//================================================
	SMF::MATH::CSIMD::shutdown();

    return 0;

}



