
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





//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{

		//==============================
	SMF::Debug::setDebugAll();
	SMF::Debug::setFilename("debug_math.txt");
	SMF::Debug::setDebugMethod(SMF::Debug::File);

    const SMF::MATH::CVec3D vec(1.0f,2.0f,3.0f);

	SMF::MATH::CSIMD::initHeap();
	//===========================================
	SMF::CCMDLineArgs args;
	SMF::MATH::CSIMD::Test3D_f( args );
	//printf("Antes do Teste");
	//float result = SMF::MATH::CSIMD::Vector3D_Length(&vec);
	SMF::MATH::CSIMD::TesteLenght3D( args );


	//================================================
	SMF::MATH::CSIMD::shutdown();

    return 0;

}



