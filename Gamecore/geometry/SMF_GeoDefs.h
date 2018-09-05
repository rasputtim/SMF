#ifndef _SMF_GEO_DEF___
#define _SMF_GEO_DEF___

#define	ON_EPSILON					0.1f
#define DEGENERATE_DIST_EPSILON		1e-4f

#define	SIDE_FRONT					0
#define	SIDE_BACK					1
#define	SIDE_ON						2
#define	SIDE_CROSS					3

// plane sides
#define PLANESIDE_FRONT				0
#define PLANESIDE_BACK				1
#define PLANESIDE_ON				2
#define PLANESIDE_CROSS				3

// plane types
#define PLANETYPE_X					0
#define PLANETYPE_Y					1
#define PLANETYPE_Z					2
#define PLANETYPE_NEGX				3
#define PLANETYPE_NEGY				4
#define PLANETYPE_NEGZ				5
#define PLANETYPE_TRUEAXIAL			6	// all types < 6 are true axial planes
#define PLANETYPE_ZEROX				6
#define PLANETYPE_ZEROY				7
#define PLANETYPE_ZEROZ				8
#define PLANETYPE_NONAXIAL			9


enum IntersectResult { NONE_INTERSECTION= 0x00, PARALLEL= 0x01, COINCIDENT=0x02, CONCURRENT =0x04, NOT_INTERESECTING=0x08, INTERESECTING=0x10, INTERESECTING_EXTREMITY_P1=0x20, INTERESECTING_EXTREMITY_P2=0x40 };
namespace SMF{
namespace GEO{
enum TInclusion {
                      eFully,
                      ePartially,
                      eOutside,
                      eUnknown
};

enum eTriangletype {
                      etEquilateral,
                      etIsosceles,
                      etRight,
                      etScalene,
                      etObtuse,
                      etUnknown
};
//**********[ orientation constants ]**********

extern const int RightHandSide;
extern const int LeftHandSide;
extern const int Clockwise;
extern const int CW; //same as Clockwise
extern const int CounterClockwise;
extern const int CCW; //same as CounterClockwise
extern const int CollinearOrientation;
extern const int AboveOrientation;
extern const int BelowOrientation;
extern const int CoplanarOrientation;

} //end GEO
} //end SMF

#endif //_SMF_GEO_DEF___