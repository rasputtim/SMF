
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

	typedef SMF::GEO::CEllipse2D Elipse;
	typedef SMF::GEO::CPoint2D Point;
	typedef SMF::MATH::CVec2D Vec2;
	typedef SMF::GEO::CRay2D Ray;
	typedef SMF::GEO::CLineSegment2D Segment;
	typedef SMF::GEO::CLine2D Line;
	typedef SMF::GEO::CTriangle2D Triangle;
	typedef SMF::GEO::CCircle2D Circle;
	typedef SMF::GEO::CConic2D  Conic;
	SMF::MATH::CSIMD::initHeap();
	//===========================================
	SMF::GEO::CPoint2D::test_point_2d();
	SMF::Debug::debug(SMF::Debug::math,"TESTING")<<"==================   CLINE TESTS ================" <<endl;
	SMF::GEO::CLine2D::test_line_2d();
	SMF::Debug::debug(SMF::Debug::math,"TESTING")<<"==================   CAABBox2D TESTS ================" <<endl;
	SMF::GEO::CAABBox2D::test_box_2d();
	SMF::Debug::debug(SMF::Debug::math,"TESTING")<<"==================   CVec2D TESTS ================" <<endl;
	SMF::MATH::CVec2D::test_vector_2d();

	SMF::GEO::CCircle2D circle(SMF::MATH::CPoint2D(1,2),9);
	SMF::Util::CMyString type = circle.real_type();
	Elipse elipse(SMF::GEO::CPoint2D(1,1),5,3);
	Point p1(4,3);
	Point p2(4,4);
	Point p3(6,1);
	Point p4(7,1);
	bool inside = elipse.contains(p1); //true
	inside = elipse.contains(p2);  //false
	bool inside2 = elipse.contains(p1); //true
	bool inside3 = elipse.contains(p2);//false
	bool inside4 = elipse.edgeContains(p1); //false
	bool inside5 = elipse.edgeContains(p3); //true
	bool inside6 = elipse.edgeContains(p1); //false
	bool inside7 = elipse.edgeContains(p3); //true
	Point f1 = elipse.getFocus1();  //5,1
	Point f2 = elipse.getFocus2();  //-3,1
	float a=(3.141516*90/180);
	double radius1 = elipse.getRadius(a);  //3
	double radius2 = elipse.getRadius(Vec2(cos(a), sin(a))); //3
	double radius3 = elipse.getRadius(Vec2(0,3)); //3

	double dist =elipse.distanceToEdge(p4);  //1
	Point x = elipse.getPoint(a); //1,4
    Point z = elipse.getPoint(Vec2(0,3));  //1,4

	//=======================closestPoint===================

	Line l1(-4,4,-16);
	Line l2(4,4,0);
	Line l3(-2,2,-4);
	Segment seg1(Point(-2,2),Point(2,6));
	Segment seg2(Point(-1,2),Point(2,6));
	Segment seg3(Point(-3,1),Point(2,6));
	Segment seg4(Point(-1,2),Point(5,-4));
	Vec2 vecx(2,2);
	Vec2 vecy(0,2);
	Vec2 vecz(0,-2);
	Vec2 vecw(-2,0);
	vecx.toNormal();
	vecy.toNormal();
	vecz.toNormal();
	vecw.toNormal();
	Ray ray1(Point(0,4),vecx);
	Ray ray2(Point(-1,2),vecy);
	Ray ray3(Point(-1,2),vecz);

	Point dist1,dist2,dist3,dist4,dist5,dist6;
	Point origin(0,0);

	//intersecção line================
	bool int1,int2,int3, int4,int5,int6;

	int1 = l1.intersects(l2);  //true
	int2 = l1.intersects(l3);  //false
	int3 = l1.intersects(seg2);  //true
	int4 = l1.intersects(seg4);  //false
	int5 = l1.intersects(ray2);  //true
	int6 = l1.intersects(ray3);  //false


	dist1 = origin.closestPoint(l1);  //(-2,2)
	dist2 = origin.closestPoint(seg1); //(-2,2)
	dist3 = origin.closestPoint(ray1); // (0,4)
	dist4 = l1.closestPoint(l2);  // (-2,2)
	dist5 = l2.closestPoint(seg2); //(-1.5 , 1.5)
	dist6 = l2.closestPoint(seg3); //-2,2  erro


//-------------------INTERSECSÃO SEGMENTOS------------
	int intresult;
	Segment seg10(Point(2,2),Point(4,4));
	Segment seg20(Point(2,1),Point(4,3));
	Segment seg30(Point(3,3),Point(5,5));
	Segment seg40(Point(2,2),Point(2,8));
	Segment seg50(Point(2,1),Point(2,-5));
	Segment seg60(Point(2,2),Point(2,-5));
	Segment seg70(Point(3,3),Point(10,-2));

	Segment seg80(Point(-6,6),Point(0,0));
	Segment seg90(Point(0,6),Point(-6,0));
	Segment seg100(Point(-2,2),Point(0,0));


	Point intersec;
	intresult= Segment::intersects_m1(seg10,seg20,intersec); //Parallel
	intresult= Segment::intersects_m1(seg10,seg30,intersec); // COincident nenhum ponto
	intresult= Segment::intersects_m1(seg40,seg50,intersec); // NOT_INTERSECTING nos pontos 2,1 e 2,1 // também são coincidentes nenhum ponto
	intresult= Segment::intersects_m1(seg40,seg60,intersec); // COINCIDENT & INTERSECTING , INTERESECTING_EXTREMITY_P1 no ponto 2,2
	intresult= Segment::intersects_m1(seg10,seg70,intersec); // INTERSECTING no ponto 3,3
	intresult= Segment::intersects_m1(seg80,seg90,intersec); // INTERSECTING no ponto -3,3
	intresult= Segment::intersects_m1(seg90,seg100,intersec); // NOT_INTERSECTING no ponto -3,3
	intresult= Segment::intersects_m1(seg80,seg100,intersec); // COINCIDENTES e INTERSECTING, varios pontos


	//using m2
	int intr;
	bool res;
	res= Segment::intersects_m2(seg10,seg20,intersec,&intr); //Parallel
	res= Segment::intersects_m2(seg10,seg30,intersec,&intr); // COincident nenhum ponto
	res= Segment::intersects_m2(seg40,seg50,intersec,&intr); // NOT_INTERSECTING nos pontos 2,1 e 2,1 // também são coincidentes ???
	res= Segment::intersects_m2(seg40,seg60,intersec,&intr); // INTERSECTING E COINCIDENT na EXTREMIDADE ponto 2,2
	res= Segment::intersects_m2(seg10,seg70,intersec,&intr); // INTERSECTING no ponto 3,3
	res= Segment::intersects_m2(seg80,seg90,intersec,&intr); // INTERSECTING no ponto -3,3
	res= Segment::intersects_m2(seg90,seg100,intersec,&intr); // NOT_INTERSECTING no ponto -3,3
	res= Segment::intersects_m2(seg80,seg100,intersec,&intr); // COINCIDENTES e INTERESECTING ...varios pontos

	//================================================
	//TEST TRIANGLE
	Triangle t1(Point(4,2),Point(9,1),Point(6,5));

	//Triangle ray
	Ray raytri1(Point(3,-2),vecx);
	Ray raytri2(Point(6,-2),vecy);
	Ray raytri3(Point(6,3),vecz);
	Ray raytri4(Point(6,3),vecw);
	int count;
	Point Pray1,Pray2;
	t1.intersection(raytri1,&count,&Pray1,&Pray2); //2 -> (6.5,1.5)(7.71,2.71)
	t1.intersection(raytri2,&count,&Pray1,&Pray2); //2 -> (6,1.6) (6,5)
	t1.intersection(raytri3,&count,&Pray1,&Pray2); //1 -> (6,1.6)
	t1.intersection(raytri4,&count,&Pray1,&Pray2); //1 -> (4.66,3)

	//triangle line
	Line linetri1(Point(3,-2),vecx);
	Line linetri2(Point(6,-2),vecy);
	Line linetri3(Point(6,3),vecz);
	Line linetri4(Point(6,3),vecw);

	t1.intersection(linetri1,&count,&Pray1,&Pray2); //2 -> (6.5,1.5)(7.71,2.71)
	t1.intersection(linetri2,&count,&Pray1,&Pray2); //2 -> (6,1.6) (6,5)
	t1.intersection(linetri3,&count,&Pray1,&Pray2); //2 -> (6,1.6) (6,5)
	t1.intersection(linetri4,&count,&Pray1,&Pray2); //2 -> (7.5,3) (4.66,3)

//INCENTER, CIRCUNCENTER, MEDIAN
	Triangle t2(Point(1,2),Point(9,2),Point(7,8));

	Point circuncenter,incenter,orthocenter;
	Segment med1,med2,med3;

	incenter=t2.getIncenter();  // (6.08,4.10)
	circuncenter=t2.getCircumcenter(); //(5,4)
	orthocenter=t2.getOrthocenter(); // (7,3.99)

	med1 = t2.getMedian(0);
	med2 = t2.getMedian(1);
	med3 = t2.getMedian(2);

//CIRC intersection
	Circle cteste1(Point(4,4),3.0f);
	Circle cteste2(Point(10,4),4.0f);
	Circle cteste3(Point(10,4),3.0f);
	Line cirline1(0,4,-4);
	Line cirline2(-2,0,8);
	Line cirline3(-4,4,0);
	Ray  cirray1(Point(5,3),Point(4,1));
	Ray  cirray2(Point(7,1),Point(4,4));
	Segment segcir1(Point(5,3),Point(4,1));
	Point pcircint1,pcircint2;
	int circCount;
	cteste1.intersection(cirline1,&circCount,&pcircint1,&pcircint2);//1 -> (4,1)
	cteste1.intersection(cirline2,&circCount,&pcircint1,&pcircint2);//2 -> (4,7)(4,1)
	cteste1.intersection(cirline3,&circCount,&pcircint1,&pcircint2);//2 -> (6.12,6.12)(1.87,1.87)

	cteste1.intersection(cirray1,&circCount,&pcircint1,&pcircint2);  //1-> (4,1)
	cteste1.intersection(cirray2,&circCount,&pcircint1,&pcircint2);  //2->(1.87,6.12)(6.12,1.87)

	cteste1.intersection(cteste2,&circCount,&pcircint1,&pcircint2);  //2-> (6.41,5,77)(6.41,2.22)
	cteste1.intersection(cteste3,&circCount,&pcircint1,&pcircint2);  //1->(7,4)


	//INTERSEC ELipse reta

	Elipse elipse2(SMF::GEO::CPoint2D(0,0),3,2);
	Line ltelip(-2,1,2);
	Line ltelip2(-2,0,4);
	Point pel1,pel2;
	intresult= Elipse::intersection(elipse2,ltelip,pel1,pel2);
	intresult= Elipse::intersection(elipse2,ltelip2,pel1,pel2);

    SMF::GEO::CConic2D::testEllipseIntersection();

//soluão de sistemas
	float rz1=0,rz2=0,rz3=0,rz4=0;

	//SMF::MATH::CPolynomial::solveQuartic(65304,0,2.31414e7,0,-6.27057e7,rz1,rz2,rz3,rz4);
	//SMF::MATH::CPolynomial::solveQuartic(3600,0,4.93517e6,0,0,rz1,rz2,rz3,rz4);
	SMF::MATH::CPolynomial::solveQuartic(295936,-155976e7,1.82944e8,7.52353e8,6.76594e8,rz1,rz2,rz3,rz4);

	SMF::MATH::CPolyConSystem system;
	//nenhum ponto
	system.setEqu1(16,0,25,-48,0,-64);
	system.setEqu2(9,3,8,0,-80,128);
	system.solveSystem();
	//dois pontos
	system.setEqu1(3,1,4,-6,0,-9);
	system.setEqu2(3,4,4,-18,0,15);
	system.solveSystem();

	system.setEqu1(9,0,10,-198,0,999);
	system.setEqu2(9,0,8,-198,-24,945);
	system.solveSystem();


	system.setEqu1(1,0,1,0,0,-1);
	system.setEqu2(0,0,0,0,1,-1);
	system.solveSystem();

	system.setEqu1(1,0,0,-6,-4,17);
	system.setEqu2(1,2,1,-8,-24,48);
	system.solveSystem();

//Intersecção de Conicas
	Conic c1(16,0,25,-48,0,-64);
	Conic c2(9,3,8,0,-80,128);
	Conic c3(3,1,4,-6,0,-9);
	Conic c4(3,4,4,-18,0,15);
	Conic c5(9,0,10,-198,0,999);
	Conic c6(9,0,8,-198,-24,945);
	Conic c7(1,0,1,0,0,-1);
	Conic c8(0,0,0,0,1,-1);
	Conic c9(1,0,0,-6,-4,17);
	Conic c10(1,2,1,-8,-24,48);
	Point pcon1,pcon2,pcon3,pcon4;
	int points;
	Conic::intersection(c1,c2,points, pcon1,pcon2,pcon3,pcon4);
	Conic::intersection(c3,c4,points, pcon1,pcon2,pcon3,pcon4);
	Conic::intersection(c5,c6,points, pcon1,pcon2,pcon3,pcon4,3);
	Conic::intersection(c7,c8,points, pcon1,pcon2,pcon3,pcon4);
	Conic::intersection(c9,c10,points, pcon1,pcon2,pcon3,pcon4);
	SMF::MATH::CSIMD::shutdown();

    return 0;

}



