/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the freeMem Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the freeMem Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#include "SMF_StringUtils.h"

#ifndef _SMF_TEST_LIB_H_
#define _SMF_TEST_LIB_H_
//:
// \file
// \brief Testing software
// \author Tim Cootes
// \verbatim
//  Modifications
//   Apr 2002, Amitha Perera: Copied from vil_test and moved into testlib in an
//                  attempt to consolidate all the repeated test functionality.
//   Sep.2004, Peter Vanroose: added testlib_test_assert_near_relative().
// \endverbatim
namespace SMF {
namespace Util{

//: initialise test counters, check test name 'name' exists
void testlib_test_start(const char* name = 0);
//: increment number of tests, then output msg
void testlib_test_begin(const char* msg);
//: increment success/failure counters
void testlib_test_perform(bool success);
//: output summary of tests performed
int  testlib_test_summary();

//: output msg, then perform test in expr
void testlib_test_assert(const CMyString& msg, bool expr);
//: output msg, then perform test to see if expr is within tol of target
void testlib_test_assert_near(const CMyString& msg, float expr,
                              float target = 0, float tol = 1e-12);
//: output msg, then perform test to see if expr is within tol of target
//void testlib_test_assert_near(const CMyString& msg, vcl_complex<float> expr,
//                              vcl_complex<float> target, float tol = 1e-12);
//: output msg, then test to see if expr is within relative tol of target
void testlib_test_assert_near_relative(const CMyString& msg, float expr,
                                       float target = 0, float tol = 1e-12);
//: output msg, then test to see if expr is within relative tol of target
//void testlib_test_assert_near_relative(const CMyString& msg,
//                                       vcl_complex<float> expr,
//                                       vcl_complex<float> target,
//                                       float tol = 1e-12);
//: output msg, then perform test to see if expr is not within tol of target
void testlib_test_assert_far(const CMyString& msg, float expr,
                             float target = 0, float tol = 1e-12);
//: output msg, then perform test to see if expr is not within tol of target
//void testlib_test_assert_far(const CMyString& msg, vcl_complex<float> expr,
//                             vcl_complex<float> target, float tol = 1e-12);
//: output msg, then perform test to see if expr is equal to target
void testlib_test_assert_equal(const CMyString& msg, long expr, long target);

#define Assert testlib_test_assert
#define AssertNear testlib_test_assert_near
#define AssertFar testlib_test_assert_far

//: initialise test
#define START(s) testlib_test_start(s)

//: TEST function, s is message, test to see if p==v
#define TEST(s,p,v) \
do { \
  testlib_test_begin(s); \
  testlib_test_perform((p)==(v)); \
} while (false)

//: TEST function, s is message, test to see if p==v for integral values
#define TEST_EQUAL(s,p,v) \
do { \
  testlib_test_begin(s); \
  testlib_test_assert_equal("",p,v); \
} while (false)

//: TEST function, s is message, test to see if p is close to v, tolerance t
#define TEST_NEAR(s,p,v,t) \
do { \
  testlib_test_begin(s); \
  testlib_test_assert_near("",p,v,t); \
} while (false)

//: TEST function, message s, test to see if (p-v)/p is small compared to t
#define TEST_NEAR_REL(s,p,v,t) \
do { \
  testlib_test_begin(s); \
  testlib_test_assert_near_relative("",p,v,t); \
} while (false)

//: TEST function, s is message, test to see if p is far from v, tolerance t
#define TEST_FAR(s,p,v,t) \
do { \
  testlib_test_begin(s); \
  testlib_test_assert_far("",p,v,t); \
} while (false)

//: run x, s is message, then test to see if p==v
#define TEST_RUN(s,x,p,v) \
do { \
  testlib_test_begin(s); \
  x; \
  testlib_test_perform((p)==(v)); \
} while (false)

//: Summarise test
#define SUMMARY() return testlib_test_summary()

//: Run a singleton test function
#define RUN_TEST_FUNC(x) \
  testlib_test_start(#x); x(); return testlib_test_summary()

//: Declare the main function.
#define MAIN( testname ) \
  int testname ## _main(int,char*[])

//: Declare the main function with parameter passing.
#define MAIN_ARGS( testname ) \
  int testname ## _main(int argc, char* argv[])

//: A simplified version of the main test, just in one line.
// Avoids compiler warnings about "unused argc and argv".
#define TESTMAIN( testname ) \
  int testname ## _main(int,char*[]) { START(#testname); testname(); SUMMARY(); }

//: A simplified version of the main test, just in one line.
// This (new) variant is to be used with the (new) CMake GENERATE_TEST_DRIVER()
#define TEST_MAIN( testname ) \
  int testname(int,char*[]) { START(#testname); testname(); SUMMARY(); }

//: A simplified version of the main test, with parameter passing.
#define TESTMAIN_ARGS( testname ) \
  int testname ## _main(int argc, char*argv[]) { START(#testname); testname(argc,argv); SUMMARY(); }

//: A simplified version of the main test, with parameter passing.
// This (new) variant is to be used with the (new) CMake GENERATE_TEST_DRIVER()
#define TEST_MAIN_ARGS( testname ) \
  int testname(int argc, char*argv[]) { START(#testname); testname(argc,argv); SUMMARY(); }

//: Another simplified main test.  To be used in a standalone executable.
#undef TESTLIB_DEFINE_MAIN
#define TESTLIB_DEFINE_MAIN(testname) \
  int main() { START(#testname); testname(); return testlib_test_summary(); }

//: A simplified main test with parameter passing.  To be used in a standalone executable.
#undef TESTLIB_DEFINE_MAIN_ARGS
#define TESTLIB_DEFINE_MAIN_ARGS(testname) \
  int main(int argc, char * argv[]) { START(#testname); testname(argc,argv); SUMMARY(); }



void testlib_test (int status, const char *test_description, ...);

void testlib_test_rel (float result, float expected, float relative_error,
              const char *test_description, ...) ;

void testlib_test_abs (float result, float expected, float absolute_error,
              const char *test_description, ...) ;

void testlib_test_factor (float result, float expected, float factor,
                 const char *test_description, ...) ;

void testlib_test_int (int result, int expected, const char *test_description, ...) ;

void testlib_test_str (const char * result, const char * expected, 
              const char *test_description, ...) ;



}//end Util
} //end SMF
#endif // _SMF_TEST_LIB_H_
