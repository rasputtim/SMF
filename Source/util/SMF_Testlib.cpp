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

#include "util/SMF_TestLib.h"
#include "util/SMF_Debug.h"
#include "math/SMF_Math.h"
#include "util/SMF_StringUtils.h"

namespace SMF {
namespace Util{
using namespace MATH;

static int num_test=0;
static int tests_passed=0;
static int tests_failed=0;
static const char* test_name;

void testlib_test_start(const char* name)
{
  num_test = 0;
  tests_passed = 0;
  tests_failed = 0;
  test_name = name;
  Debug::debug(Debug::math,__FUNCTION__) << "-----------------------------------------------------------------------------\n"
           << "Start Testing";
  if (test_name != NULL) Debug::debug(Debug::math,__FUNCTION__) << ' ' << test_name;
  Debug::debug(Debug::math,__FUNCTION__) << ":\n-----------------------------------------------------------------------------\n" << std::endl;
 }


void testlib_test_begin(const char* msg)
{
  num_test++;
  Debug::debug(Debug::math,__FUNCTION__) <<" Test "<< num_test
           <<": " << msg <<" --> "
           << std::endl;
}

// NOTE: We don't pass in the message (see test_begin) because
//       we want to ensure that the message is printed BEFORE
//       the test is executed.  This way when a test crashes
//       we can tell if it was during a test, or between tests.
void testlib_test_perform(bool success)
{
  if (success) {
    tests_passed++;
    Debug::debug(Debug::math,__FUNCTION__) << "  PASSED\n" << std::endl;
  } else {
    tests_failed++;
    Debug::debug(Debug::math,__FUNCTION__) << "**FAILED**\n" << std::endl;
  }
}

static void update (int s)
{
  num_test++;

  if (s == 0) 
    {
      testlib_test_perform(true);
    }
  else
    {
      testlib_test_perform(false);
    }
}

int testlib_test_summary()
{
  Debug::debug(Debug::math,__FUNCTION__) << "-----------------------------------------------------------------------------\n";
 
   if (num_test != tests_passed + tests_failed)
    {
       Debug::debug(Debug::math,__FUNCTION__) << "TEST RESULTS DO NOT ADD UPTotal tests: "<< num_test <<" != Passed: " <<tests_passed<< " + tests failed: " << tests_failed << endl;
       return tests_failed;
    }

  if (test_name) Debug::debug(Debug::math,__FUNCTION__) << test_name << ' ';
  Debug::debug(Debug::math,__FUNCTION__) << "Test Summary: ";
  if (tests_failed > 0)
  {
    if (tests_passed == 0)
      Debug::debug(Debug::math,__FUNCTION__) << "No tests succeeded";
    else if (tests_passed == 1)
      Debug::debug(Debug::math,__FUNCTION__) << "1 test succeeded";
    else
      Debug::debug(Debug::math,__FUNCTION__) << tests_passed <<" tests succeeded";
    if (tests_failed == 1)
      Debug::debug(Debug::math,__FUNCTION__) <<", 1 test tests_failed";
    else
      Debug::debug(Debug::math,__FUNCTION__) <<", "<< tests_failed <<" tests tests_failed";
    Debug::debug(Debug::math,__FUNCTION__)<<"\t\t*****";
  }
  else
  {
    if (tests_passed > 1)
      Debug::debug(Debug::math,__FUNCTION__) << "All "<< tests_passed <<" tests succeeded";
    else if (tests_passed == 1)
      Debug::debug(Debug::math,__FUNCTION__) << "1 test succeeded";
    else
      Debug::debug(Debug::math,__FUNCTION__) << "Test succeeded";
  }
  Debug::debug(Debug::math,__FUNCTION__) << "\n-----------------------------------------------------------------------------\n" << std::endl;
  return tests_failed;
}


void testlib_test_assert(const CMyString & msg, bool expr)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " - " << std::endl;
  testlib_test_perform(expr);
}

void testlib_test_assert_near(const CMyString& msg, float expr, float target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should be " << target << ", is " << expr << ", " << std::endl;
  float diff = CMath::abs(expr - target);
  if (target != 0.0 && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "difference " << diff << ", " << std::endl;
  testlib_test_perform(diff <= tol);
}
#if 0
void testlib_test_assert_near(const CMyString& msg, vcl_complex<float> expr, vcl_complex<float> target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should be " << target << ", is " << expr << ", " << std::endl;
  float diff = CMath::fabs(expr - target);
  if (target != vcl_complex<float>(0,0) && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "difference " << diff << ", " << std::endl;
  testlib_test_perform(diff <= tol);
}
#endif
void testlib_test_assert_near_relative(const CMyString& msg, float expr, float target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should be " << target << ", is " << expr << ", " << std::endl;
  float max = CMath::fabs(target); if (CMath::fabs(expr) > max) max = CMath::fabs(expr);
  if (max==0.0 || target==0.0) max=1.0;
  float diff = CMath::fabs(expr - target) / max;
  if (target != 0.0 && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "relative difference " << diff << ", " << std::endl;
  testlib_test_perform(diff <= tol);
}
#if 0
void testlib_test_assert_near_relative(const CMyString& msg, vcl_complex<float> expr, vcl_complex<float> target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should be " << target << ", is " << expr << ", " << std::endl;
  float max = CMath::fabs(target); if (CMath::fabs(expr) > max) max = CMath::fabs(expr);
  if (max==0.0 || target==vcl_complex<float>(0,0)) max=1.0;
  float diff = CMath::fabs(expr - target) / max;
  if (target != vcl_complex<float>(0,0) && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "relative difference " << diff << ", " << std::endl;
  testlib_test_perform(diff <= tol);
}
#endif
void testlib_test_assert_far(const CMyString& msg, float expr, float target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should not be " << target << ", is " << expr << ", " << std::endl;
  float diff = CMath::fabs(expr - target);
  if (target != 0.0 && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "difference " << diff << ", " << std::endl;
  testlib_test_perform(diff > tol);
}
#if 0
void testlib_test_assert_far(const CMyString& msg, vcl_complex<float> expr, vcl_complex<float> target, float tol)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should not be " << target << ", is " << expr << ", " << std::endl;
  float diff = CMath::fabs(expr - target);
  if (target != vcl_complex<float>(0,0) && diff != 0.0)
    Debug::debug(Debug::math,__FUNCTION__) << "difference " << diff << ", " << std::endl;
  testlib_test_perform(diff > tol);
}
#endif
void testlib_test_assert_equal(const CMyString & msg, long expr, long target)
{
  Debug::debug(Debug::math,__FUNCTION__) << msg << " should be " << target << ", is " << expr << ", " << std::endl;
  long diff = MATH::CMath::fabs(expr - target);
  if (target != 0 && diff != 0)
    Debug::debug(Debug::math,__FUNCTION__) << "difference " << diff << ", " << std::endl;
  testlib_test_perform(diff == 0);
}




void testlib_test (int status, const char *test_description,...)
{
  if (!num_test) testlib_test_start();

  update (status);
    
}


void testlib_test_rel (float result, float expected, float relative_error,
              const char *test_description,...)
{
  int status ;

  if (!num_test) testlib_test_start();

  /* Check for NaN vs inf vs number */

  if (MATH::isNan(result) || MATH::isNan(expected)) 
    {
      status = MATH::isNan(result) != MATH::isNan(expected); 
    }
  else if (MATH::isInf(result) || MATH::isInf(expected)) 
    {
      status = MATH::isInf(result) != MATH::isInf(expected); 
    }
  else if ((expected > 0 && expected < IEEE_FLT_MIN)
           || (expected < 0 && expected > -(IEEE_FLT_MIN)))
    {
      status = -1;
    }
  else if (expected != 0 ) 
    {
      status = (fabs(result-expected)/fabs(expected) > relative_error) ;
    }
  else
    {
      status = (fabs(result) > relative_error) ;
    }

  update (status);

  if (status)
    {

#if HAVE_VPRINTF
      {
        va_list ap;
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status )
        printf(" [%u]", num_test);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
testlib_test_abs (float result, float expected, float absolute_error,
              const char *test_description,...)
{
  int status ;

  if (!num_test) testlib_test_start();

  /* Check for NaN vs inf vs number */

  if (MATH::isNan(result) || MATH::isNan(expected)) 
    {
      status = MATH::isNan(result) != MATH::isNan(expected); 
    }
  else if (MATH::isInf(result) || MATH::isInf(expected)) 
    {
      status = MATH::isInf(result) != MATH::isInf(expected); 
    }
  else if ((expected > 0 && expected < IEEE_FLT_MIN)
           || (expected < 0 && expected > -(IEEE_FLT_MIN)))
    {
      status = -1;
    }
  else 
    {
      status = fabs(result-expected) > absolute_error ;
    }

  update (status);

  if (status )
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status)
        printf(" [%u]", num_test);

      printf ("\n") ;
      fflush (stdout);
    }
}


void testlib_test_factor (float result, float expected, float factor, const char *test_description,...)
{
  int status;

  if (!num_test) testlib_test_start();
  
  if ((expected > 0 && expected < IEEE_FLT_MIN)
      || (expected < 0 && expected > -(IEEE_FLT_MIN)))
    {
      status = -1;
    }
  else if (result == expected) 
    {
      status = 0;
    }
  else if (expected == 0.0) 
    {
      status = (result > expected || result < expected);
    }
  else
    {
      float u = result / expected; 
      status = (u > factor || u < 1.0 / factor) ;
    }

  update (status);

  if (status)
    {
 
#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status)
        printf(" [%u]", num_test);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
testlib_test_int (int result, int expected, const char *test_description,...)
{
  int status = (result != expected) ;

  if (!num_test) testlib_test_start();

  update (status);

  if (status)
    {

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }
      else 
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }

      if (status)
        printf(" [%u]", num_test);

      printf ("\n");
      fflush (stdout);
    }
}

void
testlib_test_str (const char * result, const char * expected, 
              const char *test_description,...)
{
  int status = strcmp(result,expected) ;

  if (!num_test) testlib_test_start();

  update (status);

  if (status)
    {

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status)
        {
          printf(" (%s observed vs %s expected)", result, expected) ;
        }

      if (status)
        printf(" [%u]", num_test);

      printf ("\n");
      fflush (stdout);
    }
}




}//end Util
}//end SMF