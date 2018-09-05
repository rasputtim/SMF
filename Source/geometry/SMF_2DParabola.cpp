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

#include "SMF_Config.h"
#include "geometry/SMF_2DParabola.h"
#include "geometry/SMF_2DLineSegment.h"
#include "geometry/SMF_2DLine.h"
#include "geometry/SMF_2DConics.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{

   CParabola::CParabola(float roc)
      : CConicBase(roc, -1.0)
    {
    }

    float CParabola::sagitta(float dist) const
    {
      return CMath::square(dist) / (2.0 * _roc);
    }

    float CParabola::derivative(float r) const
    {
      return r / _roc;
    }

    bool CParabola::intersect(CVec3D &point, const CRay &ray) const
    {
      const float      ax = ray.pos.x;
      const float      ay = ray.pos.y;
      const float      az = ray.pos.z;
      const float      bx = ray.dir.x;
      const float      by = ray.dir.y;
      const float      bz = ray.dir.z;

      /*
        find intersection point between conical section and line,
        Telescope optics, page 266
      */
      float a = (CMath::square(by) + CMath::square(bx));
      float b = ((by * ay + bx * ax) / _roc - bz) * 2.0;
      float c = (CMath::square(ay) + CMath::square(ax)) / _roc - 2.0 * az;

      float t;

      if (a == 0)
        {
          t = -c / b;
        }
      else
        {
          float d = CMath::square(b) - 4.0 * a * c / _roc;

          if (d < 0)
            return false;               // no intersection

          float s = CMath::sqrt(d);

          if (a * bz < 0)
            s = -s;

          t = (2 * c) / (s - b);
        }

      if (t <= 0)               // ignore intersection if before ray origin
        return false;

      point = ray.pos + ray.dir * t;

      return true;
    }
std::string CParabola2D::toString() const
{
	char str[256];
#if 0
        std::sprintf(str, "PARABOLA> center:(%6.3f, %6.3f)\n"
		"\tMajor Radius:%6.3f \n"
		"\tMinor Radius:%6.3f \n"
		"\tRotation: %6.3f degrees\n",
		center.x, center.y,  a(),b(),180.0*Phi/CMath::PI);
#endif
		return str;
}



namespace _2D{
std::ostream &operator <<(std::ostream &o, const CParabola2D &parabola)
{
	o << parabola.toString();
	return o;
}
}
} //end GEO
}  //end SMF
