/*
 Copyright (C) 2018
     Max-Planck-Institute for Polymer Research

 This file is part of ESPResSo++.

 ESPResSo++ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ESPResSo++ is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ACTIVITY_HPP
#define _ACTIVITY_HPP

#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"

namespace espressopp {
namespace integrator{

class Activity : public Extension
{
public:
   Activity(shared_ptr< System > _system);
   ~Activity ();

   /** Register this class so it can be used from Python. */
   static void registerPython();

   void setValue(int _value);
   int getValue();

   void setNumAP(int _numAP);
   int getNumAP();

   void setVelAP(Real3D _velAP);
   Real3D getVelAP();

   void setHustleProb(real _hustleProb);
   real getHustleProb();

   void setHustleTime(real _hustleTime);
   real getHustleTime();

   void function1();
   void function2();

private:
   int value;

   int numAP;
   Real3D velAP;

   real hustleProb, hustleTime;

   shared_ptr< esutil::RNG > rng;  //!< random number generator

   // SIGNALS
   boost::signals2::connection _befIntP;
   boost::signals2::connection _befIntV;

   void connect();
   void disconnect();

   /** Logger */
   static LOG4ESPP_DECL_LOGGER(theLogger);
};

}
}

#endif // _ACTIVITY_HPP

