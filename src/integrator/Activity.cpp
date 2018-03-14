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

#include "Activity.hpp"
#include "esutil/RNG.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
   using namespace boost;
   using namespace iterator;

namespace integrator {

LOG4ESPP_LOGGER(Activity::theLogger, "Activity");

Activity::Activity(shared_ptr<System> _system) : Extension(_system)
{

   /* setup random numbers generator */
   if (!_system->rng) {
      throw std::runtime_error("system has no RNG");
   }
   rng = _system->rng;
}

void Activity::disconnect() {
   _befIntP.disconnect();
   _befIntV.disconnect();
}

void Activity::connect() {
   _befIntP = integrator->recalc2.connect ( boost::bind(&Activity::function1, this));
   _befIntV = integrator->befIntV.connect ( boost::bind(&Activity::function2, this));
}

/* Setter and getter*/
void Activity::setValue (int _value) { value = _value;}
int Activity::getValue () { return value;}

void Activity::setNumAP(int _numAP) { numAP = _numAP; }
int Activity::getNumAP() { return numAP; }

void Activity::setVelAP(Real3D _velAP) { velAP = _velAP; }
Real3D Activity::getVelAP() { return velAP; }

void Activity::setHustleProb(real _hustleProb) { hustleProb = _hustleProb; }
real Activity::getHustleProb() { return hustleProb; }

void Activity::setHustleTime(real _hustleTime) { hustleTime = _hustleTime; }
real Activity::getHustleTime() { return hustleTime; }

void Activity::function1() {
   printf("befIntP: Hello Activity");
}

void Activity::function2() {
   printf("befIntV: Hello Activity");
}

/* Destructor of the Activity */
Activity::~Activity() {
   disconnect();
}

/******************************
 ** REGISTRATION WITH PYTHON **
 ******************************/

void Activity::registerPython() {

   using namespace espressopp::python;
   class_<Activity, shared_ptr<Activity>, bases<Extension> >

   ("integrator_Activity", init<	shared_ptr< System >>())
   .add_property("numAP", &Activity::getNumAP, &Activity::setNumAP)
   .add_property("velAP", &Activity::getVelAP, &Activity::setVelAP)
   .add_property("hustleProb", &Activity::getHustleProb, &Activity::setHustleProb)
   .add_property("hustleTime", &Activity::getHustleTime, &Activity::setHustleTime)
   .def("connect", &Activity::connect)
   .def("disconnect", &Activity::disconnect)
   ;
}

}
}
