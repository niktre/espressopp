/*
  Copyright (C) 2012-2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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

#include "python.hpp"
#include "FastLBInit.hpp"

namespace espressopp {
  namespace integrator {
    /* for abstract class there is no need to define anything here,
     * except of python registration */

    void FastLBInit::registerPython() {
      using namespace espressopp::python;

      class_<FastLBInit, boost::noncopyable>
      ("integrator_FastLBInit", no_init)
        .def("createDenVel", pure_virtual(&FastLBInit::createDenVel))
        .def("setForce", pure_virtual(&FastLBInit::setForce))
        .def("addForce", pure_virtual(&FastLBInit::addForce))
      ;
    }
  }
}
