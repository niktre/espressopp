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
#include "FastLBSite.hpp"
#include <iomanip>

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  using namespace iterator;
  namespace integrator {

/*******************************************************************************************/

        /* SET AND GET PART */
    void FastLBMom::setPhiLoc (int _i, real _phi) { phiLoc[_i] = _phi;}
    real FastLBMom::getPhiLoc (int _i) { return phiLoc[_i];}

/*******************************************************************************************/

    /* MANAGING STATIC VARIABLES */
    /* create storage for static variables */
    std::vector<real> FastLBMom::phiLoc(19, 0.);

/*******************************************************************************************/

        void FastLBMom::collision(bool _fluct, bool _extForce, 
                                  bool _coupling, Real3D _force, 
                                  std::vector<real> &_gamma) {
            real m[19];

//            calcLocalMoments(m);

            relaxMoments(m, _extForce, _force, _gamma);

//            if (_fluct) thermalFluct(m);

            // coupling counts as an external force as well
//         if (_extForce) applyForces(m, _force, _gamma);

            btranMomToPop(m);
        }

/*******************************************************************************************/

        /* RELAXATION OF THE MOMENTS TO THEIR EQUILIBRIUM VALUES */
        void FastLBMom::relaxMoments (real *m, bool _extForce, Real3D _f, std::vector<real> &_gamma) {
            // moments on the site //
            m[0] = getMom_i(0);
            m[1] = getMom_i(1);
            m[2] = getMom_i(2);
            m[3] = getMom_i(3);
            
            real _invTauLoc = 1. / FastLatticePar::getTauLoc();
            Real3D jLoc(m[1], m[2], m[3]);
            jLoc *= FastLatticePar::getALoc();
            jLoc *= _invTauLoc;

            // if we have external forces then modify the eq.fluxes //
            if (_extForce) jLoc += 0.5 * _f; // when doing coupling, the flag is set to 1!

            real _invRhoLoc = 1. / m[0];
            real pi_eq[6];

            /* relax bulk mode */
            m[4] =  jLoc.sqr()*_invRhoLoc;
            
            /* relax shear modes */
            m[5] =  (jLoc[0]*jLoc[0] - jLoc[1]*jLoc[1])*_invRhoLoc;
            m[6] =  (3.*jLoc[0]*jLoc[0] - jLoc.sqr())*_invRhoLoc;
            m[7] =  jLoc[0]*jLoc[1]*_invRhoLoc;
            m[8] =  jLoc[0]*jLoc[2]*_invRhoLoc;
            m[9] =  jLoc[1]*jLoc[2]*_invRhoLoc;
            
            m[10] = 0.; m[11] = 0.; m[12] = 0.; 
            m[13] = 0.; m[14] = 0.; m[15] = 0.; 
            m[16] = 0.; m[17] = 0.; m[18] = 0.; 
        }

/*******************************************************************************************/

        /* ADDING THERMAL FLUCTUATIONS */
        void FastLBMom::thermalFluct (real *m) {
            /* values of PhiLoc were already set in LatticeBoltzmann.cpp */
            int _numVelsLoc = FastLatticePar::getNumVelsLoc();
            real rootRhoLoc = sqrt(12.*m[0]); // factor 12. comes from usage of
            // not gaussian but uniformly distributed random numbers

//			real rootRhoLoc = sqrt(m[0]); // Gaussian version

            for (int l = 4; l < _numVelsLoc; l++) {
                m[l] += rootRhoLoc*getPhiLoc(l)*((*FastLatticePar::rng)() - 0.5);
//				m[l] += rootRhoLoc*getPhiLoc(l)*((FastLatticePar::rng)->normal()); //Gaussian
            }
        }

/*******************************************************************************************/

        void FastLBMom::applyForces (real *m, Real3D _f, std::vector<real> &_gamma) {
            // set velocity _u
            Real3D _u = 0.5 * _f;
            _u[0] += m[1]; 	_u[1] += m[2]; _u[2] += m[3];
            _u /= m[0];

            /* update momentum modes */
            m[1] += _f[0];
            m[2] += _f[1];
            m[3] += _f[2];

            /* update stress modes */
            // See def. of _sigma (Eq.198) in B.Dünweg & A.J.C.Ladd in Adv.Poly.Sci. 221, 89-166 (2009)
            real _sigma[6];
            real _gamma_sp = _gamma[1] + 1.;
            real _gamma_sph = 0.5 * _gamma_sp;

            real _scalp = _u*_f;
            real _secTerm = (1./3.)*(_gamma[0] - _gamma[1])*_scalp;

            _sigma[0] = _gamma_sp*_u[0]*_f[0] + _secTerm;
            _sigma[1] = _gamma_sp*_u[1]*_f[1] + _secTerm;
            _sigma[2] = _gamma_sp*_u[2]*_f[2] + _secTerm;
            _sigma[3] = _gamma_sph*(_u[0]*_f[1]+_u[1]*_f[0]);
            _sigma[4] = _gamma_sph*(_u[0]*_f[2]+_u[2]*_f[0]);
            _sigma[5] = _gamma_sph*(_u[1]*_f[2]+_u[2]*_f[1]);

            m[4] += _sigma[0]+_sigma[1]+_sigma[2];
            m[5] += 2.*_sigma[0]-_sigma[1]-_sigma[2];
            m[6] += _sigma[1]-_sigma[2];
            m[7] += _sigma[3];
            m[8] += _sigma[4];
            m[9] += _sigma[5];
        }

/*******************************************************************************************/

        void FastLBMom::btranMomToPop (real *m) {
            int _numVelsLoc = FastLatticePar::getNumVelsLoc();

            // scale modes with inversed coefficients
            for (int i = 0; i < _numVelsLoc; i++) {
                m[i] *= FastLatticePar::getInvBLoc(i);
            }
            
            real f[19];
            
            f[0] = m[0] -m[4] +m[16];

            f[1] = m[0] +m[1] + 2.* (m[5] -m[10] -m[16] -m[17]);
            
            f[2] = m[0] -m[1] + 2.* (m[5] +m[10] -m[16] -m[17]);
            
            f[3] = m[0] +m[2] -m[5] +m[6] - 2.* (m[11] +m[16]) +m[17] -m[18];
            f[4] = m[0] -m[2] -m[5] +m[6] + 2.* (m[11] -m[16]) +m[17] -m[18];
            f[5] = m[0] +m[3] -m[5] -m[6] - 2.* (m[12] +m[16]) +m[17] +m[18];
            f[6] = m[0] -m[3] -m[5] -m[6] + 2.* (m[12] -m[16]) +m[17] +m[18];

            f[7] = m[0] +m[1] +m[2] +m[4] +m[5] +m[6] +m[7] +m[10] +m[11]
                         +m[13] -m[14] +m[16] +m[17] +m[18];
            f[8] = m[0] -m[1] -m[2] +m[4] +m[5] +m[6] +m[7] -m[10] -m[11]
                         -m[13] +m[14] +m[16] +m[17] +m[18];
            f[9] = m[0] +m[1] -m[2] +m[4] +m[5] +m[6] -m[7] +m[10] -m[11]
                         +m[13] +m[14] +m[16] +m[17] +m[18];
            f[10] = m[0] -m[1] +m[2] +m[4] +m[5] +m[6] -m[7] -m[10] +m[11]
                         -m[13] -m[14] +m[16] +m[17] +m[18];

            f[11] = m[0] +m[1] +m[3] +m[4] +m[5] -m[6] +m[8] +m[10] +m[12]
                         -m[13] +m[15] +m[16] +m[17] -m[18];
            f[12] = m[0] -m[1] -m[3] +m[4] +m[5] -m[6] +m[8] -m[10] -m[12]
                         +m[13] -m[15] +m[16] +m[17] -m[18];
            f[13] = m[0] +m[1] -m[3] +m[4] +m[5] -m[6] -m[8] +m[10] -m[12]
                         -m[13] -m[15] +m[16] +m[17] -m[18];
            f[14] = m[0] -m[1] +m[3] +m[4] +m[5] -m[6] -m[8] -m[10] +m[12]
                         +m[13] +m[15] +m[16] +m[17] -m[18];

            f[15] = m[0] +m[2] +m[3] +m[4] - 2.*m[5] +m[9] +m[11] +m[12]
                         +m[14] -m[15] +m[16] - 2.*m[17];
            f[16] = m[0] -m[2] -m[3] +m[4] - 2.*m[5] +m[9] -m[11] -m[12]
                         -m[14] +m[15] +m[16] - 2.*m[17];
            f[17] = m[0] +m[2] -m[3] +m[4] - 2.*m[5] -m[9] +m[11] -m[12]
                         +m[14] +m[15] +m[16] - 2.*m[17];
            f[18] = m[0] -m[2] +m[3] +m[4] - 2.*m[5] -m[9] -m[11] +m[12]
                         -m[14] -m[15] +m[16] - 2.*m[17];

            /* scale populations with weights */
            for (int i = 0; i < _numVelsLoc; i++) {
                f[i] *= FastLatticePar::getEqWeightLoc(i);
            }
            
            
        }

/*******************************************************************************************/

    FastLBMom::FastLBMom () {
      mom = std::vector<real>(4, 0.);
    }

        /* SET AND GET PART */
    void FastLBMom::setMom_i (int _i, real _mom) { mom[_i] = _mom;}
    real FastLBMom::getMom_i (int _i) { return mom[_i];}

    FastLBMom::~FastLBMom() {
    }

/*******************************************************************************************/

        FastLatticePar::FastLatticePar (shared_ptr<System> system, int _numVelsLoc, real _a, real _tau) {
            setNumVelsLoc(_numVelsLoc);
            setALoc(_a);
            setTauLoc(_tau);

            initEqWeights();
            initInvBLoc();

            /* declare RNG */
            if (!system->rng) {
                throw std::runtime_error("system has no RNG");
            }
            rng = system->rng;
        }

/*******************************************************************************************/

        void FastLatticePar::setNumVelsLoc (int _numVelsLoc) {numVelsLoc = _numVelsLoc;}
        int FastLatticePar::getNumVelsLoc () {return numVelsLoc;}
        void FastLatticePar::setALoc (real _a) {aLoc = _a;}
        real FastLatticePar::getALoc () {return aLoc;}
        void FastLatticePar::setTauLoc (real _tau) {tauLoc = _tau;}
        real FastLatticePar::getTauLoc () {return tauLoc;}

        void FastLatticePar::setInvBLoc (int _i, real _b) { inv_bLoc[_i] = _b;}
        real FastLatticePar::getInvBLoc (int _i) { return inv_bLoc[_i];}

        void FastLatticePar::setEqWeightLoc (int _i, real _w) { eqWeightLoc[_i] = _w;}
        real FastLatticePar::getEqWeightLoc (int _i) { return eqWeightLoc[_i];}

/*******************************************************************************************/

        void FastLatticePar::initEqWeights () {
            real _sh1 = 1./18.;
            real _sh2 = 1./36.;

            // default D3Q19 model
            setEqWeightLoc(0, 1./3.);
            setEqWeightLoc(1, _sh1);  setEqWeightLoc(2, _sh1);  setEqWeightLoc(3, _sh1);
            setEqWeightLoc(4, _sh1);  setEqWeightLoc(5, _sh1);  setEqWeightLoc(6, _sh1);
            setEqWeightLoc(7, _sh2);  setEqWeightLoc(8, _sh2);  setEqWeightLoc(9, _sh2);
            setEqWeightLoc(10, _sh2); setEqWeightLoc(11, _sh2); setEqWeightLoc(12, _sh2);
            setEqWeightLoc(13, _sh2); setEqWeightLoc(14, _sh2); setEqWeightLoc(15, _sh2);
            setEqWeightLoc(16, _sh2); setEqWeightLoc(17, _sh2); setEqWeightLoc(18, _sh2);
        }

/*******************************************************************************************/

        /* INITIALISATION OF THE LATTICE MODEL: EQ.WEIGHTS, INVERSED COEFFICIENTS */
        void FastLatticePar::initInvBLoc () {
            setInvBLoc(0, 1.);
            setInvBLoc(1, 3.);				setInvBLoc(2, 3.);				setInvBLoc(3, 3.);
            setInvBLoc(4, 3./2.);			setInvBLoc(5, 3./4.);			setInvBLoc(6, 9./4.);
            setInvBLoc(7, 9.);				setInvBLoc(8, 9.);				setInvBLoc(9, 9.);
            setInvBLoc(10, 3./2.);		setInvBLoc(11, 3./2.);		setInvBLoc(12, 3./2.);
            setInvBLoc(13, 9./2.);		setInvBLoc(14, 9./2.);		setInvBLoc(15, 9./2.);
            // In PRE 76, 036704 (2007) the 17th and 18th Bi are swapped in comp. to Ulf's thesis. Here we use the latter one
            setInvBLoc(16, 1./2.);		setInvBLoc(17, 3./4.);		setInvBLoc(18, 9./4.);
        }

/*******************************************************************************************/

        shared_ptr< esutil::RNG > FastLatticePar::rng;		// initializer
        int FastLatticePar::numVelsLoc = 0;								// initializer
        real FastLatticePar::aLoc = 0.;										// initializer
        real FastLatticePar::tauLoc = 0.;									// initializer
        std::vector<real> FastLatticePar::eqWeightLoc(19, 0.);
        std::vector<real> FastLatticePar::inv_bLoc(19, 0.);

/*******************************************************************************************/

        FastLatticePar::~FastLatticePar() {
        }

/*******************************************************************************************/
        FastLBForce::FastLBForce () {
            extForceLoc = Real3D(0.);
            couplForceLoc = Real3D(0.);
        }

        void FastLBForce::setExtForceLoc (Real3D _extForceLoc) {extForceLoc = _extForceLoc;}
        Real3D FastLBForce::getExtForceLoc () {return extForceLoc;}

        void FastLBForce::setCouplForceLoc (Real3D _couplForceLoc) {couplForceLoc = _couplForceLoc;}
        Real3D FastLBForce::getCouplForceLoc () {return couplForceLoc;}

        void FastLBForce::addExtForceLoc (Real3D _extForceLoc) {extForceLoc += _extForceLoc;}
        void FastLBForce::addCouplForceLoc (Real3D _couplForceLoc) {couplForceLoc += _couplForceLoc;}

        FastLBForce::~FastLBForce() {
        }
    }
}
