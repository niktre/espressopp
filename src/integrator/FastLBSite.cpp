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
    FastLBSite::FastLBSite () {
            f   = std::vector<real>(FastLatticePar::getNumVelsLoc(), 0.);
    }

/*******************************************************************************************/

        /* SET AND GET PART */
    void FastLBSite::setF_i (int _i, real _f) { f[_i] = _f;}
    real FastLBSite::getF_i (int _i) { return f[_i];}

    void FastLBSite::setPhiLoc (int _i, real _phi) { phiLoc[_i] = _phi;}
    real FastLBSite::getPhiLoc (int _i) { return phiLoc[_i];}

/*******************************************************************************************/

    /* HELPFUL OPERATIONS WITH POPULATIONS AND MOMENTS */
    void FastLBSite::scaleF_i (int _i, real _value) { f[_i] *= _value;}

/*******************************************************************************************/

    /* MANAGING STATIC VARIABLES */
    /* create storage for static variables */
    std::vector<real> FastLBSite::phiLoc(19, 0.);

/*******************************************************************************************/

        void FastLBSite::collision(bool _fluct, bool _extForce,
                                                     bool _coupling, Real3D _force, std::vector<real> &_gamma) {
            real m[19];

            calcLocalMoments(m);

            relaxMoments(m, _extForce, _force, _gamma);

            if (_fluct) thermalFluct(m);

            // coupling counts as an external force as well
         if (_extForce) applyForces(m, _force, _gamma);

            btranMomToPop(m);
        }

/*******************************************************************************************/

        /* CALCULATION OF THE LOCAL MOMENTS */
        void FastLBSite::calcLocalMoments (real *m) {
            real f0,
            f1p2, f1m2, f3p4, f3m4, f5p6, f5m6, f7p8, f7m8, f9p10, f9m10,
            f11p12, f11m12, f13p14, f13m14, f15p16, f15m16, f17p18, f17m18;

            /* shorthand functions for "simplified" notation */
            f0     =  getF_i(0);
            f1p2   =  getF_i(1) +  getF_i(2);    f1m2 =  getF_i(1) -  getF_i(2);
            f3p4   =  getF_i(3) +  getF_i(4);    f3m4 =  getF_i(3) -  getF_i(4);
            f5p6   =  getF_i(5) +  getF_i(6);    f5m6 =  getF_i(5) -  getF_i(6);
            f7p8   =  getF_i(7) +  getF_i(8);    f7m8 =  getF_i(7) -  getF_i(8);
            f9p10  =  getF_i(9) + getF_i(10);   f9m10 =  getF_i(9) - getF_i(10);
            f11p12 = getF_i(11) + getF_i(12);  f11m12 = getF_i(11) - getF_i(12);
            f13p14 = getF_i(13) + getF_i(14);  f13m14 = getF_i(13) - getF_i(14);
            f15p16 = getF_i(15) + getF_i(16);  f15m16 = getF_i(15) - getF_i(16);
            f17p18 = getF_i(17) + getF_i(18);  f17m18 = getF_i(17) - getF_i(18);

            /* mass mode */
            m[0] = f0 + f1p2 + f3p4 + f5p6 + f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18;

            /* momentum modes */
            m[1] = f1m2 +   f7m8 +  f9m10 + f11m12 + f13m14;
            m[2] = f3m4 +   f7m8 -  f9m10 + f15m16 + f17m18;
            m[3] = f5m6 + f11m12 - f13m14 + f15m16 - f17m18;

            /* stress modes */
            m[4] = -f0 +   f7p8 + f9p10 + f11p12 + f13p14 + f15p16 + f17p18;
            m[5] = 2.*f1p2 -   f3p4 -  f5p6 +   f7p8 +  f9p10 + f11p12 + f13p14 - 2.* (f15p16 + f17p18);
            m[6] = f3p4 -   f5p6 +  f7p8 +  f9p10 - f11p12 - f13p14;
            m[7] = f7p8 -  f9p10;
            m[8] = f11p12 - f13p14;
            m[9] = f15p16 - f17p18;

            /* kinetic (ghost) modes */
            m[10] = -2.* f1m2 +   f7m8 +  f9m10 + f11m12 + f13m14;
            m[11] = -2.* f3m4 +   f7m8 -  f9m10 + f15m16 + f17m18;
            m[12] = -2.* f5m6 + f11m12 - f13m14 + f15m16 - f17m18;
            m[13] = f7m8 +  f9m10 - f11m12 - f13m14;
            m[14] = -f7m8 +  f9m10 + f15m16 + f17m18;
            m[15] = f11m12 - f13m14 - f15m16 + f17m18;
            m[16] = f0 - 2.* (f1p2 + f3p4 + f5p6) +   f7p8 +  f9p10
            + f11p12 + f13p14 + f15p16 + f17p18;
            m[17] = -2.* f1p2 +   f3p4 +   f5p6 +   f7p8 +  f9p10 + f11p12
            + f13p14 -    2.* (f15p16 + f17p18);
            m[18] = -f3p4 +   f5p6 +   f7p8 +  f9p10 - f11p12 - f13p14;
        }

/*******************************************************************************************/

        /* RELAXATION OF THE MOMENTS TO THEIR EQUILIBRIUM VALUES */
        void FastLBSite::relaxMoments (real *m, bool _extForce, Real3D _f, std::vector<real> &_gamma) {
            // moments on the site //
            real _invTauLoc = 1. / FastLatticePar::getTauLoc();
            Real3D jLoc(m[1], m[2], m[3]);
            jLoc *= FastLatticePar::getALoc();
            jLoc *= _invTauLoc;

            // if we have external forces then modify the eq.fluxes //
            if (_extForce) jLoc += 0.5 * _f; // when doing coupling, the flag is set to 1!

            real _invRhoLoc = 1. / m[0];
            real pi_eq[6];

            pi_eq[0] =  jLoc.sqr()*_invRhoLoc;
            pi_eq[1] =  (jLoc[0]*jLoc[0] - jLoc[1]*jLoc[1])*_invRhoLoc;
            pi_eq[2] =  (3.*jLoc[0]*jLoc[0] - jLoc.sqr())*_invRhoLoc;
            pi_eq[3] =  jLoc[0]*jLoc[1]*_invRhoLoc;
            pi_eq[4] =  jLoc[0]*jLoc[2]*_invRhoLoc;
            pi_eq[5] =  jLoc[1]*jLoc[2]*_invRhoLoc;

            /* relax bulk mode */
            m[4] = pi_eq[0] + _gamma[0] * (m[4] - pi_eq[0]);

            /* relax shear modes */
            m[5] = pi_eq[1] + _gamma[1] * (m[5] - pi_eq[1]);
            m[6] = pi_eq[2] + _gamma[1] * (m[6] - pi_eq[2]);
            m[7] = pi_eq[3] + _gamma[1] * (m[7] - pi_eq[3]);
            m[8] = pi_eq[4] + _gamma[1] * (m[8] - pi_eq[4]);
            m[9] = pi_eq[5] + _gamma[1] * (m[9] - pi_eq[5]);

            /* relax odd modes */
            m[10] *= _gamma[2]; m[11] *= _gamma[2]; 	m[12] *= _gamma[2];
            m[13] *= _gamma[2]; m[14] *= _gamma[2]; 	m[15] *= _gamma[2];

            /* relax even modes */
            m[16] *= _gamma[3]; m[17] *= _gamma[3]; m[18] *= _gamma[3];
        }

/*******************************************************************************************/

        /* ADDING THERMAL FLUCTUATIONS */
        void FastLBSite::thermalFluct (real *m) {
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

        void FastLBSite::applyForces (real *m, Real3D _f, std::vector<real> &_gamma) {
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

        void FastLBSite::btranMomToPop (real *m) {
            int _numVelsLoc = FastLatticePar::getNumVelsLoc();

            // scale modes with inversed coefficients
            for (int i = 0; i < _numVelsLoc; i++) {
                m[i] *= FastLatticePar::getInvBLoc(i);
            }

            setF_i(0,m[0] -m[4] +m[16]);
            setF_i(1,m[0] +m[1] + 2.* (m[5] -m[10] -m[16] -m[17]));
            setF_i(2,m[0] -m[1] + 2.* (m[5] +m[10] -m[16] -m[17]));
            setF_i(3,m[0] +m[2] -m[5] +m[6] - 2.* (m[11] +m[16]) +m[17] -m[18]);
            setF_i(4,m[0] -m[2] -m[5] +m[6] + 2.* (m[11] -m[16]) +m[17] -m[18]);
            setF_i(5,m[0] +m[3] -m[5] -m[6] - 2.* (m[12] +m[16]) +m[17] +m[18]);
            setF_i(6,m[0] -m[3] -m[5] -m[6] + 2.* (m[12] -m[16]) +m[17] +m[18]);

            setF_i(7,m[0] +m[1] +m[2] +m[4] +m[5] +m[6] +m[7] +m[10] +m[11]
                         +m[13] -m[14] +m[16] +m[17] +m[18]);
            setF_i(8,m[0] -m[1] -m[2] +m[4] +m[5] +m[6] +m[7] -m[10] -m[11]
                         -m[13] +m[14] +m[16] +m[17] +m[18]);
            setF_i(9,m[0] +m[1] -m[2] +m[4] +m[5] +m[6] -m[7] +m[10] -m[11]
                         +m[13] +m[14] +m[16] +m[17] +m[18]);
            setF_i(10,m[0] -m[1] +m[2] +m[4] +m[5] +m[6] -m[7] -m[10] +m[11]
                         -m[13] -m[14] +m[16] +m[17] +m[18]);

            setF_i(11,m[0] +m[1] +m[3] +m[4] +m[5] -m[6] +m[8] +m[10] +m[12]
                         -m[13] +m[15] +m[16] +m[17] -m[18]);
            setF_i(12,m[0] -m[1] -m[3] +m[4] +m[5] -m[6] +m[8] -m[10] -m[12]
                         +m[13] -m[15] +m[16] +m[17] -m[18]);
            setF_i(13,m[0] +m[1] -m[3] +m[4] +m[5] -m[6] -m[8] +m[10] -m[12]
                         -m[13] -m[15] +m[16] +m[17] -m[18]);
            setF_i(14,m[0] -m[1] +m[3] +m[4] +m[5] -m[6] -m[8] -m[10] +m[12]
                         +m[13] +m[15] +m[16] +m[17] -m[18]);

            setF_i(15,m[0] +m[2] +m[3] +m[4] - 2.*m[5] +m[9] +m[11] +m[12]
                         +m[14] -m[15] +m[16] - 2.*m[17]);
            setF_i(16,m[0] -m[2] -m[3] +m[4] - 2.*m[5] +m[9] -m[11] -m[12]
                         -m[14] +m[15] +m[16] - 2.*m[17]);
            setF_i(17,m[0] +m[2] -m[3] +m[4] - 2.*m[5] -m[9] +m[11] -m[12]
                         +m[14] +m[15] +m[16] - 2.*m[17]);
            setF_i(18,m[0] -m[2] +m[3] +m[4] - 2.*m[5] -m[9] -m[11] +m[12]
                         -m[14] -m[15] +m[16] - 2.*m[17]);

            /* scale populations with weights */
            for (int i = 0; i < _numVelsLoc; i++) {
                scaleF_i (i, FastLatticePar::getEqWeightLoc(i));
            }
        }

/*******************************************************************************************/

    FastLBSite::~FastLBSite() {
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
