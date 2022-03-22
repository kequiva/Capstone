/*******************************************************************************
Header file for a cosmology library of general use in observational astronomy
Copyright (C) 2003-2021  Joshua Kempner

Version 2.1.5

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Please send any bug fixes, enhancements, or useful comments by email to
josh@kempner.net.
*******************************************************************************/

#ifndef __COSMO_H__
#define __COSMO_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Class to implement the cosmology
////////////////////////////////////////////////////////////////////////////////
class Cosmo
{
    // cosmological parameters
    double H0_, q0_;      // Hubble constant at z=0, q_0
    double OmegaM_, OmegaL_, Omegak_; // scaled densities of matter, vacuum energy,
			// and curvature
    double dH_;           // Hubble distance
    double z_;            // redshift of source
    // distance measures to source in Mpc
    double dA_;		// Angular diameter distance
    double dL_;		// Luminosity distance
    double dC_;		// Comoving line-of-sight distance
    double dM_;		// Comoving transverse distance
    double VC_;		// Comoving volume out to redshift z_ in cubic Gpc
    // other quantities derived from z_
    double tL_;		// Lookback time to z_ in seconds
    double age_;		// Current age of the Universe in seconds
    double scale_;        // kpc/" at redshift of source
    double rhoCrit_;	// Critical density at redshift of source

    // private member functions
    void init(const double, const double, const double);// NOT exclusive to constructors
    inline void clone(const Cosmo& a) // used in copy constructor and in assignment
    {
        init(a.H0_, a.OmegaM_, a.OmegaL_);
        setRedshift(a.z_);
    }
    inline double E(const double z) // calculate expansion factor at a given redshift
    {
        return sqrt(OmegaM_ * CUBE(1 + z) + Omegak_ * SQR(1 + z) + OmegaL_);
    }
    inline double inverseOfE(const double z) { return 1.0 / E(z); }
	inline double lookbackIntegrand(const double z) { return 1.0 / (1 + z) / E(z); }
    double ageIntegrand(const double z);
    typedef double (Cosmo::*PFD)(const double);
    double romberg(PFD, double, double);
    void setDistances(); // set the distance measures
    inline double SQR(const double a) { return a*a; }
    inline double CUBE(const double a) { return a*a*a; }

public:
    // constructors etc.
    Cosmo();
    Cosmo(const double, const double, const double);
    Cosmo(const Cosmo&);
    Cosmo& operator=(const Cosmo&);

    // inspection functions
    inline double z() { return z_; }    // redshift of source
    inline double dL() { return dL_; }  // Luminosity distance (Mpc)
    inline double dA() { return dA_; }  // Angular diameter distance (Mpc)
    inline double dC() { return dC_; }  // Comoving line-of-sight distance (Mpc)
    inline double dM() { return dM_; }  // Comoving transverse distance (Mpc)
    inline double VC() { return VC_; }  // comoving volume (Gpc**3)
    inline double lookback() { return tL_; }  // lookback time to source (sec)
    inline double scale() { return scale_; }  // kpc/" at redshift of source
    inline double rhoCrit() { return rhoCrit_; }  // critial density at source
    inline double age() { return age_; }	// Current age of the Universe (sec)
    void printParams(ostream&, const char*); // print cosmological parameters
    void printParamsAsHtml(ostream&, const char*); // same as printParams but formatted in HTML
    void printLong();	      // print derived quantities to STDOUT
    void printAsHtml();     // equivalent to printLong but formatted in HTML
    void printShortHeader(ostream&);  // print header line for columns in printShort()
    void printShort(ostream&);  // print distances in columns

    // mutation functions
    void setCosmology(const double, const double, const double);
    void setRedshift(const double);
    void getCosmologyFromUser();
};

// non-member functions
int isNumeric(const string&);
double promptForParam(const char*, const double);

#endif // __COSMO_H__
