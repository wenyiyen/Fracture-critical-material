/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.5 $
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel02.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Steel02. 
// This Steel02 is based on an f2c of the FEDEAS material
// Steel02.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//-----------------------------------------------------------------------

#include <math.h>

#include <stdlib.h>
#include "SteelFractureDI.h"
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_SteelFractureDI()
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	int    iData[1];
	double dData[12];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SteelFracture tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData != 11) {
		opserr << "Invalid #args, want: uniaxialMaterial SteelFracture " << iData[0] <<
		" fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? epsPlCr?" << endln;
		return 0;
	} else {
		if (OPS_GetDoubleInput(&numData, dData) != 0) {
			opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] <<
				" fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? epsPlCr?" << endln;
			return 0;
		}
		theMaterial = new SteelFractureDI(iData[0], dData[0], dData[1], dData[2],
			dData[3], dData[4], dData[5], dData[6],
			dData[7], dData[8], dData[9], dData[10]);
	}

	

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SteelFracture Material\n";
		return 0;
	}

	return theMaterial;
}



SteelFractureDI::SteelFractureDI(int tag,
	double _Fy, double _E0, double _b,
	double _R0, double _cR1, double _cR2,
	double _a1, double _a2, double _a3, double _a4, double _epsPlCr) :
	UniaxialMaterial(tag, MAT_TAG_SteelFractureDI),
	//UniaxialMaterial(tag, 0),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4),
	epsPlCr(_epsPlCr)
{
	konP = 0;
	kon = 0;
	eP = E0;
	epsP = 0.0;
	sigP = 0.0;
	sig = 0.0;
	eps = 0.0;
	e = E0;

	epsmaxP = Fy / E0;
	epsminP = -epsmaxP;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	// ************** added for fracture ***************
	epsElN = -Fy / E0;
	epsElP = Fy / E0;
	epsCumPl = 0.0;
	epsCont = 0.0;
	konf = 0;
	konC = 0;

	epsElNP = -Fy / E0;
	epsElPP = Fy / E0;
	epsCumPlP = 0.0;
	epsContP = 0.0;
	konfP = 0;
	konCP = 0;
	// ************** added for fracture ***************
}


SteelFractureDI::SteelFractureDI(void) :
	UniaxialMaterial(0, MAT_TAG_SteelFractureDI)
	// UniaxialMaterial(0, 0)
{
	konP = 0;
}

SteelFractureDI::~SteelFractureDI(void)
{
	// Does nothing
}

UniaxialMaterial*
SteelFractureDI::getCopy(void)
{
	//Steel02 *theCopy = new Steel02(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);
	SteelFractureDI *theCopy = new SteelFractureDI(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, epsPlCr);

	return theCopy;
}

double
SteelFractureDI::getInitialTangent(void)
{
	return E0;
}

int
SteelFractureDI::setTrialStrain(double trialStrain, double strainRate)
{
	double Esh = b * E0;
	double epsy = Fy / E0;
	
	eps = trialStrain;
	double deps = eps - epsP;

	epsmax = epsmaxP;
	epsmin = epsminP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	kon = konP;
	// ************** added for fracture ***************
	epsElN = epsElNP;
	epsElP = epsElPP;
	epsCumPl = epsCumPlP;
	epsCont = epsContP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************

	if (kon == 0 || kon == 3) { // modified C-P. Lamarche 2006


		if (fabs(deps) < 10.0*DBL_EPSILON) {

			e = E0;
			//sig = sigini;                // modified C-P. Lamarche 2006
			sig = 0;
			kon = 3;                     // modified C-P. Lamarche 2006 flag to impose initial stess/strain
			return 0;

		}
		else {

			epsmax = epsy;
			epsmin = -epsy;
			if (deps < 0.0) {
				kon = 2;
				epss0 = epsmin;
				sigs0 = -Fy;
				epspl = epsmin;
			}
			else {
				kon = 1;
				epss0 = epsmax;
				sigs0 = Fy;
				epspl = epsmax;
			}
		}
	}

	// in case of load reversal from negative to positive strain increment, 
	// update the minimum previous strain, store the last load reversal 
	// point and calculate the stress and strain (sigs0 and epss0) at the 
	// new intersection between elastic and strain hardening asymptote 
	// To include isotropic strain hardening shift the strain hardening 
	// asymptote by sigsft before calculating the intersection point 
	// Constants a3 and a4 control this stress shift on the tension side 

	if (kon == 2 && deps > 0.0) {


		kon = 1;
		epsr = epsP;
		sigr = sigP;
		//epsmin = min(epsP, epsmin);
		if (epsP < epsmin)
			epsmin = epsP;
		double d1 = (epsmax - epsmin) / (2.0*(a4 * epsy));
		double shft = 1.0 + a3 * pow(d1, 0.8);
		epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
		sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
		epspl = epsmax;

	}
	else if (kon == 1 && deps < 0.0) {

		// update the maximum previous strain, store the last load reversal 
		// point and calculate the stress and strain (sigs0 and epss0) at the 
		// new intersection between elastic and strain hardening asymptote 
		// To include isotropic strain hardening shift the strain hardening 
		// asymptote by sigsft before calculating the intersection point 
		// Constants a1 and a2 control this stress shift on compression side 

		kon = 2;
		epsr = epsP;
		sigr = sigP;
		//      epsmax = max(epsP, epsmax);
		if (epsP > epsmax)
			epsmax = epsP;

		double d1 = (epsmax - epsmin) / (2.0*(a2 * epsy));
		double shft = 1.0 + a1 * pow(d1, 0.8);
		epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
		sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
		epspl = epsmin;
	}

	

	// calculate current stress sig and tangent modulus E 
	if (kon != 4) { // non-fracture case
		double xi = fabs((epspl - epss0) / epsy);
		double R = R0 * (1.0 - (cR1*xi) / (cR2 + xi));
		double epsrat = (eps - epsr) / (epss0 - epsr);
		double dum1 = 1.0 + pow(fabs(epsrat), R);
		double dum2 = pow(dum1, (1 / R));

		sig = b * epsrat + (1.0 - b)*epsrat / dum2;
		sig = sig * (sigs0 - sigr) + sigr;

		e = b + (1.0 - b) / (dum1*dum2);
		e = e * (sigs0 - sigr) / (epss0 - epsr);

		// calculate damage index (temporary: cumulative plastic strain for "tension")
		if (eps > epsElN && eps < epsElP) {
			// do nothing
		}
		else if (eps > epsElP && deps > 0) { // tensile plastic range
			epsCumPl = epsCumPl + deps / epsy;
		}
		else if (eps < epsElN && deps < 0) { // compressive plastic range
			// could add contribution in closing the voids
		}
		else if (eps > epsElP && deps < 0) { // tensile plastic reverse to elastic
			epsElP = epsP;
			epsElN = epsP - 2 * epsy;
		}
		else if (eps < epsElN && deps > 0) { // compressive plastic reverse to elastic
			epsElP = epsP + 2 * epsy;
			epsElN = epsP;
		}

		// check if fractured
		if (epsCumPl > epsPlCr) {
			kon = 4;
			konf = 1;
			epsCont = epsP - sigP / E0;
			sig = 0;
			e = 0;
		}
	}
	else if (kon == 4) { // fractured case
		if (eps >= epsCont) { // not contacted
			sig = 1e-8;
			e = 1e-8;
			if (deps > 0) {
				konf = 2;
			}
			else {
				konf = 1;
			}
			konC = 0;
		}
		else if (eps < epsCont) { // contacted
			if (!konC) { // at first contact
				konC = 1;
				epsr = epsP;
				sigr = sigP;
				double shft = 1.0; // change here
				/*d1 = (epsmax - epsmin) / (2.0*(a2 * epsy)); // shift compression asymptote
				double shft = 1.0 + a1 * d1^0.8;*/
				epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsCont) / (E0 - Esh);
				sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
				epspl = epsmin;
			}
			else if (konf == 2 && deps > 0) { // compression --> tension
				konf = 1;
				epsr = epsP;
				sigr = sigP;
				if (epsP < epsmin) {
					epsmin = epsP;
				}
				double shft = 1.0; // change here
				/*d1 = (epsmax - epsmin) / (2.0*(a4 * epsy)); // shift tension asymptote
				shft = 1.0 + a3 * d1^0.8;*/
				epss0 = (Fy * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = Fy * shft + Esh * (epss0 - epsy * shft);
				epspl = epsmax;
			}
			else if (konf == 1 && deps < 0) { // tension --> compression
				konf = 2;
				epsr = epsP;
				sigr = sigP;
				if (epsP > epsmax)
					epsmax = epsP;
				double shft = 1.0; // change here
				/*d1 = (epsmax - epsmin) / (2.0*(a2 * epsy)); // shift compression asymptote
				shft = 1.0 + a1 * d1^0.8;*/
				epss0 = (-Fy * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
				sigs0 = -Fy * shft + Esh * (epss0 + epsy * shft);
				epspl = epsmin;
			}
			// compute sig and e
			double xi = fabs((epspl - epss0) / epsy);
			double R = R0 * (1.0 - (0.9*xi) / (0.15 + xi));
			double epsrat = (eps - epsr) / (epss0 - epsr);
			double dum1 = 1.0 + pow(fabs(epsrat), R);;
			double dum2 = pow(dum1, (1 / R));

			sig = b * epsrat + (1.0 - b)*epsrat / dum2;
			sig = sig * (sigs0 - sigr) + sigr;
			// sig = min(sig, 0); // make sure stress is below 0, probably not necessary

			e = b + (1.0 - b) / (dum1*dum2);
			e = e * (sigs0 - sigr) / (epss0 - epsr);

			// update epsCont
			if (epsCont > eps - sig / E0)
				epsCont = eps - sig / E0;
			
		}		
	}

	return 0;
}



double
SteelFractureDI::getStrain(void)
{
	return eps;
}

double
SteelFractureDI::getStress(void)
{
	return sig;
}

double
SteelFractureDI::getTangent(void)
{
	return e;
}

int
SteelFractureDI::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsplP = epspl;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	konP = kon;

	eP = e;
	sigP = sig;
	epsP = eps;

	// ************** added for fracture ***************
	epsElNP = epsElN;
	epsElPP = epsElP;
	epsCumPlP = epsCumPl;
	epsContP = epsCont;
	konfP = konf;
	konCP = konC;
	// ************** added for fracture ***************

	return 0;
}

int
SteelFractureDI::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	kon = konP;

	e = eP;
	sig = sigP;
	eps = epsP;

	// ************** added for fracture ***************
	epsElN = epsElNP;
	epsElP = epsElPP;
	epsCumPl = epsCumPlP;
	epsCont = epsContP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************

	return 0;
}

int
SteelFractureDI::revertToStart(void)
{
	eP = E0;
	epsP = 0.0;
	sigP = 0.0;
	sig = 0.0;
	eps = 0.0;
	e = E0;

	konP = 0;
	epsmaxP = Fy / E0;
	epsminP = -epsmaxP;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	// ************** added for fracture ***************
	epsElN = -Fy / E0;
	epsElP = Fy / E0;
	epsCumPl = 0.0;
	epsCont = 0.0;
	konf = 0;
	konC = 0;

	epsElNP = -Fy / E0;
	epsElPP = Fy / E0;
	epsCumPlP = 0.0;
	epsContP = 0.0;
	konfP = 0;
	konCP = 0;
	// ************** added for fracture ***************

	return 0;
}

int
SteelFractureDI::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(28);
	data(0) = Fy;
	data(1) = E0;
	data(2) = b;
	data(3) = R0;
	data(4) = cR1;
	data(5) = cR2;
	data(6) = a1;
	data(7) = a2;
	data(8) = a3;
	data(9) = a4;
	data(10) = epsminP;
	data(11) = epsmaxP;
	data(12) = epsplP;
	data(13) = epss0P;
	data(14) = sigs0P;
	data(15) = epssrP;
	data(16) = sigsrP;
	data(17) = konP;
	data(18) = epsP;
	data(19) = sigP;
	data(20) = eP;
	data(21) = this->getTag();
	//data(22) = sigini;
	data(22) = epsElNP;
	data(23) = epsElPP;
	data(24) = epsCumPlP;
	data(25) = epsContP;
	data(26) = konfP;
	data(27) = konCP;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelFractureDI::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int
SteelFractureDI::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	static Vector data(28);

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelFractureDI::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	Fy = data(0);
	E0 = data(1);
	b = data(2);
	R0 = data(3);
	cR1 = data(4);
	cR2 = data(5);
	a1 = data(6);
	a2 = data(7);
	a3 = data(8);
	a4 = data(9);
	epsminP = data(10);
	epsmaxP = data(11);
	epsplP = data(12);
	epss0P = data(13);
	sigs0P = data(14);
	epssrP = data(15);
	sigsrP = data(16);
	konP = int(data(17));
	epsP = data(18);
	sigP = data(19);
	eP = data(20);
	this->setTag(int(data(21)));
	//sigini = data(22);
	epsElNP = data(22);
	epsElPP = data(23);
	epsCumPlP = data(24);
	epsContP = data(25);
	konfP = int(data(26));
	konCP = int(data(27));

	e = eP;
	sig = sigP;
	eps = epsP;

	return 0;
}

void
SteelFractureDI::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		//    s << "Steel02:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
		s << "SteelFracture tag: " << this->getTag() << endln;
		s << "  fy: " << Fy << ", ";
		s << "  E0: " << E0 << ", ";
		s << "   b: " << b << ", ";
		s << "  R0: " << R0 << ", ";
		s << " cR1: " << cR1 << ", ";
		s << " cR2: " << cR2 << ", ";
		s << "  a1: " << a1 << ", ";
		s << "  a2: " << a2 << ", ";
		s << "  a3: " << a3 << ", ";
		s << "  a4: " << a4 << ", ";
		s << "  a4: " << epsPlCr << ", ";
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"SteelFracture\", ";
		s << "\"E\": " << E0 << ", ";
		s << "\"fy\": " << Fy << ", ";
		s << "\"b\": " << b << ", ";
		s << "\"R0\": " << R0 << ", ";
		s << "\"cR1\": " << cR1 << ", ";
		s << "\"cR2\": " << cR2 << ", ";
		s << "\"a1\": " << a1 << ", ";
		s << "\"a2\": " << a2 << ", ";
		s << "\"a3\": " << a3 << ", ";
		s << "\"a4\": " << a4 << ", ";
		s << "  a4: " << epsPlCr << ", ";
		//s << "\"sigini\": " << sigini << "}";
	}
}
