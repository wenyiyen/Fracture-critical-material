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
#include <MaterialResponse.h> // for recorder
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
	double dData[13];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SteelFracture tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData != 12) {
		opserr << "Invalid #args, want: uniaxialMaterial SteelFracture " << iData[0] <<
			" fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m?" << endln;
		return 0;
	}
	else {
		if (OPS_GetDoubleInput(&numData, dData) != 0) {
			opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] <<
				" fy? E? b? R0? cR1? cR2? a1? a2? a3? a4? sigcr? m?" << endln;
			return 0;
		}
		theMaterial = new SteelFractureDI(iData[0], dData[0], dData[1], dData[2],
			dData[3], dData[4], dData[5], dData[6],
			dData[7], dData[8], dData[9], dData[10], dData[11]);
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
	double _a1, double _a2, double _a3, double _a4, double _sigcr, double _m) :
	UniaxialMaterial(tag, MAT_TAG_SteelFractureDI),
	//UniaxialMaterial(tag, 0),
	Fy(_Fy), E0(_E0), b(_b), R0(_R0), cR1(_cR1), cR2(_cR2), a1(_a1), a2(_a2), a3(_a3), a4(_a4),
	sigcr(_sigcr), m(_m)
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
	epsCont = 0.0;
	konf = 0;
	konC = 0;

	epsContP = 0.0;
	konfP = 0;
	konCP = 0;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = 0.0;
	DI_cacheP = 0.0;
	isStartP = 1;
	isSecondP = 1;
	sigPDIP = 0.0;
	slopePP = 0.0;
	rev_cntP = {};
	rev_cntP.resize(10);
	rev_cnt_endIdxP = -1;
	isReadRevP = 0;

	DI = 0.0;
	DI_cache = 0.0;
	isStart = 1;
	isSecond = 1;
	sigPDI = 0.0;
	slopeP = 0.0;
	rev_cnt = {};
	rev_cnt.resize(10);
	rev_cnt_endIdx = -1;
	isReadRev = 0;
	// ************** added for DI ********************

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
	SteelFractureDI *theCopy = new SteelFractureDI(this->getTag(), Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigcr, m);

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
	epsCont = epsContP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DI = DIP;
	DI_cache = DI_cacheP;
	isStart = isStartP;
	isSecond = isSecondP;
	sigPDI = sigPDIP;
	slopeP = slopePP;
	rev_cnt = rev_cntP;
	rev_cnt_endIdx = rev_cnt_endIdxP;
	isReadRev = isReadRevP;
	// ************** added for DI ********************

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

		// calculate DI
		// DI = DI + 0.0038;
		calcDI(sigcr, m, isStart, isSecond, sig, sigPDI, DI_cache, DI, slopeP, rev_cnt, rev_cnt_endIdx, isReadRev);

		// check if fractured
		if (DI > 1) {
			kon = 4;
			konf = 1;
			epsCont = epsP - sigP / E0;
			konC = 1;
			sig = 0.0;
			// stiffness? check here			
			e = 0.0;
		}
		else {
			DI = DI_cache;
		}

	}
	else if (kon == 4) { // fractured case
		if (eps >= epsCont) { // not contacted
			sig = 0.0;
			e = 0.0;
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
				konf = 1; // konf == 1 means tensile
				if (epsP < epsmin) {
					epsmin = epsP;
				}
				sig = sigP + deps * E0;
				e = E0;
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
			}
			else if (konf == 1) {
				sig = sigP + deps * E0;
				e = E0;
			} 
			else if (konf == 2) {
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
			}
						
			
			if (sig > 0) {
				epsCont = epsP;
				sig = 0.0;
				e = 0.0;
			}
			
		}
	}

	return 0;
}

void
SteelFractureDI::calcDI(double sigcr, double m, int &isStart, int &isSecond, double sig, double &sigPDI, double &DI_cache, double &DI, double &slopeP, std::vector<double> &rev_cnt, int &rev_cnt_endIdx, int &isReadRev)
{
	// initialize variables ****
	double sig_tmp;
	double slope;
	std::vector<double> rev_tmp(2, 0.0);
	std::vector<double> rev_cnt_cache;
	int rev_tmp_endIdx = -1;
	int rev_cnt_endIdx_cache = -1;
	int isReadRev_cache = 0;
	
	if (DI > 1) {
		return;
	}

	if (sig < 0) {
		sig_tmp = 0.0;
	}
	else {
		sig_tmp = sig;
	}

	// if is starting point
	if (isStart) {
		isStart = 0;
		sigPDI = sig_tmp;
		return;
	}

	// determine slope sign
	if (sig_tmp == sigPDI) {
		slope = slopeP;
	}
	else {
		slope = sig_tmp - sigPDI;
	}

	// if there is a reversal, include into the reversal vector
	if (isSecond) {
		rev_cnt[0] = 0.0;
		rev_cnt_endIdx = 0;
		rev_tmp[0] = sig_tmp;
		rev_tmp[1] = 0.0;
		rev_tmp_endIdx = 0;
		isReadRev = 1;
		isSecond = 0;
	}
	else if (returnSign(slope) != returnSign(slopeP) && returnSign(slopeP) != 0) {
		rev_tmp[0] = sigPDI;
		rev_tmp[1] = sig_tmp;
		rev_tmp_endIdx = 1;
	}
	else {
		rev_tmp[0] = sig_tmp;
		rev_tmp[1] = 0.0;
		rev_tmp_endIdx = 0;
	}
	
	// count cycles and calculate damage
	while (rev_cnt_endIdx != -1) {
		if (isReadRev) {
			if (rev_tmp_endIdx == -1) {
				break;
			}
			else if (rev_tmp_endIdx == 0) {
				// cache
				rev_cnt_cache = rev_cnt;
				rev_cnt_endIdx_cache = rev_cnt_endIdx;
				isReadRev_cache = isReadRev;
				DI_cache = DI;
				// read
				rev_cnt_endIdx = rev_cnt_endIdx + 1;
				if (rev_cnt_endIdx > rev_cnt.size()) {
					opserr << "Vector 'rev_cnt' overflow, check material source code.\n";
				}
				rev_cnt[rev_cnt_endIdx] = rev_tmp[0];
				rev_tmp[0] = rev_tmp[1];
				rev_tmp_endIdx = rev_tmp_endIdx - 1;
			}
			else {
				// read
				rev_cnt_endIdx = rev_cnt_endIdx + 1;
				rev_cnt[rev_cnt_endIdx] = rev_tmp[0];
				rev_tmp[0] = rev_tmp[1];
				rev_tmp_endIdx = rev_tmp_endIdx - 1;
			}
		}
		// < 3 reversals?
		if (rev_cnt_endIdx < 3 - 1) {
			isReadRev = 1;
			continue;
		}
		// form X and Y using the 3 latest nondiscarded reversals
		double X = fabs(rev_cnt[rev_cnt_endIdx] - rev_cnt[rev_cnt_endIdx - 1]);
		double Y = rev_cnt[rev_cnt_endIdx - 1] - rev_cnt[rev_cnt_endIdx - 2];
		int Y_sign = returnSign(Y);
		Y = fabs(Y);
		// r(X) < r(Y)?
		if (X < Y) {
			isReadRev = 1;
			continue;
		}
		// Y includes Z?
		if (rev_cnt_endIdx == 3 - 1) {
			if (Y_sign > 0) {
				// count Y as 0.5 cycle
				DI = DI + 0.5 / 0.5 * pow(Y / sigcr, 1 / m);
			}
			// discard the first reversal of Y
			rev_cnt[rev_cnt_endIdx - 2] = rev_cnt[rev_cnt_endIdx - 1];
			rev_cnt[rev_cnt_endIdx - 1] = rev_cnt[rev_cnt_endIdx];
			rev_cnt_endIdx = rev_cnt_endIdx - 1;
			isReadRev = 0;
		}
		else {
			// count Y as 1 cycle
			DI = DI + 0.5 / 0.5 * pow(Y / sigcr, 1 / m); // here adjusted for "tensile"
														 // discard both points of Y
			rev_cnt[rev_cnt_endIdx - 2] = rev_cnt[rev_cnt_endIdx];
			rev_cnt_endIdx = rev_cnt_endIdx - 2;
			isReadRev = 0;
		}
	}
	
	// out of data, count each remaining range as half cycles
	for (int i = 0; i < rev_cnt_endIdx; i++) {
		double range = rev_cnt[i + 1] - rev_cnt[i];
		if (range > 0) {
			DI = DI + 0.5 / 0.5 * pow(range / sigcr, 1 / m);
		}
	}

	// update variables
	sigPDI = sig_tmp;
	slopeP = slope;
	rev_cnt = rev_cnt_cache;
	rev_cnt_endIdx = rev_cnt_endIdx_cache;
	isReadRev = isReadRev_cache;
}

int
SteelFractureDI::returnSign(double v) {
	if (v < 0) return -1;
	if (v > 0) return 1;
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

double
SteelFractureDI::getDI(void)
{
	return DI;
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
	epsContP = epsCont;
	konfP = konf;
	konCP = konC;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = DI;
	DI_cacheP = DI_cache;
	isStartP = isStart;
	isSecondP = isSecond;
	sigPDIP = sigPDI;
	slopePP = slopeP;
	rev_cntP = rev_cnt;
	if (rev_cnt_endIdx + 2 >= rev_cnt.size()) {
		rev_cnt.resize(rev_cnt.size() + 10);
		rev_cntP.resize(rev_cnt.size() + 10);
	}
	rev_cnt_endIdxP = rev_cnt_endIdx;
	isReadRevP = isReadRev;
	// ************** added for DI ********************
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
	epsCont = epsContP;
	konf = konfP;
	konC = konCP;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DI = DIP;
	DI_cache = DI_cacheP;
	isStart = isStartP;
	isSecond = isSecondP;
	sigPDI = sigPDIP;
	slopeP = slopePP;
	rev_cnt = rev_cntP;
	rev_cnt_endIdx = rev_cnt_endIdxP;
	isReadRev = isReadRevP;
	// ************** added for DI ********************

	return 0;
}

int
SteelFractureDI::revertToStart(void)
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
	epsCont = 0.0;
	konf = 0;
	konC = 0;

	epsContP = 0.0;
	konfP = 0;
	konCP = 0;
	// ************** added for fracture ***************
	// ************** added for DI ********************
	DIP = 0.0;
	DI_cacheP = 0.0;
	isStartP = 1;
	isSecondP = 1;
	sigPDIP = 0.0;
	slopePP = 0.0;
	rev_cntP = {};
	rev_cntP.resize(10);
	rev_cnt_endIdxP = -1;
	isReadRevP = 0;

	DI = 0.0;
	DI_cache = 0.0;
	isStart = 1;
	isSecond = 1;
	sigPDI = 0.0;
	slopeP = 0.0;
	rev_cnt = {};
	rev_cnt.resize(10);
	rev_cnt_endIdx = -1;
	isReadRev = 0;
	// ************** added for DI ********************

	return 0;
}

int
SteelFractureDI::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
SteelFractureDI::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{

	return -1;
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
		s << "  sigcr: " << sigcr << ", ";
		s << "  m: " << m << ", ";
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
		s << "  sigcr: " << sigcr << ", ";
		s << "  m: " << m << ", ";
		//s << "\"sigini\": " << sigini << "}";
	}
}

Response*
SteelFractureDI::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
	if (argc == 0)
		return 0;

	Response *theResponse = 0;

	theOutput.tag("UniaxialMaterialOutput");
	theOutput.attr("matType", this->getClassType());
	theOutput.attr("matTag", this->getTag());


	// stress
	if (strcmp(argv[0], "stress") == 0) {
		theOutput.tag("ResponseType", "sigma11");
		theResponse = new MaterialResponse(this, 1, this->getStress());
	}
	// tangent
	else if (strcmp(argv[0], "tangent") == 0) {
		theOutput.tag("ResponseType", "C11");
		theResponse = new MaterialResponse(this, 2, this->getTangent());
	}

	// strain
	else if (strcmp(argv[0], "strain") == 0) {
		theOutput.tag("ResponseType", "eps11");
		theResponse = new MaterialResponse(this, 3, this->getStrain());
	}

	// strain
	else if ((strcmp(argv[0], "stressStrain") == 0) ||
		(strcmp(argv[0], "stressANDstrain") == 0)) {
		theOutput.tag("ResponseType", "sig11");
		theOutput.tag("ResponseType", "eps11");
		theResponse = new MaterialResponse(this, 4, Vector(2));
	}

	// damage 
	else if (strcmp(argv[0], "damage") == 0) {
		theResponse = new MaterialResponse(this, 5, this->getDI());
		theOutput.tag("ResponseType", "DI");
		// added 6/9/2006
	}

	else if (strcmp(argv[0], "failure") == 0) {
		int res = 0;
		theResponse = new MaterialResponse(this, 6, res);
		theOutput.tag("ResponseType", "Failure");
	}
	// end add


	theOutput.endTag();
	return theResponse;
}

int
SteelFractureDI::getResponse(int responseID, Information &matInfo)
{
	static Vector stressStrain(2);
	static Vector cyclesAndRange(6);

	// each subclass must implement its own stuff    
	switch (responseID) {
	case 1:
		matInfo.setDouble(this->getStress());
		return 0;

	case 2:
		matInfo.setDouble(this->getTangent());
		return 0;

	case 3:
		matInfo.setDouble(this->getStrain());
		return 0;

	case 4:
		stressStrain(0) = this->getStress();
		stressStrain(1) = this->getStrain();
		matInfo.setVector(stressStrain);
		return 0;

	case 5:
		matInfo.setDouble(this->getDI());
		return 0;

	case 6:
		if (DI > 1)
			matInfo.setInt(1);
		else
			matInfo.setInt(0);
		return 0;

	default:
		return -1;

	}
}
