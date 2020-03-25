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

// $Revision: 1.0 $
// $Date: 10/2018 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelFractureDI.h,v $

// Written: wyy
// Created: 05/2019
//


#ifndef SteelFractureDI_h
#define SteelFractureDI_h

#include <UniaxialMaterial.h>
// #include <vector>
// #include <deque>

class SteelFractureDI : public UniaxialMaterial
{
public:
	SteelFractureDI(int tag,
		double fy, double E0, double b,
		double R0, double cR1, double cR2,
		double a1, double a2, double a3, double a4, double sigcr, double m);

	SteelFractureDI(void);
	virtual ~SteelFractureDI();


	const char *getClassType(void) const { return "SteelFractureDI"; };

	double getInitialTangent(void);
	UniaxialMaterial *getCopy(void);

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

	// counting function
	void calcDI(double sigcr, double m, int &isStart, int &isSecond, double sig, double &sigPDI, double &DI_cache, double &DI, double &slopeP, double *rev_cnt, int &rev_cnt_endIdx, int &isReadRev);
	int sgn(double v);

	// override get-response function
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &matInformation);


protected:

private:
	// matpar : STEEL FIXED PROPERTIES
	double Fy;  //  = matpar(1)  : yield stress
	double E0;  //  = matpar(2)  : initial stiffness
	double b;   //  = matpar(3)  : hardening ratio (Esh/E0)
	double R0;  //  = matpar(4)  : exp transition elastic-plastic
	double cR1; //  = matpar(5)  : coefficient for changing R0 to R
	double cR2; //  = matpar(6)  : coefficient for changing R0 to R
	double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
	double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
	double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
	double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
	//double sigini; // initial
				   // hstvP : STEEL HISTORY VARIABLES
	double epsminP; //  = hstvP(1) : max eps in compression
	double epsmaxP; //  = hstvP(2) : max eps in tension
	double epsplP;  //  = hstvP(3) : plastic excursion
	double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
	double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
	double epssrP;  //  = hstvP(6) : eps at last inversion point
	double sigsrP;  //  = hstvP(7) : sig at last inversion point
	int    konP;    //  = hstvP(8) : index for loading/unloading
					// hstv : STEEL HISTORY VARIABLES   
	double epsP;  //  = strain at previous converged step
	double sigP;  //  = stress at previous converged step
	double eP;    //   stiffness modulus at last converged step;



	double epsmin;
	double epsmax;
	double epspl;
	double epss0;
	double sigs0;
	double epsr;
	double sigr;
	int    kon;
	double sig;
	double e;
	double eps;   //  = strain at current step

	// ************** added for fracture ***************
	// double epsPlCr;
	
	// double epsElNP;
	// double epsElPP;
	// double epsCumPlP;
	double epsContP;
	int konfP;
	int konCP;
	// double sigCumP;

	// double epsElN;
	// double epsElP;
	// double epsCumPl;
	double epsCont;
	int konf;
	int konC;
	// double sigCum;
	
	// ************** added for fracture ***************

	// ************** added for DI ***************
	double sigcr;
	double m;

	double DI;
	double DI_cache;
	int isStart;
	int isSecond;
	double sigPDI;
	double slopeP; // slope of previous stress history point
	// variables in counting calculations
	double rev_cnt[20] = { 0.0 };
	int rev_cnt_endIdx;
	int isReadRev;

	double DIP;
	double DI_cacheP;
	int isStartP;
	int isSecondP;
	double sigPDIP;
	double slopePP; // slope of previous stress history point
	// variables in counting calculations
	double rev_cntP[20] = { 0.0 };
	int rev_cnt_endIdxP;
	int isReadRevP;
	
	// ************** added for DI ***************

	// ***** cycle counting (ASTM one-pass algorithm) *****
	
	// ***** private function: cycle counting (ASTM 3-point algorithm) *****
};


#endif

