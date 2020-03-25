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

// Written: Chris McGann, U.Washington
//          02.2011
//
// Description: This file contains the class definition for SurfaceLoader.

#ifndef SurfaceLoader_h
#define SurfaceLoader_h

#include <ElementalLoad.h>

class SurfaceLoader : public ElementalLoad
{
	public:
		SurfaceLoader(int tag, int eleTag);
		SurfaceLoader();
		~SurfaceLoader();

		const Vector &getData(int &type, double loadFactor);

		int sendSelf(int commitTag, Channel &theChannel);
		int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
		void Print(OPS_Stream &s, int flag = 0);

	protected:

	private:
		static Vector data;
};

#endif
