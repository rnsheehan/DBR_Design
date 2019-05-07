#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in the DBR class
// R. Sheehan 7 - 5 - 2019

DBR::DBR()
{
	// Default constructor
	
}

DBR::~DBR()
{
	// Deconstructor
	neff_m.clear(); 
	neff_l.clear(); 
	neff_h.clear(); 
}